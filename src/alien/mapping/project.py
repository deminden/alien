"""Project canonical memberships onto the symbols namespace and every target namespace.

Each distinct (symbol, source Ensembl ID, NCBI-allowed) key is resolved once per target; outcomes are
then joined back to the membership rows to build term gene sets, audit tables and status counts.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np
import pandas as pd

from ..config.schema import AlienConfig
from ..io.gmt import description_from_meta
from ..model import (
    SYMBOLS_NAMESPACE,
    AnnotationCandidate,
    MappingStatus,
    Target,
    TermRecord,
    TermSet,
    ensembl_stable_id,
)
from ..resources.ensembl_archive import EnsemblArchiveResolver
from ..resources.hgnc import HgncMaps
from ..resources.ncbi import NcbiMaps
from ..schema import (
    AMBIGUOUS_COLUMNS,
    ARCHIVE_AUDIT_COLUMNS,
    MAPPING_STATUS_COLUMNS,
    NCBI_AUDIT_COLUMNS,
    NON_GENE_COLUMNS,
    SYMBOL_RESCUE_COLUMNS,
    TERM_META_COLUMNS,
    UNMAPPED_COLUMNS,
)
from .audit import (
    Fields,
    ambiguous_fields,
    archive_fields,
    audit_table,
    ncbi_fields,
    non_gene_fields,
    unmapped_fields,
)
from .resolver import MappingContext, MappingOutcome, resolve
from .symbols import (
    AMBIGUOUS_STAGES,
    SHORT_ALIAS,
    SymbolOptions,
    non_gene_matcher,
    resolve_symbol,
    short_alias_guarded,
)

LOGGER = logging.getLogger("alien")
AUDIT_TABLES = ("unmapped", "ambiguous", "non_gene", "archive", "ncbi")
#: Statuses that are not a plain exact match; each distinct token reaching one is listed in symbol_rescues.
_PLAIN_STATUSES = frozenset({MappingStatus.DIRECT_SOURCE_ENSEMBL, MappingStatus.HGNC_CURRENT_SYMBOL})
AUDIT_COLUMNS = {
    "unmapped": UNMAPPED_COLUMNS,
    "ambiguous": AMBIGUOUS_COLUMNS,
    "non_gene": NON_GENE_COLUMNS,
    "archive": ARCHIVE_AUDIT_COLUMNS,
    "ncbi": NCBI_AUDIT_COLUMNS,
}


@dataclass
class MappingResult:
    """Per-namespace term sets plus the audit tables and per-status membership counts."""

    terms: dict[str, TermSet] = field(default_factory=dict)
    unmapped: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=list(UNMAPPED_COLUMNS)))
    ambiguous: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=list(AMBIGUOUS_COLUMNS)))
    non_gene: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=list(NON_GENE_COLUMNS)))
    archive_audit: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=list(ARCHIVE_AUDIT_COLUMNS)))
    ncbi_audit: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=list(NCBI_AUDIT_COLUMNS)))
    symbol_rescues: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=list(SYMBOL_RESCUE_COLUMNS)))
    status_summary: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=list(MAPPING_STATUS_COLUMNS)))


@dataclass
class _Projection:
    terms: TermSet
    tables: dict[str, pd.DataFrame]
    status: pd.DataFrame | None = None
    rescues: pd.DataFrame | None = None


class Memberships:
    """Canonical membership rows with the derived arrays every namespace needs."""

    def __init__(self, frame: pd.DataFrame) -> None:
        columns = list(dict.fromkeys([*TERM_META_COLUMNS, "gene_symbol", "gene_ensembl_from_source"]))
        self.rows = pd.DataFrame({c: frame[c].fillna("").astype(str) if c in frame else "" for c in columns})
        self.rows = self.rows.reset_index(drop=True)
        unique_sources = pd.unique(self.rows["gene_ensembl_from_source"])
        stable = {value: ensembl_stable_id(value) for value in unique_sources}
        self.source_stable = self.rows["gene_ensembl_from_source"].map(stable).to_numpy(dtype=object)
        self.symbols = self.rows["gene_symbol"].to_numpy(dtype=object)
        self.term_ids = self.rows["term_id"].to_numpy(dtype=object)
        self.meta = {column: self.rows[column].to_numpy(dtype=object) for column in TERM_META_COLUMNS}

    def __len__(self) -> int:
        return len(self.rows)

    def collect_terms(self, genes: np.ndarray, mapped: np.ndarray) -> TermSet:
        """Group mapped rows into term records; a term's metadata comes from its first mapped row."""
        terms: TermSet = {}
        for index in np.flatnonzero(mapped):
            gene = genes[index]
            if not gene:
                continue
            term_id = self.term_ids[index]
            record = terms.get(term_id)
            if record is None:
                meta = {column: values[index] for column, values in self.meta.items()}
                record = terms[term_id] = TermRecord(term_id, description_from_meta(meta), meta)
            record.genes.add(str(gene))
            record.original_genes.add(str(self.symbols[index]))
        return terms


def _factorize(*arrays: np.ndarray) -> tuple[np.ndarray, list[tuple[Any, ...]]]:
    codes, uniques = pd.MultiIndex.from_arrays(list(arrays)).factorize()
    return np.asarray(codes), list(uniques)


def project_symbols(members: Memberships, config: AlienConfig, hgnc: HgncMaps) -> _Projection:
    """The symbols namespace: current HGNC symbols, keeping unresolved tokens when configured."""
    mapping = config.gene_mapping
    options = {guard: SymbolOptions.from_config(mapping, short_alias_guard=guard) for guard in (False, True)}
    is_non_gene = non_gene_matcher(mapping)
    codes, uniques = _factorize(members.symbols, _guards(members, config))
    outputs: list[str] = []
    per_table: dict[str, list[Fields | None]] = {"non_gene": [], "ambiguous": [], "unmapped": []}
    for raw, guard in uniques:
        symbol = str(raw)
        fields: dict[str, Fields | None] = dict.fromkeys(per_table)
        out = ""
        if is_non_gene(symbol):
            fields["non_gene"] = {"reason": "configured_non_gene_token"}
        else:
            resolution = resolve_symbol(symbol, hgnc, mapping.manual_symbol_repairs, options=options[bool(guard)])
            out = resolution.symbol
            if resolution.is_ambiguous:
                fields["ambiguous"] = {
                    "mapping_stage": AMBIGUOUS_STAGES[resolution.status],
                    "possible_symbols": ",".join(sorted(set(resolution.candidates))),
                    "possible_ensembl_ids": "",
                    "action": "drop_alias",
                }
            elif not resolution.symbol:
                short = resolution.status == SHORT_ALIAS
                fields["unmapped"] = {
                    "attempted_current_symbol": "",
                    "reason": SHORT_ALIAS if short else "unmapped_symbol",
                    "possible_matches": ",".join(resolution.candidates) if short else "",
                }
            if not out and mapping.preserve_unmapped_symbols_in_symbol_gmt:
                out = symbol
        outputs.append(out)
        for name, value in fields.items():
            per_table[name].append(value)
    codes = np.asarray(codes)
    genes = np.asarray(outputs, dtype=object)[codes] if outputs else np.asarray([], dtype=object)
    terms = members.collect_terms(genes, genes != "")
    tables = {
        name: audit_table(members.rows, SYMBOLS_NAMESPACE, codes, per_table[name], AUDIT_COLUMNS[name])
        for name in per_table
    }
    return _Projection(terms, tables)


def symbol_index(annotation: pd.DataFrame) -> dict[str, list[AnnotationCandidate]]:
    index: dict[str, list[AnnotationCandidate]] = {}
    for gene_id, symbol, biotype, in_output, source in zip(
        annotation["ensembl_gene_id"],
        annotation["gene_symbol"],
        annotation["gene_biotype"],
        annotation["is_in_output_genes"],
        annotation["target_id_metadata_source"],
        strict=True,
    ):
        if not symbol or pd.isna(symbol):
            continue
        flag = in_output is True or (isinstance(in_output, np.bool_) and bool(in_output)) or str(in_output) == "True"
        index.setdefault(str(symbol), []).append(
            AnnotationCandidate(
                str(gene_id), "" if pd.isna(biotype) else str(biotype), flag, str(source or "annotation_gtf")
            )
        )
    return index


def target_context(
    target: Target, config: AlienConfig, hgnc: HgncMaps, ncbi: NcbiMaps, resolver: EnsemblArchiveResolver | None
) -> MappingContext:
    target_ids = target.effective_ids
    known = target.annotation_ids | target_ids | target.output_gene_ids | target.filter_ids | hgnc.current_ensembl_ids
    return MappingContext(
        target_ids=target_ids,
        known_ids=known,
        symbol_index=symbol_index(target.annotation),
        hgnc=hgnc,
        ncbi=ncbi,
        archive_lookup=resolver.lookup if resolver is not None else None,
        repairs=config.gene_mapping.manual_symbol_repairs,
        use_alias=config.gene_mapping.use_hgnc_alias_mapping,
        is_non_gene=non_gene_matcher(config.gene_mapping),
        case_insensitive=config.gene_mapping.case_insensitive_symbols,
        symbol_options=SymbolOptions.from_config(config.gene_mapping),
    )


def problematic_source_ids(members: Memberships, known_ids: frozenset[str]) -> list[str]:
    """Source Ensembl IDs unknown to the annotation, target lists and HGNC: the only ones sent to the archive."""
    return sorted({stable for stable in set(members.source_stable) if stable and stable not in known_ids})


def _guards(members: Memberships, config: AlienConfig) -> np.ndarray:
    tags = members.rows["source_tag"]
    unique = list(pd.unique(tags))
    guarded = dict(zip(unique, short_alias_guarded(config.gene_mapping, unique), strict=True))
    return tags.map(guarded).to_numpy(dtype=bool)


def project_target(
    members: Memberships, target: Target, ctx: MappingContext, allowed: np.ndarray, guards: np.ndarray | None = None
) -> _Projection:
    """``guards`` marks rows whose source is subject to the short-alias rule (default: the context's own options)."""
    if guards is None:
        guards = np.ones(len(members), dtype=bool)
    open_ctx = replace(ctx, symbol_options=replace(ctx.options, min_alias_length=1))
    codes, keys = _factorize(members.symbols, members.source_stable, allowed, guards)
    outcomes: list[MappingOutcome] = [resolve((str(s), str(e), bool(a)), ctx if g else open_ctx) for s, e, a, g in keys]
    statuses = np.asarray([outcome.status.value for outcome in outcomes], dtype=object)[codes]
    genes = np.asarray([outcome.ensembl_id for outcome in outcomes], dtype=object)[codes]
    mapped = np.asarray([outcome.status.is_mapped for outcome in outcomes], dtype=bool)[codes]
    terms = members.collect_terms(genes, mapped)
    per_code: dict[str, list[Fields | None]] = {
        "unmapped": [unmapped_fields(o) for o in outcomes],
        "ambiguous": [ambiguous_fields(o) for o in outcomes],
        "non_gene": [non_gene_fields(o) for o in outcomes],
        "archive": [archive_fields(o, str(key[1])) for o, key in zip(outcomes, keys, strict=True)],
        "ncbi": [ncbi_fields(o) for o in outcomes],
    }
    tables = {
        name: audit_table(members.rows, target.name, codes, per_code[name], AUDIT_COLUMNS[name])
        for name in AUDIT_TABLES
    }
    counts = (
        pd.DataFrame(
            {"family": members.rows["family"], "source_tag": members.rows["source_tag"], "mapping_status": statuses}
        )
        .groupby(["family", "source_tag", "mapping_status"], sort=True)
        .size()
        .reset_index(name="n_gene_members")
    )
    counts.insert(0, "target_namespace", target.name)
    return _Projection(terms, tables, counts, rescue_table(members, target.name, codes, keys, outcomes))


def rescue_table(
    members: Memberships,
    namespace: str,
    codes: np.ndarray,
    keys: Sequence[tuple[Any, ...]],
    outcomes: Sequence[MappingOutcome],
) -> pd.DataFrame:
    """One row per (source tag, token) mapped by anything other than a direct ID or exact current symbol."""
    flagged = np.asarray([o.status.is_mapped and o.status not in _PLAIN_STATUSES for o in outcomes], dtype=bool)
    selected = np.flatnonzero(flagged[codes]) if len(codes) else np.asarray([], dtype=int)
    if not len(selected):
        return pd.DataFrame(columns=list(SYMBOL_RESCUE_COLUMNS))
    frame = pd.DataFrame(
        {
            "code": codes[selected],
            "source_tag": members.rows["source_tag"].to_numpy()[selected],
            "term_id": members.term_ids[selected],
        }
    )
    grouped = (
        frame.groupby(["source_tag", "code"], sort=False)
        .agg(n_memberships=("term_id", "size"), n_terms=("term_id", "nunique"), example_term_id=("term_id", "min"))
        .reset_index()
    )
    rows = []
    for tag, code, n_rows, n_terms, example in grouped.itertuples(index=False):
        outcome, key = outcomes[code], keys[code]
        rows.append(
            {
                "target_namespace": namespace,
                "source_tag": tag,
                "input_gene": str(key[0]),
                "source_ensembl_id": str(key[1]),
                "mapping_status": outcome.status.value,
                "mapped_symbol": outcome.current_symbol,
                "ensembl_gene_id": outcome.ensembl_id,
                "n_memberships": int(n_rows),
                "n_terms": int(n_terms),
                "example_term_id": example,
            }
        )
    table = pd.DataFrame(rows, columns=list(SYMBOL_RESCUE_COLUMNS))
    return table.sort_values(
        ["target_namespace", "mapping_status", "source_tag", "n_memberships", "input_gene"],
        ascending=[True, True, True, False, True],
        kind="mergesort",
    ).reset_index(drop=True)


def _ncbi_allowed(members: Memberships, config: AlienConfig, ncbi: NcbiMaps) -> np.ndarray:
    if not ncbi:
        return np.zeros(len(members), dtype=bool)
    tags = config.ncbi_gene.allowed_source_tags
    if not tags:
        return np.ones(len(members), dtype=bool)
    return members.rows["source_tag"].isin(set(tags)).to_numpy()


def _concat(frames: Sequence[pd.DataFrame], columns: Sequence[str]) -> pd.DataFrame:
    kept = [frame for frame in frames if not frame.empty]
    return pd.concat(kept, ignore_index=True)[list(columns)] if kept else pd.DataFrame(columns=list(columns))


def map_namespaces(
    memberships: pd.DataFrame,
    config: AlienConfig,
    targets: Sequence[Target],
    *,
    hgnc: HgncMaps,
    ncbi: NcbiMaps,
    resolver: EnsemblArchiveResolver | None,
) -> MappingResult:
    members = Memberships(memberships)
    projections: list[_Projection] = []
    result = MappingResult()
    if config.outputs.include_symbols:
        projection = project_symbols(members, config, hgnc)
        result.terms[SYMBOLS_NAMESPACE] = projection.terms
        projections.append(projection)
    contexts = {target.name: target_context(target, config, hgnc, ncbi, resolver) for target in targets}
    if resolver is not None:
        for target in targets:
            resolver.prefetch(problematic_source_ids(members, contexts[target.name].known_ids), target.name)
    allowed = _ncbi_allowed(members, config, ncbi)
    guards = _guards(members, config)
    for target in targets:
        LOGGER.info("Projecting %d memberships onto %s", len(members), target.name)
        projection = project_target(members, target, contexts[target.name], allowed, guards)
        result.terms[target.name] = projection.terms
        projections.append(projection)
    tables = {
        name: _concat([p.tables[name] for p in projections if name in p.tables], AUDIT_COLUMNS[name])
        for name in AUDIT_TABLES
    }
    result.unmapped, result.ambiguous, result.non_gene = tables["unmapped"], tables["ambiguous"], tables["non_gene"]
    result.archive_audit, result.ncbi_audit = tables["archive"], tables["ncbi"]
    result.symbol_rescues = _concat([p.rescues for p in projections if p.rescues is not None], SYMBOL_RESCUE_COLUMNS)
    status = _concat([p.status for p in projections if p.status is not None], MAPPING_STATUS_COLUMNS)
    result.status_summary = status.sort_values(list(MAPPING_STATUS_COLUMNS[:4]), kind="mergesort").reset_index(
        drop=True
    )
    return result


def status_counts(frame: pd.DataFrame) -> dict[MappingStatus, int]:
    """Total memberships per status (handy for tests and BuildResult)."""
    return {MappingStatus(str(k)): int(v) for k, v in frame.groupby("mapping_status")["n_gene_members"].sum().items()}
