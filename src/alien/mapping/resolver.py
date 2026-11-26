"""The mapping decision for one (symbol, source Ensembl ID, NCBI-allowed) key in one target namespace.

:func:`resolve` is pure given its context. It runs an ordered chain of small steps; each returns a final
:class:`MappingOutcome` or ``None`` to pass control on. Steps after symbol resolution share a
:class:`_State` so that the NCBI rescue step appears once, whichever earlier step failed.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from typing import Any, Final

from ..model import AnnotationCandidate, MappingStatus, ensembl_stable_id
from ..resources.ensembl_archive import archive_target_matches, unresolved_status
from ..resources.hgnc import HgncMaps
from ..resources.ncbi import ID_STAGE as NCBI_ID_STAGE
from ..resources.ncbi import RESCUE_STATUSES as NCBI_RESCUE_STATUSES
from ..resources.ncbi import STAGES as NCBI_STAGES
from ..resources.ncbi import NcbiMaps, NcbiRecord
from .symbols import AMBIGUOUS_STAGES, SHORT_ALIAS, SYMBOL_STATUSES, SymbolOptions, SymbolResolution, resolve_symbol

MappingKey = tuple[str, str, bool]  # (gene symbol, stable source Ensembl ID, NCBI rescue allowed)
CANDIDATE_PRIORITY: Final = {"annotation_gtf": 0, "metadata_fallback": 1, "output_genes": 2, "HGNC": 3}
PROTEIN_CODING: Final = "protein_coding"
_LOC_GENE_ID: Final = re.compile(r"^LOC([0-9]+)$")


@dataclass(frozen=True, slots=True)
class NcbiRescue:
    status: str
    stage: str = ""
    current_symbol: str = ""
    current_gene_ids: tuple[str, ...] = ()
    ncbi_ensembl_ids: tuple[str, ...] = ()
    possible_symbols: tuple[str, ...] = ()
    possible_ensembl_ids: tuple[str, ...] = ()
    ensembl_gene_id: str = ""

    @property
    def rescued(self) -> bool:
        return self.status in NCBI_RESCUE_STATUSES

    @property
    def checked(self) -> bool:
        return self.status != "ncbi_not_checked"


NOT_CHECKED: Final = NcbiRescue("ncbi_not_checked")


@dataclass(frozen=True, slots=True)
class AmbiguousDetail:
    stage: str
    symbols: tuple[str, ...]
    ensembl_ids: tuple[str, ...]
    action: str


@dataclass(frozen=True, slots=True)
class UnmappedDetail:
    attempted_symbol: str
    reason: str
    matches: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class MappingOutcome:
    status: MappingStatus
    ensembl_id: str = ""
    current_symbol: str = ""
    candidates: tuple[str, ...] = ()
    failure: str = ""  # archive-audit label for a failed outcome
    ambiguous: AmbiguousDetail | None = None
    unmapped: UnmappedDetail | None = None
    ncbi: tuple[NcbiRescue, str] | None = None  # (rescue attempt, audit action)
    archive: tuple[dict[str, Any], str] | None = None  # (archive record, rescue_status)


@dataclass(frozen=True)
class MappingContext:
    """Everything one target namespace needs to resolve keys."""

    target_ids: frozenset[str]
    known_ids: frozenset[str]
    symbol_index: Mapping[str, Sequence[AnnotationCandidate]]
    hgnc: HgncMaps
    ncbi: NcbiMaps
    archive_lookup: Callable[[str], dict[str, Any]] | None
    repairs: Mapping[str, str]
    use_alias: bool
    is_non_gene: Callable[[object], bool]
    case_insensitive: bool = False
    symbol_options: SymbolOptions | None = None

    @property
    def options(self) -> SymbolOptions:
        if self.symbol_options is not None:
            return self.symbol_options
        return SymbolOptions(use_alias=self.use_alias, case_insensitive=self.case_insensitive)


@dataclass(slots=True)
class _State:
    key: MappingKey
    resolution: SymbolResolution | None = None
    candidates: list[AnnotationCandidate] = field(default_factory=list)
    pending: str = ""  # failure awaiting an NCBI rescue: unmapped_symbol, not_in_output_genes, ambiguous_gencode
    rescue: NcbiRescue = NOT_CHECKED


# ------------------------------------------------------------------------------------------ candidates


def choose_candidate(candidates: Sequence[AnnotationCandidate]) -> AnnotationCandidate | None:
    """One candidate, else prefer output genes, then annotation source priority, then a unique protein-coding gene."""
    if len(candidates) <= 1:
        return candidates[0] if candidates else None
    pool = _best_sourced(candidates)
    unique = pool if len(pool) == 1 else [c for c in pool if c.gene_biotype == PROTEIN_CODING]
    return unique[0] if len(unique) == 1 else None


def _best_sourced(candidates: Sequence[AnnotationCandidate]) -> list[AnnotationCandidate]:
    """Candidates in the output genes (if any), restricted to the best annotation source."""
    pool = [c for c in candidates if c.is_in_output_genes] or list(candidates)
    best = min(_source_rank(c) for c in pool)
    return [c for c in pool if _source_rank(c) == best]


def _source_rank(candidate: AnnotationCandidate) -> int:
    return CANDIDATE_PRIORITY.get(candidate.annotation_source, len(CANDIDATE_PRIORITY))


def hgnc_tiebreak(
    symbol: str, candidates: Sequence[AnnotationCandidate], ctx: MappingContext
) -> AnnotationCandidate | None:
    """The official HGNC Ensembl ID, when it is a target ID among several same-named candidates."""
    official = ensembl_stable_id(ctx.hgnc.current_to_ensembl.get(symbol, ""))
    if not official or official not in ctx.target_ids:
        return None
    return next((c for c in candidates if c.ensembl_gene_id == official), None)


def symbol_candidates(symbol: str, ctx: MappingContext) -> list[AnnotationCandidate]:
    candidates = list(ctx.symbol_index.get(symbol, ()))
    official = ensembl_stable_id(ctx.hgnc.current_to_ensembl.get(symbol, ""))
    if official and official not in {c.ensembl_gene_id for c in candidates}:
        candidates.append(AnnotationCandidate(official, "", True, "HGNC"))
    return candidates


def _in_target(candidates: Sequence[AnnotationCandidate], ctx: MappingContext) -> list[AnnotationCandidate]:
    return [c for c in candidates if c.ensembl_gene_id in ctx.target_ids]


# ------------------------------------------------------------------------------------------ NCBI rescue


def _ncbi_candidate(
    record: NcbiRecord, ctx: MappingContext
) -> tuple[AnnotationCandidate | None, list[AnnotationCandidate]]:
    pool = [AnnotationCandidate(i, "", True, "NCBI") for i in record.ensembl_gene_ids if i in ctx.target_ids]
    pool += _in_target(symbol_candidates(record.current_symbol, ctx), ctx)
    first: dict[str, AnnotationCandidate] = {}
    for candidate in pool:
        first.setdefault(candidate.ensembl_gene_id, candidate)
    unique = list(first.values())
    chosen = choose_candidate(unique) or hgnc_tiebreak(record.current_symbol, unique, ctx)
    if chosen is None and len(unique) == 1:
        chosen = unique[0]
    return chosen, unique


def ncbi_gene_id_rescue(symbol_raw: str, ctx: MappingContext) -> NcbiRescue:
    """``LOC<GeneID>``: the gene that NCBI GeneID is now (following one replacement), mapped to a target ID."""
    match = _LOC_GENE_ID.match(symbol_raw.strip())
    if match is None or not ctx.ncbi:
        return NOT_CHECKED
    record, state = ctx.ncbi.by_gene_id(match.group(1))
    if record is None:
        return NcbiRescue(f"{NCBI_ID_STAGE}_{state}", NCBI_ID_STAGE)
    chosen, candidates = _ncbi_candidate(record, ctx)
    base = NcbiRescue(
        NCBI_ID_STAGE,
        NCBI_ID_STAGE,
        record.current_symbol,
        record.current_gene_ids,
        record.ensembl_gene_ids,
        (record.current_symbol,),
        tuple(sorted({c.ensembl_gene_id for c in candidates})),
    )
    if chosen is not None:
        return replace(base, ensembl_gene_id=chosen.ensembl_gene_id)
    return replace(base, status=f"ambiguous_{NCBI_ID_STAGE}" if candidates else f"{NCBI_ID_STAGE}_not_in_output_genes")


def ncbi_rescue(symbol_raw: str, allowed: bool, ctx: MappingContext) -> NcbiRescue:
    """NCBI Gene history, then gene_info synonyms; only after HGNC failed and only for allowed sources."""
    symbol = symbol_raw.strip()
    if not symbol or not ctx.ncbi or not allowed:
        return NOT_CHECKED
    for stage in NCBI_STAGES:
        if stage == "ncbi_gene_info_synonym" and len(symbol) < ctx.options.min_alias_length:
            continue
        if symbol in ctx.ncbi.ambiguous.get(stage, {}):
            return NcbiRescue(f"ambiguous_{stage}", stage, possible_symbols=tuple(ctx.ncbi.ambiguous[stage][symbol]))
        record = ctx.ncbi.stages.get(stage, {}).get(symbol)
        if record is None:
            continue
        chosen, candidates = _ncbi_candidate(record, ctx)
        base = NcbiRescue(
            stage,
            stage,
            record.current_symbol,
            record.current_gene_ids,
            record.ensembl_gene_ids,
            (record.current_symbol,),
            tuple(sorted({c.ensembl_gene_id for c in candidates})),
        )
        if chosen is not None:
            return replace(base, ensembl_gene_id=chosen.ensembl_gene_id)
        return replace(base, status=f"ambiguous_{stage}" if candidates else f"{stage}_not_in_output_genes")
    return NOT_CHECKED


# ------------------------------------------------------------------------------------------ steps

Step = Callable[[_State, MappingContext], MappingOutcome | None]


def _symbol_resolution(state: _State, ctx: MappingContext) -> MappingOutcome | None:
    resolution = resolve_symbol(state.key[0], ctx.hgnc, ctx.repairs, options=ctx.options)
    state.resolution = resolution
    if resolution.is_ambiguous:
        detail = AmbiguousDetail(AMBIGUOUS_STAGES[resolution.status], resolution.candidates, (), "drop")
        return MappingOutcome(MappingStatus.AMBIGUOUS_DROP, failure="ambiguous_symbol", ambiguous=detail)
    if not resolution.symbol:
        state.pending = "unmapped_symbol"
    return None


def _annotation_candidate(state: _State, ctx: MappingContext) -> MappingOutcome | None:
    resolution = state.resolution
    if state.pending or resolution is None:
        return None
    state.candidates = _in_target(symbol_candidates(resolution.symbol, ctx), ctx)
    if not state.candidates:
        state.pending = "not_in_output_genes"
        return None
    chosen = choose_candidate(state.candidates)
    if chosen is None:
        state.pending = "ambiguous_gencode"
        return None
    status = MappingStatus(SYMBOL_STATUSES[resolution.status])
    return MappingOutcome(status, chosen.ensembl_gene_id, resolution.symbol)


def _hgnc_tiebreak(state: _State, ctx: MappingContext) -> MappingOutcome | None:
    if state.pending != "ambiguous_gencode" or state.resolution is None:
        return None
    chosen = hgnc_tiebreak(state.resolution.symbol, state.candidates, ctx)
    if chosen is None:
        return None
    return MappingOutcome(MappingStatus.HGNC_ENSEMBL_TIEBREAK, chosen.ensembl_gene_id, state.resolution.symbol)


def _ncbi_gene_id(state: _State, ctx: MappingContext) -> MappingOutcome | None:
    """``LOC<GeneID>`` tokens: an ID lookup, so it runs for every source, not only NCBI-allowed ones."""
    if not state.pending:
        return None
    rescue = ncbi_gene_id_rescue(state.key[0], ctx)
    if rescue.checked:
        state.rescue = rescue
    if not rescue.rescued:
        return None
    return MappingOutcome(
        MappingStatus.NCBI_GENE_ID, rescue.ensembl_gene_id, rescue.current_symbol, ncbi=(rescue, "rescued")
    )


def _ncbi_rescue(state: _State, ctx: MappingContext) -> MappingOutcome | None:
    if not state.pending:
        return None
    rescue = ncbi_rescue(state.key[0], state.key[2], ctx)
    if rescue.checked or not state.rescue.checked:
        state.rescue = rescue
    if not state.rescue.rescued:
        return None
    return MappingOutcome(
        MappingStatus(state.rescue.status),
        state.rescue.ensembl_gene_id,
        state.rescue.current_symbol,
        ncbi=(state.rescue, "rescued"),
    )


def _failure(state: _State, ctx: MappingContext) -> MappingOutcome:
    rescue, resolution = state.rescue, state.resolution
    current = resolution.symbol if resolution else ""
    if state.pending == "ambiguous_gencode":
        detail = AmbiguousDetail(
            "gencode_gene_name", (current,), tuple(c.ensembl_gene_id for c in state.candidates), "drop"
        )
        return MappingOutcome(MappingStatus.AMBIGUOUS_DROP, failure=state.pending, ambiguous=detail)
    if state.pending == "not_in_output_genes":
        return MappingOutcome(
            MappingStatus.NOT_IN_OUTPUT_GENES,
            current_symbol=current,
            failure=state.pending,
            unmapped=UnmappedDetail(current, "not_in_output_genes", ()),
            ncbi=(rescue, "not_in_output_genes") if rescue.checked else None,
        )
    if rescue.status.startswith("ambiguous"):
        detail = AmbiguousDetail(rescue.status, rescue.possible_symbols, rescue.possible_ensembl_ids, "drop")
        return MappingOutcome(
            MappingStatus.AMBIGUOUS_DROP, failure=state.pending, ambiguous=detail, ncbi=(rescue, "ambiguous_drop")
        )
    matches = resolution.candidates if resolution else ()
    reason = SHORT_ALIAS if resolution is not None and resolution.status == SHORT_ALIAS else "unmapped_symbol"
    return MappingOutcome(
        MappingStatus.UNMAPPED_DROP,
        failure=state.pending,
        unmapped=UnmappedDetail("", reason, matches),
        ncbi=(rescue, "unmapped_drop") if rescue.checked else None,
    )


SYMBOL_STEPS: Final[tuple[Step, ...]] = (
    _symbol_resolution,
    _annotation_candidate,
    _hgnc_tiebreak,
    _ncbi_gene_id,
    _ncbi_rescue,
)


def _archive_check(key: MappingKey, ctx: MappingContext) -> tuple[dict[str, Any], list[str]] | None:
    source = key[1]
    if not source or ctx.archive_lookup is None or source in ctx.known_ids:
        return None
    info = ctx.archive_lookup(source)
    return info, archive_target_matches(info, ctx.target_ids)


def resolve(key: MappingKey, ctx: MappingContext) -> MappingOutcome:
    """Non-gene token, direct source Ensembl, archive rescue, then the symbol chain and failure classification."""
    symbol, source, _ = key
    if ctx.is_non_gene(symbol):
        return MappingOutcome(MappingStatus.NON_GENE_TOKEN_DROP)
    if source and source in ctx.target_ids:
        return MappingOutcome(MappingStatus.DIRECT_SOURCE_ENSEMBL, source)
    archive = _archive_check(key, ctx)
    if archive is not None and len(archive[1]) == 1:
        rescued = archive[1][0]
        return MappingOutcome(
            MappingStatus.ARCHIVE_ENSEMBL_RESCUE, rescued, archive=(archive[0], "archive_ensembl_rescue")
        )
    state = _State(key)
    outcome = next((result for step in SYMBOL_STEPS if (result := step(state, ctx)) is not None), None)
    outcome = outcome or _failure(state, ctx)
    if archive is None:
        return outcome
    info, matches = archive
    label = outcome.status.value if outcome.status.is_mapped else unresolved_status(info, matches, outcome.failure)
    return replace(outcome, archive=(info, label))
