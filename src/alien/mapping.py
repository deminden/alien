from __future__ import annotations

import json
import re
import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from .utils import LOGGER, description_from_record, ensure_dirs, strip_ensembl_version, write_json


def map_namespaces(
    term_df: pd.DataFrame,
    hgnc_maps: dict[str, object],
    ncbi_maps: dict[str, object],
    targets: list[dict[str, Any]],
    cfg: dict,
    source_dir: Path | None = None,
    workers: int = 1,
) -> tuple[dict[str, dict[str, dict[str, Any]]], pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    # Map terms to current symbols and each configured target namespace.
    workers = max(1, int(workers or 1))
    unmapped: list[dict[str, object]] = []
    ambiguous: list[dict[str, object]] = []
    non_gene: list[dict[str, object]] = []
    archive_audit: list[dict[str, object]] = []
    ncbi_audit: list[dict[str, object]] = []
    mapping_status: dict[tuple[str, str, str, str], int] = {}
    archive_resolver = _build_archive_resolver(source_dir, cfg)

    mapped: dict[str, dict[str, dict[str, Any]]] = {}
    if cfg.get("outputs", {}).get("include_symbols", True):
        mapped["symbols"] = map_symbols(term_df, hgnc_maps, cfg, unmapped, ambiguous, non_gene)

    if workers > 1 and len(targets) > 1:
        LOGGER.info("Mapping %d target namespaces serially; workers are used for filtering in this release.", len(targets))

    for target in targets:
        namespace = str(target["name"])
        mapped[namespace] = map_ensembl(
            term_df,
            hgnc_maps,
            ncbi_maps,
            target["annotation"],
            target.get("universe"),
            namespace,
            cfg,
            unmapped,
            ambiguous,
            non_gene,
            archive_audit,
            ncbi_audit,
            mapping_status,
            archive_resolver,
        )

    if archive_resolver is not None:
        archive_resolver.save()
    return (
        mapped,
        pd.DataFrame(unmapped),
        pd.DataFrame(ambiguous),
        pd.DataFrame(non_gene),
        pd.DataFrame(archive_audit),
        pd.DataFrame(ncbi_audit),
        _mapping_status_frame(mapping_status),
    )


def map_symbols(
    term_df: pd.DataFrame,
    hgnc_maps: dict[str, object],
    cfg: dict,
    unmapped: list[dict[str, object]],
    ambiguous: list[dict[str, object]],
    non_gene: list[dict[str, object]] | None = None,
) -> dict[str, dict[str, Any]]:
    # Resolve source symbols to current HGNC symbols
    output: dict[str, dict[str, Any]] = {}
    non_gene = non_gene if non_gene is not None else []
    preserve = cfg.get("gene_mapping", {}).get("preserve_unmapped_symbols_in_symbol_gmt", True)
    for _, row in term_df.iterrows():
        if _is_non_gene_token(row.get("gene_symbol", ""), cfg):
            non_gene.append(_non_gene_row("symbols", row, row["gene_symbol"], "configured_non_gene_token"))
            continue
        symbol, status, matches = resolve_current_symbol(str(row["gene_symbol"]), hgnc_maps)
        if status == "ambiguous_alias":
            ambiguous.append(_ambiguous_row("symbols", row, row["gene_symbol"], "hgnc_alias_symbol", matches, [], "drop_alias"))
            if not preserve:
                continue
            symbol = str(row["gene_symbol"])
        elif status == "unmapped":
            unmapped.append(_unmapped_row("symbols", row, row["gene_symbol"], "", "unmapped_symbol", []))
            if not preserve:
                continue
            symbol = str(row["gene_symbol"])
        _add_gene(output, row, symbol)
    return output


def map_ensembl(
    term_df: pd.DataFrame,
    hgnc_maps: dict[str, object],
    ncbi_maps: dict[str, object] | None,
    gencode: pd.DataFrame,
    universe: set[str] | None,
    namespace: str,
    cfg: dict,
    unmapped: list[dict[str, object]],
    ambiguous: list[dict[str, object]],
    non_gene: list[dict[str, object]] | None = None,
    archive_audit: list[dict[str, object]] | None = None,
    ncbi_audit: list[dict[str, object]] | None = None,
    mapping_status: dict[tuple[str, str, str, str], int] | None = None,
    archive_resolver: "EnsemblArchiveResolver | None" = None,
    prefetch_archive: bool = True,
) -> dict[str, dict[str, Any]]:
    # Project source genes onto one target matrix universe
    output: dict[str, dict[str, Any]] = {}
    non_gene = non_gene if non_gene is not None else []
    archive_audit = archive_audit if archive_audit is not None else []
    ncbi_audit = ncbi_audit if ncbi_audit is not None else []
    mapping_status = mapping_status if mapping_status is not None else {}
    symbol_index = _gencode_symbol_index(gencode)
    gencode_ids = {str(gene) for gene in gencode["ensembl_gene_id"]}
    target_ids = {strip_ensembl_version(gene) for gene in universe} if universe else set(gencode_ids)
    manual_repairs = cfg.get("gene_mapping", {}).get("manual_symbol_repairs", {}) or {}

    # Query archive only for source Ensembl IDs absent from the target universe
    problematic_source_ids = sorted(
        {
            strip_ensembl_version(value)
            for value in term_df.get("gene_ensembl_from_source", pd.Series(dtype=str))
            if strip_ensembl_version(value) and strip_ensembl_version(value) not in target_ids
        }
    )
    if archive_resolver is not None and prefetch_archive:
        archive_resolver.prefetch(problematic_source_ids, namespace)

    for _, row in term_df.iterrows():
        if _is_non_gene_token(row.get("gene_symbol", ""), cfg):
            _count_status(mapping_status, namespace, row, "non_gene_token_drop")
            non_gene.append(_non_gene_row(namespace, row, row["gene_symbol"], "configured_non_gene_token"))
            continue
        # Keep source Ensembl IDs that already match the target universe
        source_ensembl = strip_ensembl_version(row.get("gene_ensembl_from_source", ""))
        if source_ensembl and source_ensembl in target_ids:
            _add_gene(output, row, source_ensembl)
            _count_status(mapping_status, namespace, row, "direct_source_ensembl")
            continue
        archive_info: dict[str, Any] = {}
        archive_matches: list[str] = []
        if source_ensembl and archive_resolver is not None:
            # Rescue retired Ensembl IDs only when the archive gives one target match
            archive_info = archive_resolver.lookup(source_ensembl)
            archive_matches = _archive_target_matches(archive_info, target_ids)
            if len(archive_matches) == 1:
                rescued = archive_matches[0]
                _add_gene(output, row, rescued)
                _count_status(mapping_status, namespace, row, "archive_ensembl_rescue")
                archive_audit.append(_archive_audit_row(namespace, row, source_ensembl, archive_info, rescued, "archive_ensembl_rescue"))
                continue

        # Fall back through HGNC symbol history and target annotation helpers
        current_symbol, status, matches = resolve_current_symbol(str(row["gene_symbol"]), hgnc_maps, manual_repairs)
        if status == "ambiguous_alias":
            _count_status(mapping_status, namespace, row, "ambiguous_drop")
            ambiguous.append(_ambiguous_row(namespace, row, row["gene_symbol"], "hgnc_alias_symbol", matches, [], "drop"))
            if source_ensembl:
                archive_audit.append(_archive_audit_row(namespace, row, source_ensembl, archive_info, "", _archive_unresolved_status(archive_info, archive_matches, "ambiguous_symbol")))
            continue
        if not current_symbol:
            ncbi_rescue = _try_ncbi_rescue(row, ncbi_maps or {}, symbol_index, hgnc_maps, target_ids, cfg)
            if ncbi_rescue["status"] in {"ncbi_gene_history", "ncbi_gene_info_synonym"}:
                rescued = str(ncbi_rescue["ensembl_gene_id"])
                _add_gene(output, row, rescued)
                _count_status(mapping_status, namespace, row, str(ncbi_rescue["status"]))
                ncbi_audit.append(_ncbi_audit_row(namespace, row, ncbi_rescue, rescued, "rescued"))
                if source_ensembl:
                    archive_audit.append(_archive_audit_row(namespace, row, source_ensembl, archive_info, rescued, str(ncbi_rescue["status"])))
                continue
            if ncbi_rescue["status"].startswith("ambiguous"):
                _count_status(mapping_status, namespace, row, "ambiguous_drop")
                ambiguous.append(_ambiguous_row(namespace, row, row["gene_symbol"], str(ncbi_rescue["status"]), list(ncbi_rescue["possible_symbols"]), list(ncbi_rescue["possible_ensembl_ids"]), "drop"))
                ncbi_audit.append(_ncbi_audit_row(namespace, row, ncbi_rescue, "", "ambiguous_drop"))
            else:
                _count_status(mapping_status, namespace, row, "unmapped_drop")
                unmapped.append(_unmapped_row(namespace, row, row["gene_symbol"], "", "unmapped_symbol", matches))
                if ncbi_rescue["status"] != "ncbi_not_checked":
                    ncbi_audit.append(_ncbi_audit_row(namespace, row, ncbi_rescue, "", "unmapped_drop"))
            if source_ensembl:
                archive_audit.append(_archive_audit_row(namespace, row, source_ensembl, archive_info, "", _archive_unresolved_status(archive_info, archive_matches, "unmapped_symbol")))
            continue

        candidates = _symbol_candidates(current_symbol, symbol_index, hgnc_maps)
        candidates = [candidate for candidate in candidates if candidate["ensembl_gene_id"] in target_ids]
        if not candidates:
            ncbi_rescue = _try_ncbi_rescue(row, ncbi_maps or {}, symbol_index, hgnc_maps, target_ids, cfg)
            if ncbi_rescue["status"] in {"ncbi_gene_history", "ncbi_gene_info_synonym"}:
                rescued = str(ncbi_rescue["ensembl_gene_id"])
                _add_gene(output, row, rescued)
                _count_status(mapping_status, namespace, row, str(ncbi_rescue["status"]))
                ncbi_audit.append(_ncbi_audit_row(namespace, row, ncbi_rescue, rescued, "rescued"))
                if source_ensembl:
                    archive_audit.append(_archive_audit_row(namespace, row, source_ensembl, archive_info, rescued, str(ncbi_rescue["status"])))
                continue
            _count_status(mapping_status, namespace, row, "not_in_target_universe")
            unmapped.append(_unmapped_row(namespace, row, row["gene_symbol"], current_symbol, "not_in_target_universe", []))
            if ncbi_rescue["status"] != "ncbi_not_checked":
                ncbi_audit.append(_ncbi_audit_row(namespace, row, ncbi_rescue, "", "not_in_target_universe"))
            if source_ensembl:
                archive_audit.append(_archive_audit_row(namespace, row, source_ensembl, archive_info, "", _archive_unresolved_status(archive_info, archive_matches, "not_in_target_universe")))
            continue
        chosen = _choose_gencode_candidate(candidates)
        if chosen is None:
            chosen = _choose_hgnc_ensembl_tiebreak(current_symbol, candidates, hgnc_maps, target_ids)
            if chosen is not None:
                _add_gene(output, row, chosen["ensembl_gene_id"])
                _count_status(mapping_status, namespace, row, "hgnc_ensembl_tiebreak")
                if source_ensembl:
                    archive_audit.append(_archive_audit_row(namespace, row, source_ensembl, archive_info, chosen["ensembl_gene_id"], "hgnc_ensembl_tiebreak"))
                continue
            ncbi_rescue = _try_ncbi_rescue(row, ncbi_maps or {}, symbol_index, hgnc_maps, target_ids, cfg)
            if ncbi_rescue["status"] in {"ncbi_gene_history", "ncbi_gene_info_synonym"}:
                rescued = str(ncbi_rescue["ensembl_gene_id"])
                _add_gene(output, row, rescued)
                _count_status(mapping_status, namespace, row, str(ncbi_rescue["status"]))
                ncbi_audit.append(_ncbi_audit_row(namespace, row, ncbi_rescue, rescued, "rescued"))
                if source_ensembl:
                    archive_audit.append(_archive_audit_row(namespace, row, source_ensembl, archive_info, rescued, str(ncbi_rescue["status"])))
                continue
            _count_status(mapping_status, namespace, row, "ambiguous_drop")
            ambiguous.append(
                _ambiguous_row(
                    namespace,
                    row,
                    row["gene_symbol"],
                    "gencode_gene_name",
                    [current_symbol],
                    [candidate["ensembl_gene_id"] for candidate in candidates],
                    "drop",
                )
            )
            if source_ensembl:
                archive_audit.append(_archive_audit_row(namespace, row, source_ensembl, archive_info, "", _archive_unresolved_status(archive_info, archive_matches, "ambiguous_gencode")))
            continue
        _add_gene(output, row, chosen["ensembl_gene_id"])
        mapping_status_name = _symbol_mapping_status(status)
        _count_status(mapping_status, namespace, row, mapping_status_name)
        if source_ensembl:
            archive_audit.append(_archive_audit_row(namespace, row, source_ensembl, archive_info, chosen["ensembl_gene_id"], mapping_status_name))
    return output


def resolve_current_symbol(symbol: str, hgnc_maps: dict[str, object], manual_repairs: dict[str, str] | None = None) -> tuple[str, str, list[str]]:
    symbol = str(symbol).strip()
    if not symbol:
        return "", "unmapped", []
    manual_repairs = manual_repairs or {}
    current = hgnc_maps.get("current", {})
    previous = hgnc_maps.get("previous", {})
    alias = hgnc_maps.get("alias", {})
    ambiguous_alias = hgnc_maps.get("ambiguous", {}).get("alias", {})
    if symbol in current:
        return str(current[symbol]), "exact_current_symbol", [str(current[symbol])]
    if symbol in previous:
        return str(previous[symbol]), "hgnc_previous_symbol", [str(previous[symbol])]
    if symbol in ambiguous_alias:
        return "", "ambiguous_alias", list(ambiguous_alias[symbol])
    if symbol in alias:
        return str(alias[symbol]), "hgnc_alias_symbol", [str(alias[symbol])]
    if symbol in manual_repairs:
        repaired = str(manual_repairs[symbol]).strip()
        if repaired in current:
            return repaired, "manual_symbol_repair", [repaired]
        return "", "unmapped", [repaired]
    return "", "unmapped", []


def _gencode_symbol_index(gencode: pd.DataFrame) -> dict[str, list[dict[str, str]]]:
    index: dict[str, list[dict[str, str]]] = {}
    for _, row in gencode.iterrows():
        symbol = str(row.get("gene_symbol", ""))
        if not symbol:
            continue
        index.setdefault(symbol, []).append(
            {
                "ensembl_gene_id": str(row["ensembl_gene_id"]),
                "gene_biotype": str(row.get("gene_biotype", "")),
                "is_in_expression_universe": str(row.get("is_in_expression_universe", "")) == "True"
                or bool(row.get("is_in_expression_universe", False)),
                "annotation_source": "GENCODE",
            }
        )
    return index


def _choose_gencode_candidate(candidates: list[dict[str, Any]]) -> dict[str, Any] | None:
    if len(candidates) == 1:
        return candidates[0]
    in_universe = [candidate for candidate in candidates if candidate.get("is_in_expression_universe")]
    if len(in_universe) == 1:
        return in_universe[0]
    pool = in_universe or candidates
    coding = [candidate for candidate in pool if candidate.get("gene_biotype") == "protein_coding"]
    if len(coding) == 1:
        return coding[0]
    return None


def _choose_hgnc_ensembl_tiebreak(
    current_symbol: str,
    candidates: list[dict[str, Any]],
    hgnc_maps: dict[str, object],
    target_ids: set[str],
) -> dict[str, Any] | None:
    # Prefer the official HGNC Ensembl ID when GENCODE symbol names are duplicated
    hgnc_ensembl = strip_ensembl_version(hgnc_maps.get("current_to_ensembl", {}).get(current_symbol, ""))
    if not hgnc_ensembl or hgnc_ensembl not in target_ids:
        return None
    matches = [candidate for candidate in candidates if candidate.get("ensembl_gene_id") == hgnc_ensembl]
    if not matches:
        return None
    return {**matches[0], "annotation_source": matches[0].get("annotation_source", "HGNC")}


def _symbol_candidates(current_symbol: str, symbol_index: dict[str, list[dict[str, str]]], hgnc_maps: dict[str, object]) -> list[dict[str, Any]]:
    candidates = list(symbol_index.get(current_symbol, []))
    hgnc_ensembl = str(hgnc_maps.get("current_to_ensembl", {}).get(current_symbol, ""))
    hgnc_ensembl = strip_ensembl_version(hgnc_ensembl)
    if hgnc_ensembl and hgnc_ensembl not in {candidate["ensembl_gene_id"] for candidate in candidates}:
        candidates.append(
            {
                "ensembl_gene_id": hgnc_ensembl,
                "gene_biotype": "",
                "is_in_expression_universe": True,
                "annotation_source": "HGNC",
            }
        )
    return candidates


def _try_ncbi_rescue(
    row: pd.Series,
    ncbi_maps: dict[str, object],
    symbol_index: dict[str, list[dict[str, str]]],
    hgnc_maps: dict[str, object],
    target_ids: set[str],
    cfg: dict,
) -> dict[str, object]:
    # Use NCBI only after HGNC/manual mapping fails or cannot choose one target ID
    symbol = str(row.get("gene_symbol", "")).strip()
    if not symbol or not ncbi_maps or not _ncbi_rescue_allowed_for_source(row, cfg):
        return _empty_ncbi_rescue("ncbi_not_checked")
    stages = ncbi_maps.get("stages", {})
    ambiguous = ncbi_maps.get("ambiguous", {})
    for stage in ["ncbi_gene_history", "ncbi_gene_info_synonym"]:
        if symbol in ambiguous.get(stage, {}):
            return {
                **_empty_ncbi_rescue(f"ambiguous_{stage}"),
                "stage": stage,
                "possible_symbols": list(ambiguous.get(stage, {}).get(symbol, [])),
            }
        record = stages.get(stage, {}).get(symbol)
        if not record:
            continue
        chosen, candidates = _ncbi_target_candidate(record, symbol_index, hgnc_maps, target_ids)
        base = {
            "status": stage,
            "stage": stage,
            "current_symbol": record.get("current_symbol", ""),
            "current_gene_ids": list(record.get("current_gene_ids", [])),
            "ncbi_ensembl_ids": list(record.get("ensembl_gene_ids", [])),
            "possible_symbols": [str(record.get("current_symbol", ""))],
            "possible_ensembl_ids": sorted({candidate["ensembl_gene_id"] for candidate in candidates}),
            "ensembl_gene_id": "",
        }
        if chosen is not None:
            return {**base, "ensembl_gene_id": chosen["ensembl_gene_id"]}
        if candidates:
            return {**base, "status": f"ambiguous_{stage}"}
        return {**base, "status": f"{stage}_not_in_target_universe"}
    return _empty_ncbi_rescue("ncbi_not_checked")


def _ncbi_rescue_allowed_for_source(row: pd.Series, cfg: dict) -> bool:
    allowed = cfg.get("ncbi_gene", {}).get("allowed_source_tags", [])
    if not allowed:
        return True
    return str(row.get("source_tag", "")) in {str(tag) for tag in allowed}


def _ncbi_target_candidate(
    record: dict[str, object],
    symbol_index: dict[str, list[dict[str, str]]],
    hgnc_maps: dict[str, object],
    target_ids: set[str],
) -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
    current_symbol = str(record.get("current_symbol", ""))
    candidates = [
        {
            "ensembl_gene_id": ensembl_id,
            "gene_biotype": "",
            "is_in_expression_universe": True,
            "annotation_source": "NCBI",
        }
        for ensembl_id in record.get("ensembl_gene_ids", [])
        if ensembl_id in target_ids
    ]
    candidates.extend(candidate for candidate in _symbol_candidates(current_symbol, symbol_index, hgnc_maps) if candidate["ensembl_gene_id"] in target_ids)
    candidates = _dedupe_candidates(candidates)
    chosen = _choose_gencode_candidate(candidates)
    if chosen is None:
        chosen = _choose_hgnc_ensembl_tiebreak(current_symbol, candidates, hgnc_maps, target_ids)
    if chosen is None and len({candidate["ensembl_gene_id"] for candidate in candidates}) == 1:
        chosen = candidates[0]
    return chosen, candidates


def _dedupe_candidates(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    deduped: dict[str, dict[str, Any]] = {}
    for candidate in candidates:
        ensembl_id = str(candidate.get("ensembl_gene_id", ""))
        if ensembl_id and ensembl_id not in deduped:
            deduped[ensembl_id] = candidate
    return list(deduped.values())


def _empty_ncbi_rescue(status: str) -> dict[str, object]:
    return {
        "status": status,
        "stage": "",
        "current_symbol": "",
        "current_gene_ids": [],
        "ncbi_ensembl_ids": [],
        "possible_symbols": [],
        "possible_ensembl_ids": [],
        "ensembl_gene_id": "",
    }


def _symbol_mapping_status(status: str) -> str:
    return {
        "exact_current_symbol": "hgnc_current_symbol",
        "hgnc_previous_symbol": "hgnc_previous_symbol",
        "hgnc_alias_symbol": "hgnc_alias_symbol",
        "manual_symbol_repair": "manual_symbol_repair",
    }.get(status, status)


def _count_status(mapping_status: dict[tuple[str, str, str, str], int], namespace: str, row: pd.Series, status: str) -> None:
    key = (namespace, str(row.get("family", "")), str(row.get("source_tag", "")), status)
    mapping_status[key] = mapping_status.get(key, 0) + 1


def _mapping_status_frame(mapping_status: dict[tuple[str, str, str, str], int]) -> pd.DataFrame:
    columns = ["target_namespace", "family", "source_tag", "mapping_status", "n_gene_members"]
    if not mapping_status:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(
        [
            {
                "target_namespace": namespace,
                "family": family,
                "source_tag": source_tag,
                "mapping_status": status,
                "n_gene_members": count,
            }
            for (namespace, family, source_tag, status), count in sorted(mapping_status.items())
        ],
        columns=columns,
    )


def _add_gene(output: dict[str, dict[str, Any]], row: pd.Series, gene: str) -> None:
    if not gene:
        return
    term_id = str(row["term_id"])
    if term_id not in output:
        record = row.to_dict()
        output[term_id] = {
            "genes": set(),
            "description": description_from_record(record),
            "meta": record,
            "original_genes": set(),
        }
    output[term_id]["genes"].add(str(gene))
    output[term_id]["original_genes"].add(str(row.get("gene_symbol", "")))


def _is_non_gene_token(value: object, cfg: dict) -> bool:
    text = str(value).strip()
    if not text:
        return False
    mapping_cfg = cfg.get("gene_mapping", {})
    exact = {str(token).casefold() for token in mapping_cfg.get("non_gene_tokens", [])}
    if text.casefold() in exact:
        return True
    return any(re.match(pattern, text, re.I) for pattern in mapping_cfg.get("non_gene_token_regex", []))


def _non_gene_row(namespace: str, row: pd.Series, input_gene: object, reason: str) -> dict[str, object]:
    return {
        "target_namespace": namespace,
        "family": row.get("family", ""),
        "term_id": row.get("term_id", ""),
        "source_tag": row.get("source_tag", ""),
        "input_gene": input_gene,
        "reason": reason,
    }


def _unmapped_row(namespace: str, row: pd.Series, input_gene: object, attempted: object, reason: str, matches: list[str]) -> dict[str, object]:
    return {
        "target_namespace": namespace,
        "family": row.get("family", ""),
        "term_id": row.get("term_id", ""),
        "source_tag": row.get("source_tag", ""),
        "input_gene": input_gene,
        "attempted_current_symbol": attempted,
        "reason": reason,
        "possible_matches": ",".join(matches),
    }


def _ncbi_audit_row(namespace: str, row: pd.Series, rescue: dict[str, object], resolved_ensembl: str, action: str) -> dict[str, object]:
    return {
        "target_namespace": namespace,
        "family": row.get("family", ""),
        "term_id": row.get("term_id", ""),
        "source_tag": row.get("source_tag", ""),
        "input_gene": row.get("gene_symbol", ""),
        "ncbi_stage": rescue.get("stage", ""),
        "ncbi_status": rescue.get("status", ""),
        "ncbi_current_symbol": rescue.get("current_symbol", ""),
        "ncbi_gene_ids": ",".join(rescue.get("current_gene_ids", []) or []),
        "ncbi_ensembl_ids": ",".join(rescue.get("ncbi_ensembl_ids", []) or []),
        "possible_symbols": ",".join(rescue.get("possible_symbols", []) or []),
        "possible_ensembl_ids": ",".join(rescue.get("possible_ensembl_ids", []) or []),
        "resolved_ensembl_id": resolved_ensembl,
        "action": action,
    }


def _ambiguous_row(namespace: str, row: pd.Series, input_gene: object, stage: str, symbols: list[str], ensembl: list[str], action: str) -> dict[str, object]:
    return {
        "target_namespace": namespace,
        "family": row.get("family", ""),
        "term_id": row.get("term_id", ""),
        "source_tag": row.get("source_tag", ""),
        "input_gene": input_gene,
        "mapping_stage": stage,
        "possible_symbols": ",".join(sorted(set(symbols))),
        "possible_ensembl_ids": ",".join(sorted(set(ensembl))),
        "action": action,
    }


def mapping_json(obj: object) -> str:
    return json.dumps(obj, sort_keys=True)


class EnsemblArchiveResolver:
    def __init__(
        self,
        cache_path: Path,
        endpoint: str,
        batch_size: int,
        timeout_seconds: int,
        connect_timeout_seconds: int,
        retry_attempts: int,
        retry_backoff_seconds: float,
    ) -> None:
        self.cache_path = cache_path
        self.endpoint = endpoint
        self.batch_size = batch_size
        self.timeout_seconds = timeout_seconds
        self.connect_timeout_seconds = connect_timeout_seconds
        self.retry_attempts = retry_attempts
        self.retry_backoff_seconds = retry_backoff_seconds
        self.cache = self._load_cache()
        self.changed = False

    def prefetch(self, stable_ids: list[str], namespace: str = "") -> None:
        # Fetch missing records and retry cached transient failures
        label = f" for {namespace}" if namespace else ""
        missing = [
            stable_id for stable_id in stable_ids if stable_id not in self.cache or self.cache[stable_id].get("error")
        ]
        if not missing:
            LOGGER.info("Ensembl archive audit%s: %d problematic IDs already cached", label, len(stable_ids))
            return
        n_retry = sum(1 for stable_id in missing if stable_id in self.cache and self.cache[stable_id].get("error"))
        LOGGER.info(
            "Auditing %d Ensembl IDs%s against Ensembl archive endpoint (%d cached failures will be retried)",
            len(missing),
            label,
            n_retry,
        )
        for i in range(0, len(missing), self.batch_size):
            batch = missing[i : i + self.batch_size]
            LOGGER.info("Fetching Ensembl archive batch%s %d-%d of %d", label, i + 1, i + len(batch), len(missing))
            self._fetch_batch(batch)

    def lookup(self, stable_id: str) -> dict[str, Any]:
        return self.cache.get(stable_id, self._not_found(stable_id))

    def save(self) -> None:
        # Persist cache only when new records were added
        if self.changed:
            write_json(self.cache, self.cache_path)

    def _fetch_batch(self, stable_ids: list[str]) -> None:
        # Query Ensembl archive with retries before falling back to smaller batches
        payload: Any | None = None
        last_error: Exception | None = None
        for attempt in range(1, self.retry_attempts + 1):
            try:
                payload = self._post_archive_batch(stable_ids)
                last_error = None
                break
            except Exception as exc:
                last_error = exc
                LOGGER.warning(
                    "Ensembl archive audit failed for %d IDs on attempt %d/%d: %s",
                    len(stable_ids),
                    attempt,
                    self.retry_attempts,
                    exc,
                )
                if attempt < self.retry_attempts:
                    time.sleep(self.retry_backoff_seconds * attempt)

        if last_error is not None:
            if len(stable_ids) > 1:
                midpoint = max(1, len(stable_ids) // 2)
                LOGGER.warning(
                    "Splitting failed Ensembl archive batch of %d IDs into %d and %d IDs",
                    len(stable_ids),
                    midpoint,
                    len(stable_ids) - midpoint,
                )
                self._fetch_batch(stable_ids[:midpoint])
                self._fetch_batch(stable_ids[midpoint:])
                return
            LOGGER.warning("Ensembl archive audit failed for %s after retries: %s", stable_ids[0], last_error)
            for stable_id in stable_ids:
                self.cache[stable_id] = {"query_id": stable_id, "found": False, "error": str(last_error)}
            self.changed = True
            return

        returned: set[str] = set()
        if isinstance(payload, list):
            for item in payload:
                if not isinstance(item, dict) or "id" not in item:
                    continue
                stable_id = strip_ensembl_version(item["id"])
                returned.add(stable_id)
                self.cache[stable_id] = _normalize_archive_record(stable_id, item)
        for stable_id in set(stable_ids) - returned:
            self.cache[stable_id] = self._not_found(stable_id)
        self.changed = True

    def _post_archive_batch(self, stable_ids: list[str]) -> Any:
        # Send one archive request and return parsed JSON
        response = requests.post(
            self.endpoint,
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            json={"id": stable_ids},
            timeout=(self.connect_timeout_seconds, self.timeout_seconds),
        )
        response.raise_for_status()
        return response.json()

    def _load_cache(self) -> dict[str, dict[str, Any]]:
        if not self.cache_path.exists():
            return {}
        try:
            with self.cache_path.open("r", encoding="utf-8") as handle:
                cache = json.load(handle)
        except Exception as exc:
            LOGGER.warning("Could not read Ensembl archive cache %s: %s", self.cache_path, exc)
            return {}
        return cache if isinstance(cache, dict) else {}

    @staticmethod
    def _not_found(stable_id: str) -> dict[str, Any]:
        return {"query_id": stable_id, "found": False}


def _build_archive_resolver(source_dir: Path | None, cfg: dict) -> EnsemblArchiveResolver | None:
    archive_cfg = cfg.get("ensembl_archive", {})
    if not archive_cfg.get("enabled", True):
        return None
    source_dir = source_dir or Path("data/gmt_sources")
    cache_path = source_dir / "ensembl_archive" / "archive_id_cache.json"
    ensure_dirs(cache_path.parent)
    return EnsemblArchiveResolver(
        cache_path=cache_path,
        endpoint=archive_cfg.get("endpoint", "https://rest.ensembl.org/archive/id"),
        batch_size=int(archive_cfg.get("batch_size", 50)),
        timeout_seconds=int(archive_cfg.get("timeout_seconds", 120)),
        connect_timeout_seconds=int(archive_cfg.get("connect_timeout_seconds", 30)),
        retry_attempts=int(archive_cfg.get("retry_attempts", 3)),
        retry_backoff_seconds=float(archive_cfg.get("retry_backoff_seconds", 5)),
    )


def _normalize_archive_record(stable_id: str, item: dict[str, Any]) -> dict[str, Any]:
    # Keep only stable fields needed for rescue and audit
    possible_replacements = []
    for replacement in item.get("possible_replacement") or []:
        if isinstance(replacement, dict):
            replacement_id = replacement.get("id") or replacement.get("stable_id")
        else:
            replacement_id = replacement
        replacement_id = strip_ensembl_version(replacement_id)
        if replacement_id:
            possible_replacements.append(replacement_id)
    latest = strip_ensembl_version(item.get("latest", ""))
    return {
        "query_id": stable_id,
        "found": True,
        "archive_id": strip_ensembl_version(item.get("id", stable_id)),
        "is_current": str(item.get("is_current", "")),
        "release": str(item.get("release", "")),
        "assembly": str(item.get("assembly", "")),
        "latest": latest,
        "possible_replacement": sorted(set(possible_replacements)),
    }


def _archive_target_matches(archive_info: dict[str, Any], target_ids: set[str]) -> list[str]:
    # Intersect archive candidates with the target matrix universe
    if not archive_info or not archive_info.get("found"):
        return []
    candidates = set(archive_info.get("possible_replacement") or [])
    latest = strip_ensembl_version(archive_info.get("latest", ""))
    archive_id = strip_ensembl_version(archive_info.get("archive_id", ""))
    if latest:
        candidates.add(latest)
    if str(archive_info.get("is_current", "")) in {"1", "true", "True"} and archive_id:
        candidates.add(archive_id)
    return sorted(candidate for candidate in candidates if candidate in target_ids)


def _archive_unresolved_status(archive_info: dict[str, Any], archive_matches: list[str], fallback_status: str) -> str:
    if archive_matches:
        return f"archive_ambiguous_target_replacement__{fallback_status}"
    if archive_info.get("error"):
        return f"archive_query_failed__{fallback_status}"
    if archive_info and not archive_info.get("found"):
        return f"archive_not_found__{fallback_status}"
    if archive_info:
        return f"archive_no_target_replacement__{fallback_status}"
    return f"archive_not_checked__{fallback_status}"


def _archive_audit_row(
    namespace: str,
    row: pd.Series,
    source_ensembl: str,
    archive_info: dict[str, Any],
    rescued_ensembl: str,
    rescue_status: str,
) -> dict[str, object]:
    return {
        "target_namespace": namespace,
        "family": row.get("family", ""),
        "term_id": row.get("term_id", ""),
        "source_tag": row.get("source_tag", ""),
        "input_symbol": row.get("gene_symbol", ""),
        "source_ensembl_id": row.get("gene_ensembl_from_source", ""),
        "source_ensembl_stable": source_ensembl,
        "archive_found": bool(archive_info.get("found")) if archive_info else False,
        "archive_is_current": archive_info.get("is_current", ""),
        "archive_release": archive_info.get("release", ""),
        "archive_assembly": archive_info.get("assembly", ""),
        "archive_latest": archive_info.get("latest", ""),
        "archive_possible_replacements": ",".join(archive_info.get("possible_replacement", []) or []),
        "archive_error": archive_info.get("error", ""),
        "rescued_ensembl_id": rescued_ensembl,
        "rescue_status": rescue_status,
    }
