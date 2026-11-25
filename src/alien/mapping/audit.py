"""Audit tables derived from mapping outcomes: one builder per table, rows in membership order."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any

import numpy as np
import pandas as pd

from ..schema import (
    AMBIGUOUS_COLUMNS,
    ARCHIVE_AUDIT_COLUMNS,
    NCBI_AUDIT_COLUMNS,
    NON_GENE_COLUMNS,
    UNMAPPED_COLUMNS,
)
from .resolver import MappingOutcome

Fields = dict[str, Any]


def _joined(values: Sequence[str], *, sort: bool = False) -> str:
    return ",".join(sorted(set(values)) if sort else values)


def audit_table(
    rows: pd.DataFrame,
    namespace: str,
    codes: np.ndarray,
    per_code: Sequence[Fields | None],
    columns: Sequence[str],
) -> pd.DataFrame:
    """Rows whose key yields audit fields (``per_code[code]``), with membership identity columns prepended."""
    selected = np.flatnonzero(np.fromiter((per_code[code] is not None for code in codes), dtype=bool, count=len(codes)))
    if not len(selected):
        return pd.DataFrame(columns=list(columns))
    base = rows.iloc[selected]
    data: dict[str, Any] = {
        "target_namespace": namespace,
        "family": base["family"].to_numpy(),
        "term_id": base["term_id"].to_numpy(),
        "source_tag": base["source_tag"].to_numpy(),
    }
    derived = pd.DataFrame.from_records([per_code[code] for code in codes[selected]])
    for column in columns:
        if column in data:
            continue
        if column in derived:
            data[column] = derived[column].to_numpy()
        elif column in ("input_gene", "input_symbol"):
            data[column] = base["gene_symbol"].to_numpy()
        elif column == "source_ensembl_id":
            data[column] = base["gene_ensembl_from_source"].to_numpy()
    return pd.DataFrame(data, columns=list(columns))


def non_gene_fields(outcome: MappingOutcome) -> Fields | None:
    return {"reason": "configured_non_gene_token"} if outcome.status == "non_gene_token_drop" else None


def unmapped_fields(outcome: MappingOutcome) -> Fields | None:
    detail = outcome.unmapped
    if detail is None:
        return None
    return {
        "attempted_current_symbol": detail.attempted_symbol,
        "reason": detail.reason,
        "possible_matches": _joined(detail.matches),
    }


def ambiguous_fields(outcome: MappingOutcome) -> Fields | None:
    detail = outcome.ambiguous
    if detail is None:
        return None
    return {
        "mapping_stage": detail.stage,
        "possible_symbols": _joined(detail.symbols, sort=True),
        "possible_ensembl_ids": _joined(detail.ensembl_ids, sort=True),
        "action": detail.action,
    }


def ncbi_fields(outcome: MappingOutcome) -> Fields | None:
    if outcome.ncbi is None:
        return None
    rescue, action = outcome.ncbi
    return {
        "ncbi_stage": rescue.stage,
        "ncbi_status": rescue.status,
        "ncbi_current_symbol": rescue.current_symbol,
        "ncbi_gene_ids": _joined(rescue.current_gene_ids),
        "ncbi_ensembl_ids": _joined(rescue.ncbi_ensembl_ids),
        "possible_symbols": _joined(rescue.possible_symbols),
        "possible_ensembl_ids": _joined(rescue.possible_ensembl_ids),
        "resolved_ensembl_id": outcome.ensembl_id if action == "rescued" else "",
        "action": action,
    }


def archive_fields(outcome: MappingOutcome, source_stable: str) -> Fields | None:
    if outcome.archive is None:
        return None
    info, rescue_status = outcome.archive
    return {
        "source_ensembl_stable": source_stable,
        "archive_found": bool(info.get("found")) if info else False,
        "archive_is_current": info.get("is_current", ""),
        "archive_release": info.get("release", ""),
        "archive_assembly": info.get("assembly", ""),
        "archive_latest": info.get("latest", ""),
        "archive_possible_replacements": ",".join(info.get("possible_replacement", []) or []),
        "archive_error": info.get("error", ""),
        "rescued_ensembl_id": outcome.ensembl_id if outcome.status.is_mapped else "",
        "rescue_status": rescue_status,
    }


TABLES: dict[str, tuple[Callable[[MappingOutcome], Fields | None], Sequence[str]]] = {
    "unmapped": (unmapped_fields, UNMAPPED_COLUMNS),
    "ambiguous": (ambiguous_fields, AMBIGUOUS_COLUMNS),
    "non_gene": (non_gene_fields, NON_GENE_COLUMNS),
    "ncbi": (ncbi_fields, NCBI_AUDIT_COLUMNS),
}
ARCHIVE_COLUMNS = ARCHIVE_AUDIT_COLUMNS
