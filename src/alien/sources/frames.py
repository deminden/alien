"""Helpers for canonical membership frames."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from pathlib import Path

import pandas as pd

from ..io.gmt import parse_gmt
from ..model import sanitize_term_id
from ..schema import CANONICAL_COLUMNS


def empty_memberships() -> pd.DataFrame:
    return pd.DataFrame(columns=list(CANONICAL_COLUMNS))


def concat_memberships(frames: Iterable[pd.DataFrame | None]) -> pd.DataFrame:
    kept = [frame for frame in frames if frame is not None and not frame.empty]
    if not kept:
        return empty_memberships()
    combined = pd.concat(kept, ignore_index=True)
    for column in CANONICAL_COLUMNS:
        if column not in combined:
            combined[column] = ""
    return combined


def column_or(raw: pd.DataFrame, name: str, default: str = "") -> pd.Series:
    if name in raw:
        return raw[name].fillna("")
    return pd.Series([default] * len(raw), index=raw.index, dtype=object)


def rows_from_gmt(path: Path, term_prefix: str, constants: Mapping[str, str]) -> pd.DataFrame:
    """One canonical row per (term, gene) of a symbol GMT; term IDs are ``<prefix>__<sanitized name>``."""
    rows: list[dict[str, str]] = []
    for record in parse_gmt(path):
        term_id = f"{term_prefix}__{sanitize_term_id(record.term_id)}"
        base = {
            **constants,
            "term_id": term_id,
            "original_name": record.term_id,
            "display_name": record.term_id,
            "description": record.description,
            "gene_ensembl_from_source": "",
        }
        rows.extend({**base, "gene_symbol": gene} for gene in record.genes)
    return pd.DataFrame(rows, columns=list(CANONICAL_COLUMNS))
