"""Helpers for canonical membership frames."""

from __future__ import annotations

import logging
from collections.abc import Iterable, Mapping
from pathlib import Path

import pandas as pd

from ..io.gmt import parse_gmt
from ..model import sanitize_term_id
from ..schema import CANONICAL_COLUMNS

LOGGER = logging.getLogger("alien")


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
    seen: dict[str, int] = {}
    for record in parse_gmt(path):
        name = record.name or record.term_id
        seen[name] = seen.get(name, 0) + 1
        term_id = f"{term_prefix}__{sanitize_term_id(record.term_id)}"
        base = {
            **constants,
            "term_id": term_id,
            "original_name": name,
            "display_name": name,
            "description": record.description,
            "gene_ensembl_from_source": "",
        }
        rows.extend({**base, "gene_symbol": gene} for gene in record.genes)
    repeated = sorted(name for name, count in seen.items() if count > 1)
    if repeated:
        LOGGER.warning(
            "%s lists %d term name(s) more than once; their members are combined (first: %s).",
            path,
            len(repeated),
            repeated[0],
        )
    return pd.DataFrame(rows, columns=list(CANONICAL_COLUMNS))
