"""HGNC symbol history: current, previous and alias symbols (human only)."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

import pandas as pd

from ..cache import ResourceCache
from ..config.schema import HgncConfig
from ..errors import SourceError
from ..model import ensembl_stable_id

# HGNC moved its downloads to Google Cloud Storage; the EBI FTP path is kept as a fallback only.
HGNC_URL: Final = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
HGNC_FALLBACK_URLS: Final = ("https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt",)
HGNC_ARCHIVE_URLS: Final = (
    "https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/hgnc_complete_set_{date}.txt",
    "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_{date}.txt",
)


@dataclass(frozen=True)
class HgncMaps:
    """Symbol lookups. Ambiguous previous and alias symbols map to the sorted list of candidates."""

    current: frozenset[str] = frozenset()
    previous: dict[str, str] = field(default_factory=dict)
    alias: dict[str, str] = field(default_factory=dict)
    ambiguous_previous: dict[str, list[str]] = field(default_factory=dict)
    ambiguous_alias: dict[str, list[str]] = field(default_factory=dict)
    current_to_ensembl: dict[str, str] = field(default_factory=dict)

    @property
    def current_ensembl_ids(self) -> frozenset[str]:
        return frozenset(stable for value in self.current_to_ensembl.values() if (stable := ensembl_stable_id(value)))


def hgnc_urls(config: HgncConfig) -> list[str]:
    if config.url:
        return [config.url]
    if config.archive_date:
        return [url.format(date=config.archive_date) for url in HGNC_ARCHIVE_URLS]
    return [HGNC_URL, *HGNC_FALLBACK_URLS]


def hgnc_path(source_dir: Path, config: HgncConfig) -> Path:
    name = f"hgnc_complete_set_{config.archive_date}.txt" if config.archive_date else "hgnc_complete_set.txt"
    return source_dir / "hgnc" / name


def load_hgnc_table(cache: ResourceCache, config: HgncConfig) -> pd.DataFrame:
    path = cache.fetch(
        "hgnc",
        hgnc_urls(config),
        hgnc_path(cache.root, config),
        upstream_version=config.archive_date or "",
    )
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if "symbol" not in table:
        raise SourceError(f"{path} is not an HGNC complete set (no 'symbol' column).")
    return table


def _split(value: object) -> list[str]:
    text = "" if value is None else str(value).strip()
    return [item.strip() for item in text.split("|") if item.strip()] if text else []


def _column(table: pd.DataFrame, name: str) -> pd.Series:
    return table[name] if name in table else pd.Series([""] * len(table), index=table.index)


def _split_unambiguous(lookup: dict[str, set[str]]) -> tuple[dict[str, str], dict[str, list[str]]]:
    unique = {key: next(iter(values)) for key, values in lookup.items() if len(values) == 1}
    ambiguous = {key: sorted(values) for key, values in lookup.items() if len(values) > 1}
    return unique, ambiguous


def build_hgnc_maps(table: pd.DataFrame) -> HgncMaps:
    current: set[str] = set()
    history: dict[str, dict[str, set[str]]] = {"prev_symbol": {}, "alias_symbol": {}}
    ensembl: dict[str, str] = {}
    rows = zip(*(_column(table, name) for name in ("symbol", "ensembl_gene_id", *history)), strict=True)
    for symbol_raw, ensembl_raw, *lists in rows:
        symbol = str(symbol_raw).strip()
        if not symbol:
            continue
        current.add(symbol)
        if ensembl_raw:
            ensembl[symbol] = str(ensembl_raw).strip()
        for lookup, raw in zip(history.values(), lists, strict=True):
            for name in _split(raw):
                lookup.setdefault(name, set()).add(symbol)
    previous, ambiguous_previous = _split_unambiguous(history["prev_symbol"])
    alias, ambiguous_alias = _split_unambiguous(history["alias_symbol"])
    return HgncMaps(frozenset(current), previous, alias, ambiguous_previous, ambiguous_alias, ensembl)
