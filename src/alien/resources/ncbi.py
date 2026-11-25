"""NCBI Gene history and synonyms: a lower-confidence rescue layer used only after HGNC fails (human only)."""

from __future__ import annotations

import json
import logging
import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

import pandas as pd

from ..cache import ResourceCache, derived_cache_path, file_sha256
from ..config.schema import NcbiGeneConfig
from ..io.atomic import atomic_write_json

LOGGER = logging.getLogger("alien")
HUMAN_TAXON: Final = "9606"
STAGES: Final = ("ncbi_gene_history", "ncbi_gene_info_synonym")
RESCUE_COLUMNS: Final = ("input_symbol", "stage", "current_symbol", "current_gene_id", "ensembl_gene_ids")


@dataclass(frozen=True)
class NcbiRecord:
    current_symbol: str
    current_gene_ids: tuple[str, ...]
    ensembl_gene_ids: tuple[str, ...]


@dataclass(frozen=True)
class NcbiMaps:
    """Per stage: unambiguous input symbol -> record, and ambiguous input symbol -> candidate symbols."""

    stages: dict[str, dict[str, NcbiRecord]] = field(default_factory=dict)
    ambiguous: dict[str, dict[str, list[str]]] = field(default_factory=dict)

    def __bool__(self) -> bool:
        return bool(self.stages or self.ambiguous)

    def to_json(self) -> dict[str, Any]:
        return {
            "stages": {
                stage: {
                    symbol: {
                        "current_symbol": record.current_symbol,
                        "current_gene_ids": list(record.current_gene_ids),
                        "ensembl_gene_ids": list(record.ensembl_gene_ids),
                    }
                    for symbol, record in records.items()
                }
                for stage, records in self.stages.items()
            },
            "ambiguous": self.ambiguous,
        }

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> NcbiMaps:
        stages = {
            stage: {
                symbol: NcbiRecord(
                    str(record["current_symbol"]),
                    tuple(record["current_gene_ids"]),
                    tuple(record["ensembl_gene_ids"]),
                )
                for symbol, record in records.items()
            }
            for stage, records in data["stages"].items()
        }
        return cls(stages=stages, ambiguous={stage: dict(values) for stage, values in data["ambiguous"].items()})


def load_ncbi_maps(config: NcbiGeneConfig, cache: ResourceCache) -> NcbiMaps:
    if not config.enabled:
        return NcbiMaps()
    ncbi_dir = cache.root / "ncbi_gene"
    info_path = cache.fetch("ncbi_gene:gene_info", config.gene_info_url, ncbi_dir / "Homo_sapiens.gene_info.gz")
    history_path = cache.fetch("ncbi_gene:gene_history", config.gene_history_url, ncbi_dir / "gene_history.gz")
    cache_path = derived_cache_path(
        ncbi_dir,
        "ncbi_gene_rescue_maps",
        {
            "gene_info_sha256": file_sha256(info_path),
            "gene_history_sha256": file_sha256(history_path),
            "use_gene_info_synonyms": config.use_gene_info_synonyms,
        },
        suffix=".json",
    )
    if cache_path.exists() and not cache.force:
        try:
            return NcbiMaps.from_json(json.loads(cache_path.read_text(encoding="utf-8")))
        except (OSError, ValueError, KeyError, TypeError) as exc:
            LOGGER.warning("Ignoring unreadable NCBI Gene map cache %s: %s", cache_path, exc)
    LOGGER.info("Parsing NCBI Gene info and history rescue tables")
    gene_info = read_gene_info(info_path)
    rows = build_rescue_rows(gene_info, history_path, config.use_gene_info_synonyms)
    maps = maps_from_rows(rows)
    with cache.file_lock(cache_path):
        atomic_write_json(cache_path, maps.to_json())
    return maps


def _ensembl_ids_from_dbxrefs(value: object) -> str:
    matches = re.findall(r"(?:^|\|)Ensembl:(ENSG[0-9]+)", "" if value is None else str(value))
    return ",".join(sorted(set(matches)))


def _split_synonyms(value: object) -> list[str]:
    text = "" if value is None else str(value).strip()
    if not text or text == "-":
        return []
    return [item.strip() for item in text.split("|") if item.strip() and item.strip() != "-"]


def _human(frame: pd.DataFrame) -> pd.DataFrame:
    tax_column = "#tax_id" if "#tax_id" in frame else "tax_id"
    return frame[frame[tax_column].eq(HUMAN_TAXON)]


def read_gene_info(path: Path) -> pd.DataFrame:
    raw = _human(pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)).copy()
    xrefs = raw["dbXrefs"] if "dbXrefs" in raw else pd.Series([""] * len(raw), index=raw.index)
    raw["ensembl_gene_ids"] = xrefs.map(_ensembl_ids_from_dbxrefs)
    return raw


def _current_records(gene_info: pd.DataFrame) -> Iterator[tuple[str, dict[str, str], str]]:
    for gene_id_raw, symbol_raw, ensembl, synonyms in zip(
        gene_info.get("GeneID", pd.Series(dtype=str)),
        gene_info.get("Symbol", pd.Series(dtype=str)),
        gene_info["ensembl_gene_ids"],
        gene_info.get("Synonyms", pd.Series([""] * len(gene_info), index=gene_info.index)),
        strict=True,
    ):
        gene_id, symbol = str(gene_id_raw).strip(), str(symbol_raw).strip()
        if symbol and symbol != "-":
            yield (
                gene_id,
                {"current_symbol": symbol, "current_gene_id": gene_id, "ensembl_gene_ids": str(ensembl)},
                synonyms,
            )


def _usable(text: str) -> bool:
    return bool(text) and text != "-"


def _history_rows(history_path: Path, current_by_gene_id: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for chunk in pd.read_csv(history_path, sep="\t", dtype=str, keep_default_na=False, chunksize=500_000):
        human = _human(chunk)
        for symbol_raw, gene_id_raw in zip(human["Discontinued_Symbol"], human["GeneID"], strict=True):
            symbol, gene_id = str(symbol_raw).strip(), str(gene_id_raw).strip()
            current = current_by_gene_id.get(gene_id) if _usable(symbol) and _usable(gene_id) else None
            if current:
                rows.append({"input_symbol": symbol, "stage": "ncbi_gene_history", **current})
    return rows


def build_rescue_rows(gene_info: pd.DataFrame, history_path: Path, use_synonyms: bool) -> pd.DataFrame:
    current_by_gene_id = {gene_id: record for gene_id, record, _ in _current_records(gene_info) if gene_id}
    rows = _history_rows(history_path, current_by_gene_id)
    if use_synonyms:
        for _, current, synonyms in _current_records(gene_info):
            rows.extend(
                {"input_symbol": synonym, "stage": "ncbi_gene_info_synonym", **current}
                for synonym in _split_synonyms(synonyms)
                if synonym != current["current_symbol"]
            )
    return pd.DataFrame(rows, columns=list(RESCUE_COLUMNS)).drop_duplicates()


def maps_from_rows(rows: pd.DataFrame) -> NcbiMaps:
    stages: dict[str, dict[str, NcbiRecord]] = {}
    ambiguous: dict[str, dict[str, list[str]]] = {}
    for stage, stage_rows in rows.groupby("stage"):
        stage_map: dict[str, NcbiRecord] = {}
        stage_ambiguous: dict[str, list[str]] = {}
        for input_symbol, symbol_rows in stage_rows.groupby("input_symbol"):
            targets = symbol_rows[["current_symbol", "current_gene_id", "ensembl_gene_ids"]].drop_duplicates()
            symbols = sorted(set(targets["current_symbol"]))
            if len(symbols) != 1:
                stage_ambiguous[str(input_symbol)] = symbols
                continue
            ensembl_ids = sorted(
                {item for value in targets["ensembl_gene_ids"] for item in str(value).split(",") if item}
            )
            stage_map[str(input_symbol)] = NcbiRecord(
                symbols[0], tuple(sorted(set(targets["current_gene_id"]))), tuple(ensembl_ids)
            )
        stages[str(stage)] = stage_map
        ambiguous[str(stage)] = stage_ambiguous
    return NcbiMaps(stages=stages, ambiguous=ambiguous)
