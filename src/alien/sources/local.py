"""Local gene-set files: symbol GMTs and canonical membership tables."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from ..cache import ResourceCache
from ..config.schema import CanonicalTsvSource, LocalSourceBase, SymbolGmtSource
from ..errors import SourceError
from ..model import is_ensembl_gene_id, sanitize_term_id
from ..schema import CANONICAL_COLUMNS
from .frames import concat_memberships, rows_from_gmt
from .registry import register

LOGGER = logging.getLogger("alien")
REQUIRED_CANONICAL = ("term_id", "gene_symbol")
NAME_FALLBACKS = (("original_name", "term_id"), ("display_name", "original_name"), ("description", "display_name"))


def source_paths(spec: LocalSourceBase, source_dir: Path) -> list[Path]:
    """Configured paths; a symbol GMT without paths reads every ``*.gmt`` in ``source_dir/<name>``."""
    paths = spec.input_paths()
    if not paths and isinstance(spec, SymbolGmtSource):
        directory = source_dir / (spec.name or "symbol_gmt")
        paths = sorted(directory.glob("*.gmt"))
        LOGGER.info("Source %s has no paths; using %d GMT file(s) from %s", spec.name, len(paths), directory)
    for path in paths:
        if not path.exists():
            raise SourceError(f"Source file {path} does not exist.")
    return paths


@register("symbol_gmt", "enrichr_gmt")
def load_symbol_gmt(spec: SymbolGmtSource, cache: ResourceCache) -> pd.DataFrame:
    return concat_memberships(read_symbol_gmt(path, spec) for path in source_paths(spec, cache.root))


def read_symbol_gmt(path: Path, spec: SymbolGmtSource) -> pd.DataFrame:
    source = spec.source or spec.name or "symbol_gmt"
    source_tag = spec.source_tag or sanitize_term_id(source).upper()
    constants = {
        "source": source,
        "source_tag": source_tag,
        "collection": spec.collection or path.stem,
        "subcollection": spec.subcollection or "",
        "family": spec.family or "biology_process_pathway",
        "aspect": spec.aspect or "library",
        "db_version": spec.db_version or "",
        "source_url": spec.source_url or "",
        "source_license_note": spec.source_license_note or "",
        "metadata_json": "{}",
    }
    return rows_from_gmt(path, sanitize_term_id(spec.term_prefix or source_tag), constants)


@register("canonical_tsv")
def load_canonical_tsv(spec: CanonicalTsvSource, cache: ResourceCache) -> pd.DataFrame:
    return concat_memberships(read_canonical_memberships(path, spec) for path in source_paths(spec, cache.root))


def _canonical_defaults(spec: CanonicalTsvSource, path: Path, first_source: str) -> dict[str, str]:
    """Values for empty canonical columns: the spec's settings, else built-in defaults."""
    configured = {
        "source_tag": spec.source_tag,
        "family": spec.family,
        "aspect": spec.aspect,
        "collection": spec.collection,
        "subcollection": spec.subcollection,
    }
    fallback = {
        "source_tag": sanitize_term_id(first_source).upper(),
        "family": "biology_process_pathway",
        "aspect": "library",
        "collection": path.stem,
        "subcollection": "",
    }
    return {column: value if value else fallback[column] for column, value in configured.items()}


def _source_ensembl_ids(raw: pd.DataFrame, current: pd.Series) -> pd.Series:
    """Source Ensembl IDs: gene_id rows whose namespace says Ensembl, or that look like Ensembl IDs."""
    if "gene_id" not in raw:
        return current
    if "gene_id_namespace" in raw:
        ensembl = raw["gene_id_namespace"].str.lower().str.contains("ensembl|ensg", regex=True, na=False)
        return raw["gene_id"].where(ensembl, current)
    return raw["gene_id"].map(lambda value: value if is_ensembl_gene_id(value) else "")


def read_canonical_memberships(path: Path, spec: CanonicalTsvSource) -> pd.DataFrame:
    raw = pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)
    missing = [column for column in REQUIRED_CANONICAL if column not in raw]
    if missing:
        raise SourceError(f"{path} is missing canonical columns: {', '.join(missing)}")
    frame = pd.DataFrame({column: raw.get(column, "") for column in CANONICAL_COLUMNS})
    frame["term_id"] = frame["term_id"].map(sanitize_term_id)
    frame["source"] = frame["source"].replace("", spec.source or spec.name or "canonical")
    first_source = str(frame["source"].iloc[0]) if len(frame) else "canonical"
    for column, value in _canonical_defaults(spec, path, first_source).items():
        frame[column] = frame[column].replace("", value)
    for column, fallback in (
        ("original_name", "term_id"),
        ("display_name", "original_name"),
        ("description", "display_name"),
    ):
        frame[column] = frame[column].where(frame[column].ne(""), frame[fallback])
    frame["gene_ensembl_from_source"] = _source_ensembl_ids(raw, frame["gene_ensembl_from_source"])
    frame["metadata_json"] = frame["metadata_json"].replace("", "{}")
    return frame[list(CANONICAL_COLUMNS)]
