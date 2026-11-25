"""Target preparation: annotation GTFs, metadata fallbacks, output-gene symbol rows and the final Target."""

from __future__ import annotations

import logging
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Final

import pandas as pd

from ..cache import ResourceCache, derived_cache_path, file_sha256
from ..config.schema import AnnotationFile, TargetSpec
from ..errors import ConfigError, ResourceError
from ..io.gtf import parse_gtf_genes
from ..io.tables import write_table
from ..model import GeneIdList, Target, sanitize_term_id
from ..schema import ANNOTATION_COLUMNS
from .output_genes import read_gene_id_list

LOGGER = logging.getLogger("alien")
GENCODE_HUMAN_BASE_URL: Final = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"
SOURCE_PRIORITY: Final = {"annotation_gtf": 0, "metadata_fallback": 1, "output_genes": 2, "HGNC": 3}


def normalize_gencode_version(version: object) -> str:
    text = str(version).strip()
    text = text.removeprefix("v").removeprefix("V")
    if not re.fullmatch(r"[0-9]+", text):
        raise ConfigError(f"GENCODE version must be numeric, got {version!r}.")
    return text


def gencode_annotation_url(version: object) -> str:
    version = normalize_gencode_version(version)
    return f"{GENCODE_HUMAN_BASE_URL}/release_{version}/gencode.v{version}.annotation.gtf.gz"


def annotation_location(
    target_name: str, spec: AnnotationFile, role: str, source_dir: Path
) -> tuple[Path, str, str, str]:
    """(local path, download URL or "", label, version) for one annotation file."""
    version = spec.version or ""
    source = spec.source or "GTF"
    if source.upper() == "GENCODE" and version:
        version = normalize_gencode_version(version)
    label = f"{source} v{version}" if version else source
    if spec.path is not None:
        return spec.path, "", label, version
    if spec.url:
        file_name = spec.file_name or Path(spec.url.split("?", 1)[0]).name or f"{sanitize_term_id(target_name)}.gtf.gz"
        path = source_dir / "targets" / sanitize_term_id(target_name) / sanitize_term_id(role) / file_name
        return path, spec.url, label, version
    if source.upper() != "GENCODE" or not version:
        raise ConfigError(
            f"Target {target_name} uses annotation source {source!r}; provide annotation.path or annotation.url "
            "for non-GENCODE target annotations."
        )
    path = source_dir / "gencode" / f"gencode.v{version}.annotation.gtf.gz"
    return path, gencode_annotation_url(version), f"GENCODE v{version}", version


def load_annotation_file(
    target_name: str, spec: AnnotationFile, role: str, cache: ResourceCache
) -> tuple[pd.DataFrame, str, Path]:
    path, url, label, version = annotation_location(target_name, spec, role, cache.root)
    if url:
        cache.fetch(f"target:{target_name}:{role}", url, path, upstream_version=label)
    elif not path.exists():
        raise ResourceError(f"Target {target_name} annotation file {path} does not exist.")
    LOGGER.info("Parsing %s target annotation genes from %s.", label, path)
    return parse_gtf_cached(path, version, spec, cache), label, path


def parse_gtf_cached(path: Path, version: str, spec: AnnotationFile, cache: ResourceCache) -> pd.DataFrame:
    """Parse a GTF once per machine: the gene table is cached by GTF sha256 and attribute mapping."""
    attributes = spec.attributes.model_dump()
    key = {"gtf_sha256": file_sha256(path), "version": version, "attributes": attributes}
    cache_path = derived_cache_path(cache.root / "gtf_cache", "genes", key)
    if cache_path.exists() and not cache.force:
        frame = pd.read_parquet(cache_path)
        frame["is_in_output_genes"] = frame["is_in_output_genes"].astype(bool)
    else:
        frame = parse_gtf_genes(path, version, attributes)
        with cache.file_lock(cache_path):
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            frame.to_parquet(cache_path, index=False)
    frame["source_gtf"] = str(path)
    return frame


def load_metadata_fallbacks(
    target: TargetSpec, output_ids: frozenset[str], primary_ids: frozenset[str], cache: ResourceCache
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    """Fill measured output IDs that the primary annotation lacks from fallback annotations, in order."""
    records: list[dict[str, Any]] = []
    frames: list[pd.DataFrame] = []
    covered = set(primary_ids)
    for index, spec in enumerate(target.annotation.metadata_fallbacks, start=1):
        frame, label, path = load_annotation_file(target.name, spec, f"metadata_fallback_{index}", cache)
        fallback_ids = set(frame["ensembl_gene_id"].astype(str))
        missing_before = output_ids - covered if output_ids else set()
        add_ids = missing_before & fallback_ids
        selected = frame[frame["ensembl_gene_id"].isin(add_ids)].copy() if add_ids else frame.iloc[0:0].copy()
        if not selected.empty:
            selected["target_id_metadata_source"] = "metadata_fallback"
            frames.append(selected)
            covered.update(add_ids)
        records.append(
            {
                "target_namespace": target.name,
                "fallback_index": index,
                "mode": spec.mode,
                "annotation": label,
                "annotation_path": str(path),
                "n_fallback_genes": len(fallback_ids),
                "n_missing_before": len(missing_before),
                "n_rows_added": len(selected),
                "n_missing_after": len(output_ids - covered) if output_ids else 0,
            }
        )
    combined = pd.concat(frames, ignore_index=True).drop_duplicates() if frames else pd.DataFrame()
    return combined, records


def add_output_gene_symbols(annotation: pd.DataFrame, output_genes: GeneIdList) -> pd.DataFrame:
    """Add rows for output genes the annotation lacks, using symbols supplied with the output-gene list."""
    table = output_genes.table
    if not output_genes.ids or table.empty or "gene_symbol" not in table:
        return annotation
    existing = set(annotation["ensembl_gene_id"].astype(str))
    rows: list[dict[str, object]] = []
    seen: set[str] = set()
    source = f"output_genes:{output_genes.source_path}" if output_genes.source_path else "output_genes"
    for input_id, stable, symbol_raw in zip(
        table["input_gene_id"], table["ensembl_gene_id"], table["gene_symbol"], strict=True
    ):
        symbol = str(symbol_raw).strip()
        if not stable or not symbol or stable in existing or stable in seen:
            continue
        seen.add(stable)
        record: dict[str, object] = dict.fromkeys(map(str, annotation.columns), "")
        record.update(
            {
                "gencode_version": "",
                "ensembl_gene_id_versioned": str(input_id or stable),
                "ensembl_gene_id": stable,
                "gene_symbol": symbol,
                "gene_biotype": "",
                "source_gtf": source,
                "target_id_metadata_source": "output_genes",
                "is_in_output_genes": True,
            }
        )
        rows.append(record)
    if not rows:
        return annotation
    return pd.concat([annotation, pd.DataFrame(rows, columns=annotation.columns)], ignore_index=True)


def mark_output_genes(annotation: pd.DataFrame, ids: frozenset[str] | None) -> pd.DataFrame:
    annotation = annotation.copy()
    annotation["is_in_output_genes"] = annotation["ensembl_gene_id"].isin(ids) if ids else True
    return annotation


def prepare_target(target: TargetSpec, cache: ResourceCache, warnings: list[str]) -> Target:
    output_genes = read_gene_id_list(target.output_genes)
    gene_filter = read_gene_id_list(target.gene_filter, with_symbols=False)
    if output_genes.ids is None and gene_filter.ids is None:
        LOGGER.info("No output_genes or gene_filter supplied for %s; using all annotation genes.", target.name)
    LOGGER.info("%s output_genes ID type: %s", target.name, output_genes.id_type)
    LOGGER.info("%s gene_filter ID type: %s", target.name, gene_filter.id_type)
    annotation, label, path = load_annotation_file(target.name, target.annotation, "annotation", cache)
    primary_ids = frozenset(annotation["ensembl_gene_id"].astype(str))
    fallback, fallback_records = load_metadata_fallbacks(target, output_genes.ids or frozenset(), primary_ids, cache)
    if not fallback.empty:
        annotation = pd.concat([annotation, fallback], ignore_index=True).drop_duplicates()
    helper_ids = frozenset(annotation["ensembl_gene_id"].astype(str))
    if output_genes.ids and (missing := output_genes.ids - helper_ids):
        with_symbols = {
            stable
            for stable, symbol in zip(
                output_genes.table["ensembl_gene_id"], output_genes.table["gene_symbol"], strict=True
            )
            if str(symbol).strip()
        }
        warnings.append(
            f"Target {target.name}: {len(missing)} output_genes lack annotation-helper metadata; "
            f"{len(missing & with_symbols)} have output_genes symbol metadata."
        )
    annotation = add_output_gene_symbols(annotation, output_genes)
    annotation = mark_output_genes(annotation, output_genes.ids or gene_filter.ids)
    return Target(
        name=target.name,
        annotation=annotation.reset_index(drop=True),
        annotation_label=label,
        annotation_path=path,
        output_genes=output_genes,
        gene_filter=gene_filter,
        primary_annotation_ids=primary_ids,
        helper_ids=helper_ids,
        metadata_fallbacks=tuple(fallback_records),
    )


def prepare_targets(targets: Sequence[TargetSpec], cache: ResourceCache, warnings: list[str]) -> list[Target]:
    return [prepare_target(target, cache, warnings) for target in targets]


def write_target_annotation(target: Target, path: Path) -> None:
    write_table(target.annotation, path, ANNOTATION_COLUMNS)
