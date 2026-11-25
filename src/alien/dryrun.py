"""Dry-run planning: validate inputs and describe what a build would read, download and write."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .outputs import MINIMAL_FILES
from .sources import (
    HGNC_URL,
    NCBI_GENE_HISTORY_URL,
    NCBI_GENE_INFO_URL,
    _msigdb_release_info,
    _source_paths,
    gencode_annotation_url,
    normalize_gencode_version,
)
from .utils import sanitize_id

LOCAL_SOURCE_TYPES = {"canonical_tsv", "msigdb_tsv", "symbol_gmt", "enrichr_gmt"}


class DryRunError(FileNotFoundError):
    """A dry run found local inputs that do not exist."""


def _status(path: Path) -> str:
    return "cached" if path.exists() else "missing, would download"


def _local_inputs(cfg: dict[str, Any], source_dir: Path) -> list[tuple[str, Path]]:
    inputs: list[tuple[str, Path]] = []
    for index, spec in enumerate(cfg.get("sources", [])):
        if not spec.get("enabled", True):
            continue
        if spec.get("type") in LOCAL_SOURCE_TYPES:
            inputs.extend((f"sources[{index}]", path) for path in _source_paths(spec, source_dir))
        if spec.get("type") == "msigdb_cache":
            inputs.append((f"sources[{index}].path", Path(spec.get("path") or source_dir / "msigdb")))
    for index, target in enumerate(cfg.get("targets", [])):
        annotation = target.get("annotation") or {}
        local = annotation.get("path") or target.get("annotation_gtf")
        if local:
            inputs.append((f"targets[{index}].annotation.path", Path(str(local))))
        fallbacks = annotation.get("metadata_fallbacks") or []
        for fb_index, fallback in enumerate(fallbacks if isinstance(fallbacks, list) else [fallbacks]):
            if isinstance(fallback, dict) and fallback.get("path"):
                inputs.append(
                    (f"targets[{index}].annotation.metadata_fallbacks[{fb_index}].path", Path(fallback["path"]))
                )
        for key in ("output_genes", "gene_filter"):
            value = target.get(key)
            path = value.get("path") or value.get("file") if isinstance(value, dict) else value
            if isinstance(path, str) and path.strip() and path.strip() != "null":
                inputs.append((f"targets[{index}].{key}", Path(path)))
    return inputs


def _mapping_resources(cfg: dict[str, Any], source_dir: Path) -> list[str]:
    hgnc = source_dir / "hgnc" / "hgnc_complete_set.txt"
    lines = [f"HGNC {HGNC_URL} -> {hgnc} ({_status(hgnc)})"]
    ncbi = cfg.get("ncbi_gene", {})
    if ncbi.get("enabled", False):
        for url, name in [
            (ncbi.get("gene_info_url", NCBI_GENE_INFO_URL), "Homo_sapiens.gene_info.gz"),
            (ncbi.get("gene_history_url", NCBI_GENE_HISTORY_URL), "gene_history.gz"),
        ]:
            path = source_dir / "ncbi_gene" / name
            lines.append(f"NCBI Gene {url} -> {path} ({_status(path)})")
    if cfg.get("ensembl_archive", {}).get("enabled", True):
        path = source_dir / "ensembl_archive" / "archive_id_cache.json"
        lines.append(f"Ensembl archive lookups cached in {path} ({'present' if path.exists() else 'empty'})")
    return lines


def _msigdb_resource(spec: dict[str, Any], source_dir: Path) -> list[str]:
    release = _msigdb_release_info(spec)
    cache_dir = Path(spec.get("cache_dir") or source_dir / "msigdb_remote")
    path = cache_dir / release["zip_name"]
    extracted = cache_dir / Path(release["zip_name"]).stem / release["summary_rds"]
    state = "cached" if path.exists() or extracted.exists() else "missing, would download"
    return [f"MSigDB {release['zip_url']} -> {path} ({state})"]


def _enrichr_resources(spec: dict[str, Any], source_dir: Path) -> list[str]:
    cache_dir = Path(spec.get("cache_dir") or source_dir / "enrichr")
    metadata = cache_dir / "datasetStatistics.json"
    lines = [f"Enrichr library index -> {metadata} ({_status(metadata)})"]
    for library in spec.get("libraries") or []:
        exact = isinstance(library, str) or bool(library.get("name"))
        name = library if isinstance(library, str) else library.get("name") or library.get("match")
        path = cache_dir / f"{sanitize_id(str(name))}.gmt"
        lines.append(f"Enrichr library {name} -> {path} ({_status(path) if exact else 'resolved at build time'})")
    return lines


def _annotation_resource(target: dict[str, Any], source_dir: Path) -> list[str]:
    annotation = target.get("annotation") or {}
    if annotation.get("path") or target.get("annotation_gtf"):
        return []
    name = sanitize_id(str(target.get("name", "")))
    if annotation.get("url"):
        file_name = annotation.get("file_name") or Path(str(annotation["url"]).split("?", 1)[0]).name
        path = source_dir / "targets" / name / "annotation" / str(file_name)
        return [f"Target {name} annotation {annotation['url']} -> {path} ({_status(path)})"]
    version = annotation.get("version") or target.get("gencode_version")
    if not version:
        return []
    version = normalize_gencode_version(version)
    path = source_dir / "gencode" / f"gencode.v{version}.annotation.gtf.gz"
    return [f"Target {name} GENCODE {gencode_annotation_url(version)} -> {path} ({_status(path)})"]


def _remote_resources(cfg: dict[str, Any], source_dir: Path) -> list[str]:
    lines = _mapping_resources(cfg, source_dir)
    for spec in cfg.get("sources", []):
        if not spec.get("enabled", True):
            continue
        if spec.get("type") == "msigdb_remote":
            lines += _msigdb_resource(spec, source_dir)
        if spec.get("type") == "enrichr_remote":
            lines += _enrichr_resources(spec, source_dir)
    for target in cfg.get("targets", []):
        lines += _annotation_resource(target, source_dir)
    return lines


def _planned_outputs(cfg: dict[str, Any], outdir: Path, mode: str) -> list[str]:
    namespaces = (["symbols"] if cfg.get("outputs", {}).get("include_symbols", True) else []) + [
        str(target["name"]) for target in cfg.get("targets", [])
    ]
    files = [f"gmt/{name}.gmt" for name in namespaces]
    if mode == "minimal":
        files += list(MINIMAL_FILES)
    elif mode == "full":
        files += ["metadata/...", "qc/..."]
    return [f"{outdir / name}" for name in [*files, "MANIFEST.json"]]


def dry_run_plan(cfg: dict[str, Any], source_dir: Path, outdir: Path, mode: str, workers: int) -> list[str]:
    """Raise DryRunError for missing local inputs, else return a human-readable build plan."""
    inputs = _local_inputs(cfg, source_dir)
    missing = [f"{label}: {path}" for label, path in inputs if not path.exists()]
    if missing:
        raise DryRunError("Configured local inputs do not exist:\n  " + "\n  ".join(missing))
    plan = [f"Dry run: output mode {mode}, {workers} worker process(es); nothing is downloaded or written."]
    plan += ["Local inputs:"] + [f"  {label}: {path} (exists)" for label, path in inputs]
    plan += ["Remote resources:"] + [f"  {line}" for line in _remote_resources(cfg, source_dir)]
    plan += ["Outputs:"] + [f"  {line}" for line in _planned_outputs(cfg, outdir, mode)]
    return plan
