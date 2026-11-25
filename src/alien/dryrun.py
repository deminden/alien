"""Dry-run planning: validate inputs and describe what a build would read, download and write."""

from __future__ import annotations

from pathlib import Path

from .config.schema import (
    AlienConfig,
    EnrichrRemoteSource,
    LocalSourceBase,
    MsigdbCacheSource,
    MsigdbRemoteSource,
    SymbolGmtSource,
    TargetSpec,
)
from .errors import ConfigError
from .model import sanitize_term_id
from .outputs import MINIMAL_FILES
from .resources.annotation import annotation_location
from .resources.hgnc import hgnc_path, hgnc_urls
from .sources.msigdb import release_info


def _status(path: Path) -> str:
    return "cached" if path.exists() else "missing, would download"


def _source_inputs(config: AlienConfig) -> list[tuple[str, Path]]:
    inputs: list[tuple[str, Path]] = []
    for index, spec in enumerate(config.sources):
        if isinstance(spec, LocalSourceBase) and spec.enabled:
            paths = spec.input_paths()
            if not paths and isinstance(spec, SymbolGmtSource):
                paths = [config.project.source_dir / (spec.name or "symbol_gmt")]
            inputs += [(f"sources[{index}]", path) for path in paths]
        elif isinstance(spec, MsigdbCacheSource) and spec.enabled:
            inputs.append((f"sources[{index}].path", spec.path or config.project.source_dir / "msigdb"))
    return inputs


def _target_inputs(index: int, target: TargetSpec) -> list[tuple[str, Path]]:
    prefix = f"targets[{index}]"
    candidates = [(f"{prefix}.annotation.path", target.annotation.path)]
    candidates += [
        (f"{prefix}.annotation.metadata_fallbacks[{i}].path", fallback.path)
        for i, fallback in enumerate(target.annotation.metadata_fallbacks)
    ]
    candidates += [
        (f"{prefix}.{key}.path", ids.path if ids is not None else None)
        for key, ids in (("output_genes", target.output_genes), ("gene_filter", target.gene_filter))
    ]
    return [(label, path) for label, path in candidates if path is not None]


def local_inputs(config: AlienConfig) -> list[tuple[str, Path]]:
    """Every local file the config refers to, labelled with its YAML path."""
    inputs = _source_inputs(config)
    for index, target in enumerate(config.targets):
        inputs += _target_inputs(index, target)
    return inputs


def _msigdb_line(spec: MsigdbRemoteSource, source_dir: Path) -> str:
    release = release_info(spec)
    cache_dir = spec.cache_dir or source_dir / "msigdb_remote"
    path = cache_dir / release["zip_name"]
    extracted = cache_dir / Path(release["zip_name"]).stem / release["summary_rds"]
    state = "cached" if path.exists() or extracted.exists() else "missing, would download"
    return f"MSigDB {release['zip_url']} -> {path} ({state})"


def _enrichr_lines(spec: EnrichrRemoteSource, source_dir: Path) -> list[str]:
    cache_dir = spec.cache_dir or source_dir / "enrichr"
    index = cache_dir / "datasetStatistics.json"
    lines = [f"Enrichr library index -> {index} ({_status(index)})"]
    for library in spec.libraries:
        if library.name and not library.match:
            path = cache_dir / f"{sanitize_term_id(library.name)}.gmt"
            lines.append(f"Enrichr library {library.name} -> {path} ({_status(path)})")
        else:
            lines.append(f"Enrichr library matching {library.match!r} (resolved at build time)")
    return lines


def _remote_sources(config: AlienConfig, source_dir: Path) -> list[str]:
    lines: list[str] = []
    for spec in (spec for spec in config.sources if spec.enabled):
        if isinstance(spec, MsigdbRemoteSource):
            lines.append(_msigdb_line(spec, source_dir))
        elif isinstance(spec, EnrichrRemoteSource):
            lines += _enrichr_lines(spec, source_dir)
    return lines


def _annotation_resource(target: TargetSpec, source_dir: Path) -> list[str]:
    files = [("annotation", target.annotation)] + [
        (f"metadata_fallback_{i}", fallback) for i, fallback in enumerate(target.annotation.metadata_fallbacks, start=1)
    ]
    lines = []
    for role, spec in files:
        path, url, label, _ = annotation_location(target.name, spec, role, source_dir)
        if url:
            lines.append(f"Target {target.name} {role} {label} {url} -> {path} ({_status(path)})")
    return lines


def remote_resources(config: AlienConfig) -> list[str]:
    source_dir = config.project.source_dir
    hgnc = hgnc_path(source_dir, config.hgnc)
    lines = [f"HGNC {hgnc_urls(config.hgnc)[0]} -> {hgnc} ({_status(hgnc)})"]
    if config.ncbi_gene.enabled:
        for url, name in (
            (config.ncbi_gene.gene_info_url, "Homo_sapiens.gene_info.gz"),
            (config.ncbi_gene.gene_history_url, "gene_history.gz"),
        ):
            lines.append(
                f"NCBI Gene {url} -> {source_dir / 'ncbi_gene' / name} ({_status(source_dir / 'ncbi_gene' / name)})"
            )
    lines += _remote_sources(config, source_dir)
    for target in config.targets:
        lines += _annotation_resource(target, source_dir)
    if config.ensembl_archive.enabled:
        path = source_dir / "ensembl_archive" / "archive_id_cache.json"
        lines.append(f"Ensembl archive lookups cached in {path} ({'present' if path.exists() else 'empty'})")
    return lines


def planned_outputs(config: AlienConfig) -> list[str]:
    files = [f"gmt/{name}.gmt" for name in config.namespaces()]
    if config.outputs.mode == "minimal":
        files += list(MINIMAL_FILES)
    elif config.outputs.mode == "full":
        files += ["metadata/...", "qc/..."]
    return [str(config.project.outdir / name) for name in [*files, "MANIFEST.json"]]


def dry_run_plan(config: AlienConfig) -> list[str]:
    """Raise ConfigError for missing local inputs, else return a human-readable build plan."""
    inputs = local_inputs(config)
    missing = [f"{label}: {path}" for label, path in inputs if not path.exists()]
    if missing:
        raise ConfigError("Configured local inputs do not exist:\n  " + "\n  ".join(missing))
    mode, workers = config.outputs.mode, config.runtime.workers
    plan = [f"Dry run: output mode {mode}, {workers} worker process(es); nothing is downloaded or written."]
    plan += ["Local inputs:"] + [f"  {label}: {path} (exists)" for label, path in inputs]
    plan += ["Remote resources:"] + [f"  {line}" for line in remote_resources(config)]
    plan += ["Outputs:"] + [f"  {line}" for line in planned_outputs(config)]
    return plan
