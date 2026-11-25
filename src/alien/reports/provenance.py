"""source_provenance.json, assembled from what the build actually read: nothing is hard-coded."""

from __future__ import annotations

import hashlib
import os
import platform
import subprocess
from collections.abc import Sequence
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from typing import Any, Final

from ..cache import ResourceCache, file_sha256
from ..config.schema import AlienConfig, LocalSourceBase
from ..model import Target
from ..schema import OUTPUT_SCHEMA_VERSION
from ..sources.local import source_paths
from ..version import __version__

RUNTIME_PACKAGES: Final = ("pandas", "pyarrow", "PyYAML", "rdata", "requests", "urllib3", "pydantic", "filelock")


def build_timestamp() -> str:
    """UTC now, or SOURCE_DATE_EPOCH when set (reproducible builds)."""
    epoch = os.environ.get("SOURCE_DATE_EPOCH", "").strip()
    moment = datetime.fromtimestamp(int(epoch), UTC) if epoch.isdigit() else datetime.now(UTC)
    return moment.isoformat(timespec="seconds")


def python_packages() -> dict[str, str]:
    packages: dict[str, str] = {}
    for package in RUNTIME_PACKAGES:
        try:
            packages[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            packages[package] = "not installed"
    return packages


def git_describe() -> str:
    root = Path(__file__).resolve().parents[3]
    if not (root / ".git").exists():
        return ""
    try:
        result = subprocess.run(
            ["git", "describe", "--always", "--dirty", "--tags"],
            cwd=root,
            capture_output=True,
            text=True,
            check=False,
            timeout=5,
        )
    except (OSError, subprocess.SubprocessError):
        return ""
    return result.stdout.strip()


def _file_record(role: str, path: Path) -> dict[str, Any]:
    record: dict[str, Any] = {"role": role, "path": str(path)}
    if path.is_file():
        record.update({"sha256": file_sha256(path), "size": path.stat().st_size})
    return record


def input_files(config: AlienConfig, targets: Sequence[Target]) -> list[dict[str, Any]]:
    """Every local (non-downloaded) input with its sha256 and size."""
    records: list[dict[str, Any]] = []
    for index, spec in enumerate(config.sources):
        if spec.enabled and isinstance(spec, LocalSourceBase):
            records += [
                _file_record(f"sources[{index}]", path) for path in source_paths(spec, config.project.source_dir)
            ]
    for target in targets:
        if not target.annotation_path.resolve().is_relative_to(config.project.source_dir.resolve()):
            records.append(_file_record(f"target:{target.name}:annotation", target.annotation_path))
        for fallback in target.metadata_fallbacks:
            role = f"target:{target.name}:metadata_fallback_{fallback['fallback_index']}"
            records.append(_file_record(role, Path(str(fallback["annotation_path"]))))
        for kind, gene_list in (("output_genes", target.output_genes), ("gene_filter", target.gene_filter)):
            if gene_list.source_path is not None:
                records.append(_file_record(f"target:{target.name}:{kind}", gene_list.source_path))
    return records


def _target_record(target: Target) -> dict[str, Any]:
    return {
        "name": target.name,
        "type": "ensembl_gtf",
        "annotation": target.annotation_label,
        "annotation_path": str(target.annotation_path),
        "metadata_fallbacks": list(target.metadata_fallbacks),
        "output_genes_path": str(target.output_genes.source_path or ""),
        "output_genes_id_column": target.output_genes.id_column or "",
        "output_genes_symbol_column": target.output_genes.symbol_column or "",
        "gene_filter_path": str(target.gene_filter.source_path or ""),
        "gene_filter_column": target.gene_filter.id_column or "",
    }


def build_provenance(
    config: AlienConfig,
    targets: Sequence[Target],
    cache: ResourceCache,
    warnings: Sequence[str],
    effective_config: str,
) -> dict[str, Any]:
    return {
        "tool": "ALIEN",
        "alien_version": __version__,
        "git_describe": git_describe(),
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "created_at": build_timestamp(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "packages": python_packages(),
        "config_sha256": hashlib.sha256(effective_config.encode("utf-8")).hexdigest(),
        "configured_sources": [spec.model_dump(mode="json", exclude_none=True) for spec in config.sources],
        "configured_targets": [_target_record(target) for target in targets],
        "input_files": input_files(config, targets),
        "resources": dict(sorted(cache.used.items())),
        "cache_mode": "frozen" if cache.frozen else "offline" if cache.offline else "online",
        "redundancy": config.redundancy.model_dump(mode="json"),
        "workers": config.runtime.workers,
        "warnings": list(warnings),
    }
