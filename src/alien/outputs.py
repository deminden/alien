"""The output-directory contract: staging, output modes, MANIFEST.json and verification."""

from __future__ import annotations

import hashlib
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

from .io.atomic import atomic_dir_swap, atomic_write_json, unique_sibling

MANIFEST_NAME: Final = "MANIFEST.json"
FAILED_DIR: Final = "failed"
OUTPUT_MODES: Final = ("full", "minimal", "gmt")
LEGACY_TOP_LEVEL: Final = frozenset({"gmt", "metadata", "qc", FAILED_DIR})
MINIMAL_FILES: Final = (
    "metadata/effective_config.yml",
    "metadata/source_manifest.tsv",
    "metadata/source_provenance.json",
    "metadata/term_manifest.tsv.gz",
    "metadata/target_metadata_fallbacks.tsv",
    "metadata/target_output_genes.tsv.gz",
    "metadata/target_gene_filter.tsv.gz",
    "qc/target_namespace_summary.tsv",
    "qc/mapping_summary.tsv",
    "qc/redundancy_summary.tsv",
    "qc/warnings.txt",
)


class OutputDirectoryError(RuntimeError):
    """The output directory cannot be used safely."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_alien_output_dir(path: Path) -> bool:
    """True when ``path`` is absent, empty, or looks like an earlier ALIEN build."""
    if not path.exists():
        return True
    if not path.is_dir():
        return False
    entries = {entry.name for entry in path.iterdir()}
    return not entries or MANIFEST_NAME in entries or entries <= LEGACY_TOP_LEVEL


def selected_files(root: Path, mode: str) -> list[str]:
    """Relative paths kept for an output mode."""
    gmts = sorted(path.relative_to(root).as_posix() for path in (root / "gmt").glob("*.gmt"))
    if mode == "gmt":
        return gmts
    if mode == "minimal":
        return gmts + [name for name in MINIMAL_FILES if (root / name).exists()]
    return sorted(path.relative_to(root).as_posix() for path in root.rglob("*") if path.is_file())


def prune_to_mode(root: Path, mode: str) -> None:
    keep = set(selected_files(root, mode))
    for path in sorted(root.rglob("*"), reverse=True):
        rel = path.relative_to(root).as_posix()
        if path.is_file() and rel not in keep:
            path.unlink()
        elif path.is_dir() and not any(path.iterdir()):
            path.rmdir()


def write_manifest(root: Path, *, mode: str, alien_version: str, schema_version: int) -> Path:
    files = {
        rel: {"sha256": sha256_file(root / rel), "size": (root / rel).stat().st_size}
        for rel in selected_files(root, "full")
    }
    manifest = {
        "tool": "ALIEN",
        "alien_version": alien_version,
        "output_schema_version": schema_version,
        "output_mode": mode,
        "files": files,
    }
    path = root / MANIFEST_NAME
    atomic_write_json(path, manifest)
    return path


@dataclass
class StagedOutput:
    """A build that writes into a hidden sibling directory and is swapped into place on success."""

    outdir: Path
    force: bool = False
    staging: Path = field(init=False)

    def __post_init__(self) -> None:
        if not self.force and not is_alien_output_dir(self.outdir):
            raise OutputDirectoryError(
                f"Refusing to overwrite {self.outdir}: it is not empty and has no ALIEN {MANIFEST_NAME}. "
                "Choose another --outdir or pass --force."
            )
        self.staging = unique_sibling(self.outdir, "staging")

    @property
    def failed_dir(self) -> Path:
        return self.outdir / FAILED_DIR

    def commit(self) -> None:
        if self.outdir.exists() and not self.outdir.is_dir():
            self.outdir.unlink()
        atomic_dir_swap(self.staging, self.outdir)

    def fail(self) -> Path:
        """Keep the staged diagnostics under ``outdir/failed`` and leave any previous build intact."""
        self.outdir.mkdir(parents=True, exist_ok=True)
        if self.failed_dir.exists():
            shutil.rmtree(self.failed_dir)
        self.staging.replace(self.failed_dir)
        return self.failed_dir

    def discard(self) -> None:
        shutil.rmtree(self.staging, ignore_errors=True)
