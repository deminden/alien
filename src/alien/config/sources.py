"""Source specifications: a discriminated union on ``type``."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Annotated, Any, Literal

from pydantic import Field, field_validator, model_validator

from .base import Model, PositiveInt


class SourceBase(Model):
    enabled: bool = True
    """Include this source; `false` keeps it in the file but skips it."""

    @model_validator(mode="before")
    @classmethod
    def _no_optional(cls, data: Any) -> Any:
        if isinstance(data, dict) and "optional" in data:
            raise ValueError(
                "source option 'optional' is not supported; configured sources are required. "
                "Use enabled: false to exclude a source."
            )
        return data


class LocalSourceBase(SourceBase):
    path: Path | None = None
    """Input file."""
    paths: list[Path] | None = None
    """Several input files (instead of `path`)."""
    name: str | None = None
    """Source name; also the directory under `source_dir` scanned for GMTs when no path is given (`symbol_gmt`)."""
    source: str | None = None
    """Source label written to the manifest (`msigdb_tsv` always records MSigDB)."""
    source_tag: str | None = None
    """Short tag used in term IDs and `source_priority` (default: derived from `source` or `name`)."""
    family: str | None = None
    """Term family, which groups terms for redundancy filtering (default `biology_process_pathway`)."""
    aspect: str | None = None
    """Finer label within the family (default `library`)."""
    collection: str | None = None
    """Collection label (GMT and canonical sources default to the file name)."""
    subcollection: str | None = None
    """Subcollection label."""

    def input_paths(self) -> list[Path]:
        if self.paths is not None:
            return list(self.paths)
        return [self.path] if self.path is not None else []


class SymbolGmtSource(LocalSourceBase):
    """A GMT whose members are gene symbols (local or Enrichr-style)."""

    type: Literal["symbol_gmt", "enrichr_gmt"]
    """`symbol_gmt` or `enrichr_gmt`."""
    term_prefix: str | None = None
    """Prefix for term IDs (default: the source tag)."""
    db_version: str | None = None
    """Upstream version recorded in the manifest."""
    source_url: str | None = None
    """Upstream URL recorded in the manifest."""
    source_license_note: str | None = None
    """License note recorded in the manifest."""


class MsigdbTsvSource(LocalSourceBase):
    """An msigdbr-like long table with term-name and gene-symbol columns."""

    type: Literal["msigdb_tsv"]
    """`msigdb_tsv`."""
    term_name_column: str = "gs_name"
    """Column holding the term name."""
    symbol_column: str = "gene_symbol"
    """Column holding the gene symbol."""
    gene_id_column: str | None = None
    """Optional column holding a source Ensembl gene ID."""


class CanonicalTsvSource(LocalSourceBase):
    """A table already in ALIEN's canonical membership layout (needs term_id and gene_symbol)."""

    type: Literal["canonical_tsv"]
    """`canonical_tsv`."""


class MsigdbCacheSource(SourceBase):
    type: Literal["msigdb_cache"]
    """`msigdb_cache`."""
    path: Path | None = None
    """Directory with `msigdbr_<TAG>.tsv.gz` files."""
    include_c4_cm: bool = False
    """Include the C4 cancer-module collection."""


class MsigdbRelease(Model):
    zip_url: str | None = None
    """Archive URL; when set, no built-in release fields are inherited."""
    zip_md5: str | None = None
    """Expected archive MD5."""
    zip_sha256: str | None = None
    """Expected archive sha256."""
    zip_name: str | None = None
    """Local file name for the archive (default: from the URL)."""
    summary_rds: str | None = None
    """Summary table inside the archive (default: `<zip name>.summary.rds`)."""


class MsigdbRemoteSource(SourceBase):
    type: Literal["msigdb_remote"]
    """`msigdb_remote`."""
    version: str = "2026.1"
    """Bundled MSigDB release to download, or a label for a custom `release`."""
    db_species: str = "HS"
    """Database species: `HS` (human) or `MM` (mouse)."""
    collection: str | None = None
    """MSigDB collection, such as `C2` or `C5`."""
    subcollection: str | None = None
    """Subcollection, full (`GO:BP`) or short (`BP`)."""
    family: str | None = None
    """Override the term family."""
    aspect: str | None = None
    """Override the aspect."""
    cache_dir: Path | None = None
    """Cache directory (default `source_dir/msigdb_remote`)."""
    force: bool = False
    """Re-download and re-extract the archive."""
    timeout_seconds: PositiveInt = 600
    """Download timeout."""
    release: MsigdbRelease = MsigdbRelease()
    """Custom release archive instead of a bundled version."""
    zip_url: str | None = None
    """Shorthand for `release.zip_url`, used with `zip_md5` when no bundled release matches `version`."""
    zip_md5: str | None = None
    """Shorthand for `release.zip_md5`."""

    @field_validator("version", mode="before")
    @classmethod
    def _version_text(cls, value: Any) -> Any:
        return str(value) if isinstance(value, int | float) else value


class EnrichrLibrary(Model):
    name: str | None = None
    """Exact Enrichr library name, tried first."""
    match: str | None = None
    """Regex selecting the latest matching library when the exact name is absent."""
    source_tag: str | None = None
    """Tag for this library's terms."""
    family: str | None = None
    """Term family."""
    aspect: str | None = None
    """Aspect label."""

    @model_validator(mode="before")
    @classmethod
    def _from_string(cls, data: Any) -> Any:
        if isinstance(data, dict) and "optional" in data:
            raise ValueError("library option 'optional' is not supported; configured libraries are required.")
        return {"name": data} if isinstance(data, str) else data

    @model_validator(mode="after")
    def _name_or_match(self) -> EnrichrLibrary:
        if not self.name and not self.match:
            raise ValueError("every Enrichr library must define name or match")
        if self.match:
            re.compile(self.match)
        return self


class EnrichrRemoteSource(SourceBase):
    type: Literal["enrichr_remote"]
    """`enrichr_remote`."""
    libraries: list[EnrichrLibrary] = Field(min_length=1)
    """Libraries to download; a string is shorthand for `{name: ...}`."""
    source_tag: str | None = None
    """Default tag for all libraries."""
    family: str | None = None
    """Default family for all libraries."""
    aspect: str | None = None
    """Default aspect for all libraries."""
    cache_dir: Path | None = None
    """Cache directory (default `source_dir/enrichr`)."""
    force: bool = False
    """Refresh the library list and downloads."""
    timeout_seconds: PositiveInt = 120
    """Download timeout."""
    dataset_statistics_endpoint: str | None = None
    """Override the Enrichr library-list URL."""
    library_download_endpoint: str | None = None
    """Override the Enrichr download URL template (`{library}` is substituted)."""

    @model_validator(mode="before")
    @classmethod
    def _single_library(cls, data: Any) -> Any:
        if not isinstance(data, dict) or "libraries" in data:
            return data
        data = dict(data)
        if "library" in data:
            data["libraries"] = [{"name": data.pop("library")}]
        elif "name" in data or "match" in data:
            data["libraries"] = [{key: data.pop(key) for key in ("name", "match") if key in data}]
        return data


class PluginSource(SourceBase):
    """A source provided by an installed package through the ``alien.sources`` entry-point group."""

    type: Literal["plugin"]
    """`plugin`."""
    plugin: str = Field(min_length=1)
    """Entry-point name in the `alien.sources` group."""
    options: dict[str, Any] = Field(default={})
    """Passed to the plugin unchanged."""


SourceSpec = Annotated[
    SymbolGmtSource
    | MsigdbTsvSource
    | CanonicalTsvSource
    | MsigdbCacheSource
    | MsigdbRemoteSource
    | EnrichrRemoteSource
    | PluginSource,
    Field(discriminator="type"),
]
SOURCE_TYPES = (
    "symbol_gmt",
    "enrichr_gmt",
    "msigdb_tsv",
    "canonical_tsv",
    "msigdb_cache",
    "msigdb_remote",
    "enrichr_remote",
    "plugin",
)
