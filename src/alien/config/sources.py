"""Source specifications: a discriminated union on ``type``."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Annotated, Any, Literal

from pydantic import Field, field_validator, model_validator

from .base import Model, PositiveInt


class SourceBase(Model):
    enabled: bool = True

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
    paths: list[Path] | None = None
    name: str | None = None
    source: str | None = None
    source_tag: str | None = None
    family: str | None = None
    aspect: str | None = None
    collection: str | None = None
    subcollection: str | None = None

    def input_paths(self) -> list[Path]:
        if self.paths is not None:
            return list(self.paths)
        return [self.path] if self.path is not None else []


class SymbolGmtSource(LocalSourceBase):
    """A GMT whose members are gene symbols (local or Enrichr-style)."""

    type: Literal["symbol_gmt", "enrichr_gmt"]
    term_prefix: str | None = None
    db_version: str | None = None
    source_url: str | None = None
    source_license_note: str | None = None


class MsigdbTsvSource(LocalSourceBase):
    """An msigdbr-like long table with term-name and gene-symbol columns."""

    type: Literal["msigdb_tsv"]
    term_name_column: str = "gs_name"
    symbol_column: str = "gene_symbol"
    gene_id_column: str | None = None


class CanonicalTsvSource(LocalSourceBase):
    """A table already in ALIEN's canonical membership layout (needs term_id and gene_symbol)."""

    type: Literal["canonical_tsv"]


class MsigdbCacheSource(SourceBase):
    type: Literal["msigdb_cache"]
    path: Path | None = None
    include_c4_cm: bool = False


class MsigdbRelease(Model):
    zip_url: str | None = None
    zip_md5: str | None = None
    zip_sha256: str | None = None
    zip_name: str | None = None
    summary_rds: str | None = None


class MsigdbRemoteSource(SourceBase):
    type: Literal["msigdb_remote"]
    version: str = "2026.1"
    db_species: str = "HS"
    collection: str | None = None
    subcollection: str | None = None
    family: str | None = None
    aspect: str | None = None
    cache_dir: Path | None = None
    force: bool = False
    timeout_seconds: PositiveInt = 600
    release: MsigdbRelease = MsigdbRelease()
    zip_url: str | None = None
    zip_md5: str | None = None

    @field_validator("version", mode="before")
    @classmethod
    def _version_text(cls, value: Any) -> Any:
        return str(value) if isinstance(value, int | float) else value


class EnrichrLibrary(Model):
    name: str | None = None
    match: str | None = None
    source_tag: str | None = None
    family: str | None = None
    aspect: str | None = None

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
    libraries: list[EnrichrLibrary] = Field(min_length=1)
    source_tag: str | None = None
    family: str | None = None
    aspect: str | None = None
    cache_dir: Path | None = None
    force: bool = False
    timeout_seconds: PositiveInt = 120
    dataset_statistics_endpoint: str | None = None
    library_download_endpoint: str | None = None

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
    plugin: str = Field(min_length=1)
    options: dict[str, Any] = Field(default={})


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
