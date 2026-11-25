"""Typed ALIEN configuration. Defaults live here and nowhere else.

Every model forbids unknown keys, so a typo such as ``jacard_cutoff`` fails with the YAML path of the
offending key instead of being silently ignored.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Annotated, Any, Literal

from pydantic import Field, field_validator, model_validator

from .base import Model, NonNegativeInt, PositiveInt, _reject_removed
from .removed import REMOVED_OPTIONS
from .sources import (
    SOURCE_TYPES,
    CanonicalTsvSource,
    EnrichrLibrary,
    EnrichrRemoteSource,
    LocalSourceBase,
    MsigdbCacheSource,
    MsigdbRelease,
    MsigdbRemoteSource,
    MsigdbTsvSource,
    PluginSource,
    SourceBase,
    SourceSpec,
    SymbolGmtSource,
)
from .targets import (
    DEFAULT_TARGETS,
    AnnotationFile,
    AnnotationSpec,
    GeneIdSpec,
    GtfAttributes,
    MetadataFallbackSpec,
    TargetSpec,
)

DEFAULT_BROAD_DISEASE_REGEX = [
    "^disease$",
    "^diseases$",
    "^cancer$",
    "^neoplasm$",
    "^neoplasms$",
    "^malignant neoplasm$",
    "^malignant neoplasms$",
    "^carcinoma$",
    "^tumor$",
    "^tumour$",
    "^genetic disease$",
    "^genetic disorder$",
]
DEFAULT_SOURCE_PRIORITY = {
    "biology_process_pathway": ["REACTOME", "WIKIPATHWAYS", "KEGG_MEDICUS", "GOBP"],
    "biology_function_location": ["GOMF", "GOCC"],
    "disease_phenotype": [
        "HPO",
        "DISGENET",
        "CLINVAR",
        "GWAS_CATALOG",
        "JENSEN_DISEASES_CURATED",
        "JENSEN_DISEASES_EXPERIMENTAL",
        "JENSEN_DISEASES_TEXTMINING",
    ],
    "cancer_dependency_state": ["C9", "C6", "C4_3CA", "C4_CGN", "C4_CM", "DEPMAP", "CCLE", "NCI60"],
}
DEFAULT_NON_GENE_TOKENS = ["disease", "diseases", "phenotype", "phenotypes", "group", "groups"]
NCBI_GENE_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
NCBI_GENE_HISTORY_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_history.gz"
ENSEMBL_ARCHIVE_ENDPOINT = "https://rest.ensembl.org/archive/id"


class ProjectConfig(Model):
    name: str = "alien"
    source_dir: Path = Path("data/alien_sources")
    outdir: Path = Path("data/alien_gmt")


class OutputsConfig(Model):
    include_symbols: bool = True
    mode: Literal["full", "minimal", "gmt"] = "full"


class TermIdCollisionsConfig(Model):
    action: Literal["error", "warn", "merge"] = "error"


# --------------------------------------------------------------------------------------------- sections


class FilteringConfig(Model):
    min_size_default: NonNegativeInt = 10
    max_size_default: PositiveInt = 500
    max_size_disease: PositiveInt = 500
    max_size_cancer: PositiveInt = 500
    drop_broad_disease_terms: bool = True
    broad_disease_regex: list[str] = Field(default=DEFAULT_BROAD_DISEASE_REGEX)

    @model_validator(mode="before")
    @classmethod
    def _removed(cls, data: Any) -> Any:
        return _reject_removed("filtering", data)

    @field_validator("broad_disease_regex")
    @classmethod
    def _compiles(cls, value: list[str]) -> list[str]:
        for pattern in value:
            re.compile(pattern)
        return value


class RedundancyConfig(Model):
    exact_duplicate_removal: bool = True
    jaccard_cutoff: Annotated[float, Field(ge=0.0, le=1.0)] = 0.85
    apply_per_family: bool = True
    linkage: Literal["single", "representative"] = "single"

    @model_validator(mode="before")
    @classmethod
    def _removed(cls, data: Any) -> Any:
        return _reject_removed("redundancy", data)


class GeneMappingConfig(Model):
    use_hgnc_alias_mapping: bool = True
    preserve_unmapped_symbols_in_symbol_gmt: bool = True
    non_gene_tokens: list[str] = Field(default=DEFAULT_NON_GENE_TOKENS)
    non_gene_token_regex: list[str] = Field(default=[])
    manual_symbol_repairs: dict[str, str] = Field(default={})

    @model_validator(mode="before")
    @classmethod
    def _removed(cls, data: Any) -> Any:
        if isinstance(data, dict) and data.get("manual_symbol_repairs") is None and "manual_symbol_repairs" in data:
            data = {**data, "manual_symbol_repairs": {}}
        return _reject_removed("gene_mapping", data)

    @field_validator("non_gene_token_regex")
    @classmethod
    def _compiles(cls, value: list[str]) -> list[str]:
        for pattern in value:
            re.compile(pattern)
        return value


class NcbiGeneConfig(Model):
    enabled: bool = True
    allowed_source_tags: list[str] = Field(default=["CCLE", "NCI60", "DEPMAP"])
    gene_info_url: str = NCBI_GENE_INFO_URL
    gene_history_url: str = NCBI_GENE_HISTORY_URL
    use_gene_info_synonyms: bool = True


class EnsemblArchiveConfig(Model):
    enabled: bool = True
    endpoint: str = ENSEMBL_ARCHIVE_ENDPOINT
    batch_size: PositiveInt = 50
    connect_timeout_seconds: PositiveInt = 30
    timeout_seconds: PositiveInt = 120
    retry_attempts: PositiveInt = 3
    retry_backoff_seconds: Annotated[float, Field(ge=0.0)] = 5
    unresolved_print_limit: NonNegativeInt = 20

    @model_validator(mode="before")
    @classmethod
    def _removed(cls, data: Any) -> Any:
        return _reject_removed("ensembl_archive", data)


class HgncConfig(Model):
    """Where the HGNC complete set comes from. Pin ``archive_date`` (YYYY-MM-DD) for a dated monthly archive."""

    url: str | None = None
    archive_date: Annotated[str, Field(pattern=r"^[0-9]{4}-[0-9]{2}-[0-9]{2}$")] | None = None


class RuntimeConfig(Model):
    workers: PositiveInt = 1


class AlienConfig(Model):
    """The complete, validated configuration for one build."""

    project: ProjectConfig = ProjectConfig()
    outputs: OutputsConfig = OutputsConfig()
    term_id_collisions: TermIdCollisionsConfig = TermIdCollisionsConfig()
    sources: list[SourceSpec]
    targets: list[TargetSpec] = Field(default_factory=lambda: [TargetSpec.model_validate(t) for t in DEFAULT_TARGETS])
    filtering: FilteringConfig = FilteringConfig()
    redundancy: RedundancyConfig = RedundancyConfig()
    source_priority: dict[str, list[str]] = Field(default=DEFAULT_SOURCE_PRIORITY)
    gene_mapping: GeneMappingConfig = GeneMappingConfig()
    hgnc: HgncConfig = HgncConfig()
    ncbi_gene: NcbiGeneConfig = NcbiGeneConfig()
    ensembl_archive: EnsemblArchiveConfig = EnsemblArchiveConfig()
    runtime: RuntimeConfig = RuntimeConfig()

    @model_validator(mode="before")
    @classmethod
    def _shorthands(cls, data: Any) -> Any:
        if isinstance(data, dict) and isinstance(data.get("term_id_collisions"), str):
            data = {**data, "term_id_collisions": {"action": data["term_id_collisions"]}}
        return data

    @field_validator("sources")
    @classmethod
    def _has_sources(cls, value: list[Any]) -> list[Any]:
        if not value:
            raise ValueError("a config must define at least one source")
        return value

    @model_validator(mode="after")
    def _namespaces(self) -> AlienConfig:
        if not self.outputs.include_symbols and not self.targets:
            raise ValueError("enable outputs.include_symbols or define at least one target namespace")
        names = [target.name for target in self.targets]
        duplicates = sorted({name for name in names if names.count(name) > 1})
        if duplicates:
            raise ValueError(f"target names must be unique; repeated: {', '.join(duplicates)}")
        if "symbols" in names and self.outputs.include_symbols:
            raise ValueError("target name 'symbols' is reserved for the symbols namespace")
        return self

    def namespaces(self) -> list[str]:
        return (["symbols"] if self.outputs.include_symbols else []) + [target.name for target in self.targets]


def default_layer() -> dict[str, Any]:
    """Defaults as a plain mapping, the bottom layer that presets, files and overrides merge onto."""
    sections = {
        "project": ProjectConfig(),
        "outputs": OutputsConfig(),
        "term_id_collisions": TermIdCollisionsConfig(),
        "filtering": FilteringConfig(),
        "redundancy": RedundancyConfig(),
        "gene_mapping": GeneMappingConfig(),
        "hgnc": HgncConfig(),
        "ncbi_gene": NcbiGeneConfig(),
        "ensembl_archive": EnsemblArchiveConfig(),
        "runtime": RuntimeConfig(),
    }
    layer: dict[str, Any] = {name: model.model_dump(mode="json") for name, model in sections.items()}
    layer["source_priority"] = {family: list(tags) for family, tags in DEFAULT_SOURCE_PRIORITY.items()}
    layer["sources"] = []
    layer["targets"] = [dict(target) for target in DEFAULT_TARGETS]
    return layer


__all__ = [
    "DEFAULT_TARGETS",
    "REMOVED_OPTIONS",
    "SOURCE_TYPES",
    "AlienConfig",
    "AnnotationFile",
    "AnnotationSpec",
    "CanonicalTsvSource",
    "EnrichrLibrary",
    "EnrichrRemoteSource",
    "EnsemblArchiveConfig",
    "FilteringConfig",
    "GeneIdSpec",
    "GeneMappingConfig",
    "GtfAttributes",
    "HgncConfig",
    "LocalSourceBase",
    "MetadataFallbackSpec",
    "Model",
    "MsigdbCacheSource",
    "MsigdbRelease",
    "MsigdbRemoteSource",
    "MsigdbTsvSource",
    "NcbiGeneConfig",
    "OutputsConfig",
    "PluginSource",
    "ProjectConfig",
    "RedundancyConfig",
    "RuntimeConfig",
    "SourceBase",
    "SourceSpec",
    "SymbolGmtSource",
    "TargetSpec",
    "TermIdCollisionsConfig",
    "default_layer",
]
