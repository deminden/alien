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
    """Project label recorded in provenance."""
    source_dir: Path = Path("data/alien_sources")
    """Cache directory for downloaded and derived resources (relative to the config file)."""
    outdir: Path = Path("data/alien_gmt")
    """Output directory; a build is staged beside it and swapped in atomically."""


class OutputsConfig(Model):
    include_symbols: bool = True
    """Also write `gmt/symbols.gmt` with HGNC-resolved symbols."""
    mode: Literal["full", "minimal", "gmt"] = "full"
    """`full` writes every audit table, `minimal` GMTs plus reproducibility metadata, `gmt` only GMTs and
    the manifest."""


class TermIdCollisionsConfig(Model):
    action: Literal["error", "warn", "merge"] = "error"
    """What to do when one `term_id` carries conflicting source metadata: `error` stops, `warn` and `merge` combine the
    memberships and write an audit table."""


# --------------------------------------------------------------------------------------------- sections


class FilteringConfig(Model):
    min_size_default: NonNegativeInt = 10
    """Remove terms with fewer mapped genes than this, per namespace."""
    max_size_default: PositiveInt = 500
    """Remove terms with more mapped genes than this (families without a specific limit)."""
    max_size_disease: PositiveInt = 500
    """Upper size limit for the `disease_phenotype` family."""
    max_size_cancer: PositiveInt = 500
    """Upper size limit for the `cancer_dependency_state` family."""
    drop_broad_disease_terms: bool = True
    """Remove disease terms whose name matches `broad_disease_regex` (such as \"cancer\" or \"disease\")."""
    broad_disease_regex: list[str] = Field(default=DEFAULT_BROAD_DISEASE_REGEX)
    """Case-insensitive patterns for broad disease term names; `extra_broad_disease_regex` appends to the defaults."""

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
    """Collapse terms with identical mapped gene sets before the Jaccard pass."""
    jaccard_cutoff: Annotated[float, Field(ge=0.0, le=1.0)] = 0.85
    """Terms with Jaccard similarity strictly above this are redundant."""
    apply_per_family: bool = True
    """Compare terms only within the same family; `false` compares across all families in a namespace."""
    linkage: Literal["single", "representative"] = "single"
    """`single` clusters similar terms transitively; `representative` removes a term only when it is similar to a kept,
    better-ranked term."""

    @model_validator(mode="before")
    @classmethod
    def _removed(cls, data: Any) -> Any:
        return _reject_removed("redundancy", data)


class GeneMappingConfig(Model):
    use_hgnc_alias_mapping: bool = True
    """Resolve source symbols through HGNC alias symbols (unambiguous aliases only)."""
    preserve_unmapped_symbols_in_symbol_gmt: bool = True
    """Keep source symbols that HGNC cannot resolve in `symbols.gmt` (they are always dropped from Ensembl
    namespaces)."""
    non_gene_tokens: list[str] = Field(default=DEFAULT_NON_GENE_TOKENS)
    """Member names that are not genes (such as \"disease\") and are dropped with an audit row; `extra_non_gene_tokens`
    appends."""
    non_gene_token_regex: list[str] = Field(default=[])
    """Patterns for non-gene member names; `extra_non_gene_token_regex` appends."""
    manual_symbol_repairs: dict[str, str] = Field(default={})
    """Source symbol -> HGNC symbol corrections; they override previous-symbol and alias matches but never an exact
    current symbol."""

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
    """Use NCBI Gene info and history to rescue symbols HGNC cannot resolve."""
    allowed_source_tags: list[str] = Field(default=["CCLE", "NCI60", "DEPMAP"])
    """Source tags whose memberships may use the NCBI rescue (legacy cell-line libraries by default)."""
    gene_info_url: str = NCBI_GENE_INFO_URL
    """NCBI `Homo_sapiens.gene_info.gz` download URL."""
    gene_history_url: str = NCBI_GENE_HISTORY_URL
    """NCBI `gene_history.gz` download URL."""
    use_gene_info_synonyms: bool = True
    """Also match NCBI gene synonyms, not only official and historical symbols."""


class EnsemblArchiveConfig(Model):
    enabled: bool = True
    """Look up retired source Ensembl IDs in the Ensembl archive to find their current replacement."""
    endpoint: str = ENSEMBL_ARCHIVE_ENDPOINT
    """Ensembl REST archive endpoint."""
    batch_size: PositiveInt = 50
    """IDs per archive request."""
    connect_timeout_seconds: PositiveInt = 30
    """Connection timeout per request."""
    timeout_seconds: PositiveInt = 120
    """Read timeout per request."""
    retry_attempts: PositiveInt = 3
    """Attempts per batch before the lookup fails."""
    retry_backoff_seconds: Annotated[float, Field(ge=0.0)] = 5
    """Base delay between retries (grows exponentially)."""
    unresolved_print_limit: NonNegativeInt = 20
    """How many unresolved archive IDs to log as examples."""

    @model_validator(mode="before")
    @classmethod
    def _removed(cls, data: Any) -> Any:
        return _reject_removed("ensembl_archive", data)


class HgncConfig(Model):
    """Where the HGNC complete set comes from. Pin ``archive_date`` (YYYY-MM-DD) for a dated monthly archive."""

    url: str | None = None
    """Override the HGNC complete-set download URL."""
    archive_date: Annotated[str, Field(pattern=r"^[0-9]{4}-[0-9]{2}-[0-9]{2}$")] | None = None
    """Use the dated monthly HGNC archive (YYYY-MM-DD) instead of the current file."""


class RuntimeConfig(Model):
    workers: PositiveInt = 1
    """Worker processes for redundancy filtering; `--workers` overrides."""


class AlienConfig(Model):
    """The complete, validated configuration for one build."""

    project: ProjectConfig = ProjectConfig()
    """Names and locations."""
    outputs: OutputsConfig = OutputsConfig()
    """What the build writes."""
    term_id_collisions: TermIdCollisionsConfig = TermIdCollisionsConfig()
    """Handling of conflicting term IDs."""
    sources: list[SourceSpec]
    """Gene-set sources; at least one is required. Each entry is selected by its `type`."""
    targets: list[TargetSpec] = Field(default_factory=lambda: [TargetSpec.model_validate(t) for t in DEFAULT_TARGETS])
    """Output namespaces (one GMT each)."""
    filtering: FilteringConfig = FilteringConfig()
    """Size and broad-term filters."""
    redundancy: RedundancyConfig = RedundancyConfig()
    """Duplicate and overlap removal."""
    source_priority: dict[str, list[str]] = Field(default=DEFAULT_SOURCE_PRIORITY)
    """Per family, the order in which source tags win when terms are redundant; `extra.<family>` appends."""
    gene_mapping: GeneMappingConfig = GeneMappingConfig()
    """Symbol resolution options."""
    hgnc: HgncConfig = HgncConfig()
    """HGNC resource location."""
    ncbi_gene: NcbiGeneConfig = NcbiGeneConfig()
    """NCBI Gene rescue."""
    ensembl_archive: EnsemblArchiveConfig = EnsemblArchiveConfig()
    """Ensembl archive lookups for retired IDs."""
    runtime: RuntimeConfig = RuntimeConfig()
    """Parallelism."""

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
