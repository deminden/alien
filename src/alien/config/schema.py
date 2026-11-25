"""Typed ALIEN configuration. Defaults live here and nowhere else.

Every model forbids unknown keys, so a typo such as ``jacard_cutoff`` fails with the YAML path of the
offending key instead of being silently ignored.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path
from typing import Annotated, Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

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

#: Options that existed in 0.1.x but never had an effect; using one is an error naming the replacement.
REMOVED_OPTIONS: dict[str, dict[str, str]] = {
    "redundancy": {
        "apply_after_output_gene_projection": "redundancy always runs after target projection",
        "nested_overlap_filter": "the nested-overlap filter was never implemented; remove the block",
    },
    "filtering": {
        "remove_empty_after_output_gene_projection": "empty terms are always dropped; use min_size_default",
    },
    "gene_mapping": {
        "strip_ensembl_version": "Ensembl versions are always stripped",
        "allow_ambiguous_alias_mapping": "ambiguous aliases are always dropped and audited",
        "prefer_gencode_exact_gene_name": "annotation gene names are always preferred",
        "prefer_current_hgnc_symbol": "current HGNC symbols are always preferred",
        "drop_unmapped_from_ensg_gmt": "unmapped genes are never written to Ensembl GMTs",
    },
    "ensembl_archive": {"unresolved_example_threshold": "use ensembl_archive.unresolved_print_limit"},
}

PositiveInt = Annotated[int, Field(ge=1)]
NonNegativeInt = Annotated[int, Field(ge=0)]


def _deprecated(old: str, new: str) -> None:
    warnings.warn(f"Config key {old!r} is deprecated; use {new!r}.", DeprecationWarning, stacklevel=4)


def _rename(data: Any, aliases: dict[str, str], where: str) -> Any:
    """Map deprecated key spellings to their canonical names, warning once per key."""
    if not isinstance(data, dict):
        return data
    data = dict(data)
    for old, new in aliases.items():
        if old in data:
            _deprecated(f"{where}.{old}", f"{where}.{new}")
            value = data.pop(old)
            data.setdefault(new, value)
    return data


def _reject_removed(section: str, data: Any) -> Any:
    if isinstance(data, dict):
        for key, replacement in REMOVED_OPTIONS.get(section, {}).items():
            if key in data:
                raise ValueError(f"option {section}.{key} was removed: {replacement}")
    return data


def _blank_to_none(value: Any) -> Any:
    if isinstance(value, str) and value.strip() in {"", "null"}:
        return None
    return value


class Model(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True, populate_by_name=True)


class ProjectConfig(Model):
    name: str = "alien"
    source_dir: Path = Path("data/alien_sources")
    outdir: Path = Path("data/alien_gmt")


class OutputsConfig(Model):
    include_symbols: bool = True
    mode: Literal["full", "minimal", "gmt"] = "full"


class TermIdCollisionsConfig(Model):
    action: Literal["error", "warn", "merge"] = "error"


# --------------------------------------------------------------------------------------------- sources


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


SourceSpec = Annotated[
    SymbolGmtSource
    | MsigdbTsvSource
    | CanonicalTsvSource
    | MsigdbCacheSource
    | MsigdbRemoteSource
    | EnrichrRemoteSource,
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
)

# --------------------------------------------------------------------------------------------- targets


class GtfAttributes(Model):
    gene_id: list[str] = ["gene_id"]
    gene_name: list[str] = ["gene_name"]
    gene_biotype: list[str] = ["gene_type", "gene_biotype"]

    @field_validator("gene_id", "gene_name", "gene_biotype", mode="before")
    @classmethod
    def _listify(cls, value: Any) -> Any:
        return [value] if isinstance(value, str) else value


class AnnotationFile(Model):
    source: str | None = None
    version: str | None = None
    path: Path | None = None
    url: str | None = None
    file_name: str | None = None
    attributes: GtfAttributes = GtfAttributes()

    @field_validator("version", mode="before")
    @classmethod
    def _version_text(cls, value: Any) -> Any:
        return str(value) if isinstance(value, int | float) else value

    @field_validator("path", mode="before")
    @classmethod
    def _blank_path(cls, value: Any) -> Any:
        return _blank_to_none(value)


class MetadataFallbackSpec(AnnotationFile):
    mode: Literal["fill_missing_output_metadata"] = "fill_missing_output_metadata"


class AnnotationSpec(AnnotationFile):
    metadata_fallbacks: list[MetadataFallbackSpec] = []

    @model_validator(mode="before")
    @classmethod
    def _checks(cls, data: Any) -> Any:
        if isinstance(data, dict):
            if "supplements" in data:
                raise ValueError("annotation.supplements is not supported; use annotation.metadata_fallbacks")
            fallbacks = data.get("metadata_fallbacks")
            if isinstance(fallbacks, dict):
                data = {**data, "metadata_fallbacks": [fallbacks]}
            elif fallbacks is None and "metadata_fallbacks" in data:
                data = {**data, "metadata_fallbacks": []}
        return data

    @model_validator(mode="after")
    def _has_origin(self) -> AnnotationSpec:
        if not (self.path or self.url or self.version):
            raise ValueError("annotation must define path, url, or version")
        if not (self.path or self.url) and (self.source or "GENCODE").upper() != "GENCODE":
            raise ValueError(
                f"annotation source {self.source!r} needs annotation.path or annotation.url; "
                "only GENCODE annotations can be downloaded from a version alone"
            )
        return self


class GeneIdSpec(Model):
    """A list of Ensembl gene IDs, inline (ids) or read from a file column (path + id_column)."""

    path: Path | None = None
    id_column: str | None = None
    symbol_column: str | None = None
    ids: list[str] | None = None

    @model_validator(mode="before")
    @classmethod
    def _coerce(cls, data: Any) -> Any:
        if isinstance(data, str | Path):
            return {"path": data}
        if isinstance(data, list | tuple | set):
            return {"ids": [str(item) for item in data]}
        data = _rename(
            data,
            {
                "file": "path",
                "column": "id_column",
                "gene_column": "id_column",
                "field": "id_column",
                "gene_symbol_column": "symbol_column",
                "genes": "ids",
            },
            "output_genes/gene_filter",
        )
        if isinstance(data, dict) and isinstance(data.get("ids"), str):
            data["ids"] = [data["ids"]]
        return data

    @model_validator(mode="after")
    def _path_or_ids(self) -> GeneIdSpec:
        if self.path is None and self.ids is None:
            raise ValueError("define either path or ids")
        return self


class TargetSpec(Model):
    name: str = Field(min_length=1)
    type: Literal["ensembl_gtf"] = "ensembl_gtf"
    annotation: AnnotationSpec
    output_genes: GeneIdSpec | None = None
    gene_filter: GeneIdSpec | None = None

    @model_validator(mode="before")
    @classmethod
    def _legacy(cls, data: Any) -> Any:
        if not isinstance(data, dict):
            return data
        data = dict(data)
        replacements = {
            "restrict_to": "output_genes or gene_filter",
            "gene_universe": "output_genes",
            "output_ids": "output_genes",
            "annotation_supplements": "annotation.metadata_fallbacks",
        }
        for key, replacement in replacements.items():
            if key in data:
                raise ValueError(f"target key {key!r} is not supported; use {replacement}")
        annotation = dict(data.get("annotation") or {})
        if "annotation_gtf" in data:
            _deprecated("annotation_gtf", "annotation.path")
            annotation.setdefault("path", data.pop("annotation_gtf"))
        if "gencode_version" in data:
            _deprecated("gencode_version", "annotation.version")
            annotation.setdefault("version", data.pop("gencode_version"))
            annotation.setdefault("source", "GENCODE")
        if annotation or "annotation" in data:
            data["annotation"] = annotation
        for key in ("output_genes", "gene_filter"):
            data[key] = _blank_to_none(data.get(key))
        return data

    @model_validator(mode="after")
    def _one_gene_list(self) -> TargetSpec:
        if self.output_genes is not None and self.gene_filter is not None:
            raise ValueError(f"target {self.name} must define only one of output_genes or gene_filter")
        return self


DEFAULT_TARGETS = [{"name": "human_gencode49", "annotation": {"source": "GENCODE", "version": "49"}}]

# --------------------------------------------------------------------------------------------- sections


class FilteringConfig(Model):
    min_size_default: NonNegativeInt = 10
    max_size_default: PositiveInt = 500
    max_size_disease: PositiveInt = 500
    max_size_cancer: PositiveInt = 500
    drop_broad_disease_terms: bool = True
    broad_disease_regex: list[str] = DEFAULT_BROAD_DISEASE_REGEX

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
    non_gene_tokens: list[str] = DEFAULT_NON_GENE_TOKENS
    non_gene_token_regex: list[str] = []
    manual_symbol_repairs: dict[str, str] = {}

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
    allowed_source_tags: list[str] = ["CCLE", "NCI60", "DEPMAP"]
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
    source_priority: dict[str, list[str]] = DEFAULT_SOURCE_PRIORITY
    gene_mapping: GeneMappingConfig = GeneMappingConfig()
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
        "ncbi_gene": NcbiGeneConfig(),
        "ensembl_archive": EnsemblArchiveConfig(),
        "runtime": RuntimeConfig(),
    }
    layer: dict[str, Any] = {name: model.model_dump(mode="json") for name, model in sections.items()}
    layer["source_priority"] = {family: list(tags) for family, tags in DEFAULT_SOURCE_PRIORITY.items()}
    layer["sources"] = []
    layer["targets"] = [dict(target) for target in DEFAULT_TARGETS]
    return layer
