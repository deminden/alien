"""Target namespace specifications: annotation files and output-gene lists."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from pydantic import Field, field_validator, model_validator

from .base import Model, _blank_to_none, _deprecated, _rename


class GtfAttributes(Model):
    gene_id: list[str] = Field(default=["gene_id"])
    gene_name: list[str] = Field(default=["gene_name"])
    gene_biotype: list[str] = Field(default=["gene_type", "gene_biotype"])

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
    metadata_fallbacks: list[MetadataFallbackSpec] = Field(default=[])

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
        if not (self.path or self.url) and (self.source or "").upper() != "GENCODE":
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
