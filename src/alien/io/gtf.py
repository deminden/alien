"""GTF gene-record parsing."""

from __future__ import annotations

import gzip
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import IO, Final

import pandas as pd

from ..model import ensembl_stable_id
from ..schema import ANNOTATION_COLUMNS

GTF_FIELDS: Final = 9
KEY_VALUE_PARTS: Final = 2
DEFAULT_GTF_ATTRIBUTES: Final[Mapping[str, Sequence[str]]] = {
    "gene_id": ("gene_id",),
    "gene_name": ("gene_name",),
    "gene_biotype": ("gene_type", "gene_biotype"),
}


def open_text(path: Path) -> IO[str]:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_attributes(text: str) -> dict[str, str]:
    """GTF (``key "value";``) or GFF3 (``key=value;``) attributes."""
    attrs: dict[str, str] = {}
    for raw_item in text.split(";"):
        item = raw_item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        else:
            parts = item.split(None, 1)
            if len(parts) != KEY_VALUE_PARTS:
                continue
            key, value = parts
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def _first(attrs: Mapping[str, str], names: Sequence[str]) -> str:
    for name in names:
        value = attrs.get(name, "")
        if value:
            return value
    return ""


def parse_gtf_genes(path: Path, version: str, attributes: Mapping[str, Sequence[str]] | None = None) -> pd.DataFrame:
    """One row per ``gene`` feature, with stable and versioned Ensembl IDs."""
    names = attributes or DEFAULT_GTF_ATTRIBUTES
    rows: list[dict[str, object]] = []
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < GTF_FIELDS or parts[2] != "gene":
                continue
            attrs = parse_attributes(parts[8])
            versioned = _first(attrs, names["gene_id"])
            rows.append(
                {
                    "gencode_version": f"v{version}",
                    "ensembl_gene_id_versioned": versioned,
                    "ensembl_gene_id": ensembl_stable_id(versioned),
                    "gene_symbol": _first(attrs, names["gene_name"]),
                    "gene_biotype": _first(attrs, names["gene_biotype"]),
                    "seqname": parts[0],
                    "start": parts[3],
                    "end": parts[4],
                    "strand": parts[6],
                    "source_gtf": str(path),
                    "target_id_metadata_source": "annotation_gtf",
                    "is_in_output_genes": False,
                }
            )
    return pd.DataFrame(rows, columns=list(ANNOTATION_COLUMNS)).drop_duplicates()
