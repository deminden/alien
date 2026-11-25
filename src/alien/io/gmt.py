"""GMT reading and deterministic writing."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from ..errors import SourceError
from ..model import TermRecord, clean_field, sanitize_term_id
from .atomic import atomic_write_text

MIN_GMT_FIELDS: Final = 3


@dataclass(frozen=True, slots=True)
class GmtRecord:
    term_id: str
    description: str
    genes: tuple[str, ...]


def parse_gmt(path: Path) -> list[GmtRecord]:
    records: list[GmtRecord] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n\r")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < MIN_GMT_FIELDS:
                raise SourceError(f"GMT line {line_number} in {path} has fewer than {MIN_GMT_FIELDS} fields.")
            genes = tuple(gene.strip() for gene in parts[2:] if gene.strip())
            records.append(GmtRecord(sanitize_term_id(parts[0]), clean_field(parts[1]), genes))
    return records


def _gmt_term_id(term_id: str) -> str:
    text = clean_field(term_id)
    if not text:
        return "UNNAMED"
    return "".join(ch if ch.isalnum() or ch in "_.-" else "_" for ch in text)


def gmt_text(terms: Mapping[str, TermRecord]) -> str:
    """Terms sorted by ID, genes sorted and de-duplicated, empty terms omitted."""
    lines: list[str] = []
    for term_id in sorted(terms):
        genes = sorted({str(gene) for gene in terms[term_id].genes if str(gene)})
        if genes:
            lines.append("\t".join([_gmt_term_id(term_id), clean_field(terms[term_id].description), *genes]))
    return "".join(f"{line}\n" for line in lines)


def write_gmt(terms: Mapping[str, TermRecord], path: Path) -> None:
    atomic_write_text(path, gmt_text(terms))


def description_from_meta(meta: Mapping[str, object]) -> str:
    fields: Iterable[tuple[str, object]] = (
        ("source", meta.get("source_tag") or meta.get("source")),
        ("collection", meta.get("collection")),
        ("aspect", meta.get("aspect")),
        ("original_name", meta.get("original_name")),
        ("db_version", meta.get("db_version")),
    )
    return ";".join(f"{key}={clean_field(value)}" for key, value in fields if clean_field(value))
