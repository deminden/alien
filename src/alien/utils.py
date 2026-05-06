from __future__ import annotations

import gzip
import json
import logging
import multiprocessing as mp
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


LOGGER = logging.getLogger("alien")


def setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stderr,
    )


def process_pool_context() -> mp.context.BaseContext:
    # Avoid raw fork from IDE/test parents that may already have helper threads
    for method in ["forkserver", "spawn"]:
        try:
            return mp.get_context(method)
        except ValueError:
            continue
    return mp.get_context()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def ensure_dirs(*paths: Path) -> None:
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def strip_ensembl_version(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if not text or text.upper() in {"NA", "NAN", "NULL"}:
        return ""
    match = re.match(r"^(ENSG[0-9]+)(?:\.[0-9]+)?(?:_PAR_Y)?$", text)
    return match.group(1) if match else text


def sanitize_id(value: Any) -> str:
    text = "" if value is None else str(value)
    text = text.strip()
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    text = re.sub(r"_+", "_", text)
    return text.strip("_") or "UNNAMED"


def clean_field(value: Any) -> str:
    text = "" if value is None or pd.isna(value) else str(value)
    return re.sub(r"[\t\r\n]+", " ", text).strip()


def read_tsv(path: Path, **kwargs: Any) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, **kwargs)


def write_tsv_gz(df: pd.DataFrame, path: Path) -> None:
    ensure_dirs(path.parent)
    df.to_csv(path, sep="\t", index=False, compression="gzip")


def parse_gmt(path: Path) -> list[dict[str, object]]:
    # Read GMT records
    records: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n\r")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                raise ValueError(f"GMT line {line_number} in {path} has fewer than 3 fields.")
            records.append(
                {
                    "term_id": sanitize_id(parts[0]),
                    "description": clean_field(parts[1]),
                    "genes": [gene.strip() for gene in parts[2:] if gene.strip()],
                }
            )
    return records


def write_gmt(term_genes: dict[str, dict[str, object]], path: Path) -> None:
    # Write deterministic GMT files
    ensure_dirs(path.parent)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for term_id in sorted(term_genes):
            description = clean_field(term_genes[term_id].get("description", ""))
            genes = sorted({str(gene) for gene in term_genes[term_id].get("genes", []) if str(gene)})
            if not genes:
                continue
            fields = [_safe_existing_term_id(term_id), description] + genes
            handle.write("\t".join(fields) + "\n")


def description_from_record(record: dict[str, object]) -> str:
    fields = [
        ("source", record.get("source_tag") or record.get("source")),
        ("collection", record.get("collection")),
        ("aspect", record.get("aspect")),
        ("original_name", record.get("original_name")),
        ("db_version", record.get("db_version")),
    ]
    return ";".join(f"{key}={clean_field(value)}" for key, value in fields if clean_field(value))


def write_json(obj: Any, path: Path) -> None:
    ensure_dirs(path.parent)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2, sort_keys=True)
        handle.write("\n")


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def download_file(url: str, path: Path, force: bool = False, timeout: int = 3600) -> bool:
    # Download source files with retry support
    ensure_dirs(path.parent)
    if path.exists() and not force:
        return False
    LOGGER.info("Downloading %s", url)
    session = requests.Session()
    retry = Retry(total=3, connect=3, read=3, backoff_factor=2, status_forcelist=[429, 500, 502, 503, 504])
    session.mount("http://", HTTPAdapter(max_retries=retry))
    session.mount("https://", HTTPAdapter(max_retries=retry))
    with session.get(url, stream=True, timeout=(15, timeout)) as response:
        response.raise_for_status()
        tmp = path.with_suffix(path.suffix + ".tmp")
        with tmp.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
        shutil.move(tmp, path)
    return True


def empty_term_gene_frame() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "term_id",
            "original_name",
            "display_name",
            "description",
            "source",
            "source_tag",
            "collection",
            "subcollection",
            "family",
            "aspect",
            "db_version",
            "source_url",
            "source_license_note",
            "gene_symbol",
            "gene_ensembl_from_source",
            "metadata_json",
        ]
    )


def concat_frames(frames: list[pd.DataFrame]) -> pd.DataFrame:
    frames = [frame for frame in frames if frame is not None and not frame.empty]
    if not frames:
        return empty_term_gene_frame()
    return pd.concat(frames, ignore_index=True)


def first_existing_column(columns: list[str], candidates: list[str]) -> str | None:
    lowered = {col.lower(): col for col in columns}
    for candidate in candidates:
        if candidate in columns:
            return candidate
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]
    return None


def _safe_existing_term_id(term_id: str) -> str:
    text = clean_field(term_id)
    if not text:
        return "UNNAMED"
    return "".join(ch if ch.isalnum() or ch in "_.-" else "_" for ch in text)
