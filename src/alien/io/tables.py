"""Deterministic TSV / TSV.GZ tables with declared column order."""

from __future__ import annotations

import gzip
import io
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from .atomic import atomic_write_bytes


def with_columns(frame: pd.DataFrame | None, columns: Sequence[str]) -> pd.DataFrame:
    """Return a copy of ``frame`` with exactly ``columns`` (missing ones filled with ``""``)."""
    if frame is None or frame.empty:
        return pd.DataFrame(columns=list(columns))
    result = frame.copy()
    for column in columns:
        if column not in result:
            result[column] = ""
    return result[list(columns)]


def gzip_bytes(data: bytes) -> bytes:
    """Gzip with a zero mtime and no embedded file name, so equal input gives equal bytes."""
    buffer = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buffer, mtime=0) as handle:
        handle.write(data)
    return buffer.getvalue()


def table_bytes(frame: pd.DataFrame) -> bytes:
    return frame.to_csv(sep="\t", index=False, lineterminator="\n").encode("utf-8")


def write_table(frame: pd.DataFrame, path: Path, columns: Sequence[str] | None = None) -> None:
    """Write ``frame`` as TSV (or deterministic TSV.GZ when ``path`` ends in ``.gz``), atomically."""
    if columns is not None:
        frame = with_columns(frame, columns)
    data = table_bytes(frame)
    atomic_write_bytes(path, gzip_bytes(data) if path.name.endswith(".gz") else data)


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
