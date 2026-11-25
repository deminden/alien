"""Atomic file and directory writes.

Every writer goes through a unique temporary file in the destination directory followed by
``os.replace``, so readers never observe a half-written file and concurrent writers never share a
temporary name.
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any


def atomic_write_bytes(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
        tmp.replace(path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def atomic_write_text(path: Path, text: str) -> None:
    atomic_write_bytes(path, text.encode("utf-8"))


def dumps_json(obj: Any) -> str:
    """Deterministic JSON: sorted keys, two-space indent, trailing newline."""
    return json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=False, default=str) + "\n"


def atomic_write_json(path: Path, obj: Any) -> None:
    atomic_write_text(path, dumps_json(obj))


def unique_sibling(path: Path, label: str) -> Path:
    """A fresh, not-yet-existing directory name next to ``path`` on the same filesystem."""
    path.parent.mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(prefix=f".{path.name}.{label}-", dir=path.parent))


def atomic_dir_swap(staging: Path, destination: Path) -> None:
    """Replace ``destination`` with ``staging``; the previous tree is removed only after the swap."""
    if not destination.exists():
        staging.replace(destination)
        return
    previous = unique_sibling(destination, "previous")
    previous.rmdir()
    destination.replace(previous)
    try:
        staging.replace(destination)
    except BaseException:
        previous.replace(destination)
        raise
    shutil.rmtree(previous, ignore_errors=True)
