from __future__ import annotations

import hashlib
from pathlib import Path

import pytest
import requests

from alien.errors import ResourceError
from alien.io import http
from alien.io.atomic import atomic_dir_swap, atomic_write_text, unique_sibling
from alien.io.tables import gzip_bytes, read_table, with_columns, write_table

from conftest import FakeServer


def test_fetch_records_hashes_and_replaces_atomically(tmp_path: Path, fake_http: FakeServer) -> None:
    fake_http.routes["https://example.test/a"] = b"hello"
    record = http.fetch("https://example.test/a", tmp_path / "a", sha256=hashlib.sha256(b"hello").hexdigest())
    assert (tmp_path / "a").read_bytes() == b"hello"
    assert record.size == 5 and record.md5 == hashlib.md5(b"hello").hexdigest() and record.etag
    with pytest.raises(ResourceError, match="sha256 verification"):
        http.fetch("https://example.test/a", tmp_path / "b", sha256="0" * 64)
    with pytest.raises(ResourceError, match="Could not download"):
        http.fetch("https://example.test/missing", tmp_path / "c")
    assert sorted(p.name for p in tmp_path.iterdir()) == ["a"]


def test_json_helpers(fake_http: FakeServer) -> None:
    fake_http.routes["https://example.test/j"] = {"a": 1}
    fake_http.routes["https://example.test/p"] = lambda payload: {"echo": payload}
    assert http.get_json("https://example.test/j") == {"a": 1}
    assert http.post_json("https://example.test/p", {"x": 2}) == {"echo": {"x": 2}}
    fake_http.routes["https://example.test/e"] = requests.ConnectionError("down")
    with pytest.raises(ResourceError, match="down"):
        http.get_json("https://example.test/e")


def test_session_retries_and_user_agent() -> None:
    session = http.build_session(retries=4, backoff=1.5)
    adapter = session.get_adapter("https://x")
    assert adapter.max_retries.total == 4  # type: ignore[attr-defined]
    assert adapter.max_retries.respect_retry_after_header  # type: ignore[attr-defined]
    assert session.headers["User-Agent"].startswith("bioalien/")


def test_atomic_writes_and_swap(tmp_path: Path) -> None:
    target = tmp_path / "dir" / "f.txt"
    atomic_write_text(target, "x")
    assert target.read_text(encoding="utf-8") == "x"
    old = tmp_path / "out"
    old.mkdir()
    (old / "old.txt").write_text("old", encoding="utf-8")
    staging = unique_sibling(old, "staging")
    (staging / "new.txt").write_text("new", encoding="utf-8")
    atomic_dir_swap(staging, old)
    assert [p.name for p in old.iterdir()] == ["new.txt"]
    assert sorted(p.name for p in tmp_path.iterdir()) == ["dir", "out"]


def test_tables_are_deterministic(tmp_path: Path) -> None:
    import pandas as pd

    frame = pd.DataFrame({"b": [1, 2], "a": ["x", "y"]})
    write_table(frame, tmp_path / "t.tsv.gz", ["a", "b", "c"])
    first = (tmp_path / "t.tsv.gz").read_bytes()
    write_table(frame, tmp_path / "t.tsv.gz", ["a", "b", "c"])
    assert (tmp_path / "t.tsv.gz").read_bytes() == first == gzip_bytes(b"a\tb\tc\nx\t1\t\ny\t2\t\n")
    assert list(read_table(tmp_path / "t.tsv.gz").columns) == ["a", "b", "c"]
    assert list(with_columns(None, ["z"]).columns) == ["z"]
    assert "c" not in frame  # with_columns never mutates its input
