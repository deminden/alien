from __future__ import annotations

from pathlib import Path

import pytest

from alien.errors import SourceError
from alien.io.gmt import description_from_meta, parse_gmt, write_gmt
from alien.model import TermRecord


def test_gmt_roundtrip_and_description_sanitization(tmp_path: Path) -> None:
    path = tmp_path / "in.gmt"
    path.write_text("TERM 1\tdesc\twith newline\tA\tB\n\n", encoding="utf-8")
    records = parse_gmt(path)
    assert records[0].term_id == "TERM_1"
    assert records[0].genes == ("with newline", "A", "B")

    out = tmp_path / "out.gmt"
    write_gmt(
        {
            "TERM__1": TermRecord("TERM__1", "a\tb\nc", {}, {"B", "A"}),
            "EMPTY": TermRecord("EMPTY", "nothing", {}, set()),
            "": TermRecord("", "unnamed", {}, {"C"}),
        },
        out,
    )
    assert out.read_text(encoding="utf-8") == "UNNAMED\tunnamed\tC\nTERM__1\ta b c\tA\tB\n"


def test_short_gmt_line_is_an_error(tmp_path: Path) -> None:
    path = tmp_path / "bad.gmt"
    path.write_text("ONLY\tTWO\n", encoding="utf-8")
    with pytest.raises(SourceError, match="fewer than 3 fields"):
        parse_gmt(path)


def test_description_from_meta() -> None:
    meta = {"source_tag": "REACTOME", "collection": "C2", "aspect": "", "original_name": "x\ty"}
    assert description_from_meta(meta) == "source=REACTOME;collection=C2;original_name=x y"
