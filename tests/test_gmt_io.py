from alien.utils import parse_gmt, write_gmt


def test_gmt_roundtrip_and_description_sanitization(tmp_path):
    path = tmp_path / "in.gmt"
    path.write_text("TERM 1\tdesc\twith newline\tA\tB\n", encoding="utf-8")
    records = parse_gmt(path)
    assert records[0]["term_id"] == "TERM_1"
    assert records[0]["genes"] == ["with newline", "A", "B"]

    out = tmp_path / "out.gmt"
    write_gmt({"TERM__1": {"description": "a\tb\nc", "genes": {"B", "A"}}}, out)
    assert out.read_text(encoding="utf-8") == "TERM__1\ta b c\tA\tB\n"
