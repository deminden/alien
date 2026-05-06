from pathlib import Path


def test_no_old_project_identity_references():
    root = Path(__file__).resolve().parents[1]
    ignored_dirs = {".git", ".pytest_cache", "dist", "build"}
    hits = []
    for path in root.rglob("*"):
        if not path.is_file() or any(part in ignored_dirs for part in path.parts):
            continue
        if path.suffix in {".gz", ".parquet", ".pyc"}:
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        legacy = "corr" + "sea"
        if legacy in text.casefold():
            hits.append(str(path.relative_to(root)))
    assert hits == []
