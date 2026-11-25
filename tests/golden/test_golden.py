"""Golden snapshot tests: every output file of the offline fixture build is compared byte-for-byte.

Set ``ALIEN_UPDATE_GOLDEN=1`` to rewrite the expected files after an intended output change.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import shutil
from pathlib import Path

import pytest

from alien import build

FIXTURE = Path(__file__).with_name("fixture")
EXPECTED = Path(__file__).with_name("expected")
UPDATE = os.environ.get("ALIEN_UPDATE_GOLDEN") == "1"
VOLATILE_JSON_KEYS = {"created_at", "packages", "python_version", "platform", "alien_version", "git_describe"}


def _normalize_json(value: object) -> object:
    if isinstance(value, dict):
        return {key: _normalize_json(item) for key, item in value.items() if key not in VOLATILE_JSON_KEYS}
    if isinstance(value, list):
        return [_normalize_json(item) for item in value]
    return value


def _drop_column(text: str, column: str) -> str:
    lines = text.splitlines()
    if not lines or column not in lines[0].split("\t"):
        return text
    index = lines[0].split("\t").index(column)
    kept = ["\t".join(cell for i, cell in enumerate(line.split("\t")) if i != index) for line in lines]
    return "\n".join(kept) + "\n"


def _normalized_text(path: Path, root: Path) -> str:
    data = gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes()
    text = data.decode("utf-8").replace(str(root), "<ROOT>")
    name = path.name[:-3] if path.name.endswith(".gz") else path.name
    text = re.sub(r"\bworkers: [0-9]+", "workers: <N>", text)
    if name == "term_manifest.tsv":
        text = _drop_column(text, "created_at")
    if path.name == "source_provenance.json":
        text = json.dumps(_normalize_json(json.loads(text)), indent=2, sort_keys=True) + "\n"
    if path.name == "MANIFEST.json":
        manifest = json.loads(text)
        manifest["files"] = {
            name: entry
            for name, entry in manifest.get("files", {}).items()
            if not name.endswith("source_provenance.json")
        }
        text = json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    return text


def _build(tmp_path: Path, mode: str, workers: int, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    root = tmp_path / "fixture"
    shutil.copytree(FIXTURE, root)
    monkeypatch.chdir(root)
    outdir = root / "out"
    build(config=Path("config.yml"), outdir=outdir, workers=workers, output_mode=mode)
    return root, outdir


@pytest.mark.parametrize("workers", [1, 4])
@pytest.mark.parametrize("mode", ["full", "minimal", "gmt"])
def test_golden_outputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mode: str, workers: int) -> None:
    root, outdir = _build(tmp_path, mode, workers, monkeypatch)
    expected_dir = EXPECTED / mode
    produced = {path.relative_to(outdir).as_posix(): path for path in sorted(outdir.rglob("*")) if path.is_file()}
    if UPDATE and workers == 1:
        if expected_dir.exists():
            shutil.rmtree(expected_dir)
        for rel, path in produced.items():
            target = expected_dir / (rel[:-3] if rel.endswith(".gz") else rel)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(_normalized_text(path, root), encoding="utf-8")
        return
    expected = {
        path.relative_to(expected_dir).as_posix(): path for path in sorted(expected_dir.rglob("*")) if path.is_file()
    }
    produced_names = {rel[:-3] if rel.endswith(".gz") else rel for rel in produced}
    assert produced_names == set(expected)
    for rel, path in produced.items():
        name = rel[:-3] if rel.endswith(".gz") else rel
        assert _normalized_text(path, root) == expected[name].read_text(encoding="utf-8"), name
