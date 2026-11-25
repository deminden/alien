from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_outputs_doc_is_generated_from_schema() -> None:
    result = subprocess.run(
        [sys.executable, "scripts/generate_outputs_doc.py", "--check"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


def test_every_written_file_is_documented(tmp_path: Path, monkeypatch: object) -> None:
    import shutil

    from alien import build
    from alien.schema import OUTPUT_FILES

    fixture = ROOT / "tests" / "golden" / "fixture"
    shutil.copytree(fixture, tmp_path / "fixture")
    result = build(tmp_path / "fixture" / "config.yml", outdir=tmp_path / "out")
    documented = [f.path for f in OUTPUT_FILES]
    for path in result.outdir.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(result.outdir).as_posix()
        generic = rel.replace(f"gene_mapping_{rel.split('gene_mapping_')[-1].split('.')[0]}", "gene_mapping_<target>")
        generic = generic if not rel.startswith("gmt/") else "gmt/<namespace>.gmt"
        assert generic in documented, rel
