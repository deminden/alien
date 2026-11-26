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


def test_config_doc_is_generated_from_models() -> None:
    result = subprocess.run(
        [sys.executable, "scripts/generate_config_doc.py", "--check"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


def test_every_config_option_is_described() -> None:
    from pydantic import BaseModel

    from alien.config import schema

    missing = []
    seen: set[type[BaseModel]] = set()
    stack: list[type[BaseModel]] = [schema.AlienConfig, *schema.SourceSpec.__args__[0].__args__]
    while stack:
        model = stack.pop()
        if model in seen:
            continue
        seen.add(model)
        for name, field in model.model_fields.items():
            if not field.description:
                missing.append(f"{model.__name__}.{name}")
            for arg in [field.annotation, *getattr(field.annotation, "__args__", ())]:
                inner = [arg, *getattr(arg, "__args__", ())]
                stack += [t for t in inner if isinstance(t, type) and issubclass(t, BaseModel)]
    assert missing == []


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
