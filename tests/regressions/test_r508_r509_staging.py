"""5.8 and 5.9: failed builds keep their diagnostics; successful builds replace outdir atomically."""

from __future__ import annotations

import json

import pytest

from alien import build

from conftest import MiniWorld


def _colliding_config(mini_world: MiniWorld) -> dict[str, object]:
    gmt_a = mini_world.gmt("a.gmt", {"T1": ["G1", "G2"]})
    gmt_b = mini_world.gmt("b.gmt", {"T1": ["G3", "G4"]})
    cfg = mini_world.config(gmt_a)
    cfg["sources"] = [
        {"type": "symbol_gmt", "path": str(gmt_a), "source_tag": "LOCAL", "term_prefix": "X", "collection": "a"},
        {"type": "symbol_gmt", "path": str(gmt_b), "source_tag": "LOCAL", "term_prefix": "X", "collection": "b"},
    ]
    cfg["term_id_collisions"] = {"action": "error"}
    return cfg


@pytest.mark.parametrize("mode", ["full", "minimal", "gmt"])
def test_failed_build_keeps_diagnostics(mini_world: MiniWorld, mode: str) -> None:
    outdir = mini_world.root / "out"
    with pytest.raises(Exception, match=r"term_id_collisions\.tsv") as info:
        build(config=_colliding_config(mini_world), outdir=outdir, output_mode=mode)
    audit = outdir / "failed" / "metadata" / "term_id_collisions.tsv"
    assert str(audit) in str(info.value)
    assert audit.exists()


def test_rebuild_removes_stale_outputs_and_writes_manifest(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"T1": ["G1", "G2", "G3"]})
    outdir = mini_world.root / "out"
    build(config=mini_world.config(gmt, targets=[mini_world.target("a"), mini_world.target("b")]), outdir=outdir)
    assert (outdir / "gmt" / "b.gmt").exists()
    build(config=mini_world.config(gmt, targets=[mini_world.target("a")]), outdir=outdir)
    assert not (outdir / "gmt" / "b.gmt").exists()
    manifest = json.loads((outdir / "MANIFEST.json").read_text(encoding="utf-8"))
    assert "gmt/a.gmt" in manifest["files"]
    assert len(manifest["files"]["gmt/a.gmt"]["sha256"]) == 64
    assert not [path for path in mini_world.root.iterdir() if "staging" in path.name]


def test_refuses_to_overwrite_foreign_directory(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"T1": ["G1", "G2", "G3"]})
    outdir = mini_world.root / "precious"
    outdir.mkdir()
    (outdir / "notes.txt").write_text("keep me", encoding="utf-8")
    with pytest.raises(Exception, match="--force"):
        build(config=mini_world.config(gmt), outdir=outdir)
    assert (outdir / "notes.txt").read_text(encoding="utf-8") == "keep me"
    build(config=mini_world.config(gmt), outdir=outdir, force=True)
    assert not (outdir / "notes.txt").exists()
