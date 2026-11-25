"""5.11: builds are byte-for-byte reproducible."""

from __future__ import annotations

from pathlib import Path

import pytest

from alien import build

from conftest import MiniWorld


def _files(root: Path) -> dict[str, bytes]:
    return {path.relative_to(root).as_posix(): path.read_bytes() for path in root.rglob("*") if path.is_file()}


def test_two_builds_differ_only_in_provenance(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"T1": ["G1", "G2", "G3"], "T2": ["G4", "G5"]})
    one = _files(build(config=mini_world.config(gmt), outdir=mini_world.root / "one").outdir)
    two = _files(build(config=mini_world.config(gmt), outdir=mini_world.root / "two").outdir)
    assert set(one) == set(two)
    differing = {name for name in one if one[name] != two[name]}
    assert differing <= {"metadata/source_provenance.json", "MANIFEST.json"}


def test_source_date_epoch_makes_every_byte_identical(mini_world: MiniWorld, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("SOURCE_DATE_EPOCH", "1700000000")
    gmt = mini_world.gmt("lib.gmt", {"T1": ["G1", "G2", "G3"]})
    one = _files(build(config=mini_world.config(gmt), outdir=mini_world.root / "one").outdir)
    two = _files(build(config=mini_world.config(gmt), outdir=mini_world.root / "two").outdir)
    assert one == two
    assert b"2023-11-14T22:13:20" in one["metadata/source_provenance.json"]
