"""5.10: --dry-run validates the configuration and inputs without building."""

from __future__ import annotations

import pytest

from alien import build

from conftest import MiniWorld


def test_dry_run_fails_on_missing_local_path(mini_world: MiniWorld) -> None:
    cfg = mini_world.config(mini_world.root / "does_not_exist.gmt")
    with pytest.raises(Exception, match=r"does_not_exist\.gmt"):
        build(config=cfg, dry_run=True)


def test_dry_run_reports_plan_without_writing(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"T1": ["G1", "G2"]})
    result = build(config=mini_world.config(gmt), dry_run=True)
    assert not result.outdir.exists()
    plan = "\n".join(result.plan)
    assert "hgnc_complete_set.txt" in plan
    assert "cached" in plan
    assert "gmt/mini.gmt" in plan
