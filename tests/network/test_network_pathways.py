"""Nightly integration test (``pytest -m network``): fetch the pathways preset, then build it frozen."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from alien import build, fetch
from alien.outputs import verify_output_dir

pytestmark = pytest.mark.network


def test_fetch_then_frozen_pathways_build(tmp_path: Path) -> None:
    config = tmp_path / "pathways.yml"
    config.write_text(
        yaml.safe_dump(
            {
                "extends": "pathways",
                "project": {"source_dir": str(tmp_path / "sources"), "outdir": str(tmp_path / "out")},
                "runtime": {"workers": 2},
            }
        ),
        encoding="utf-8",
    )
    lock = fetch(config)
    assert lock.exists()
    result = build(config, frozen=True, output_mode="minimal")
    assert result.term_counts["gtex_v11_gencode47"] > 1000
    assert verify_output_dir(result.outdir) == []
