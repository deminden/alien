"""5.7: relative paths in a config file resolve against that file's directory."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from alien import build

from conftest import MiniWorld, read_gmt


def test_relative_paths_resolve_against_config_file(mini_world: MiniWorld, monkeypatch: pytest.MonkeyPatch) -> None:
    project = mini_world.root / "proj"
    project.mkdir()
    (project / "a.gmt").write_text("T1\tdesc\tG1\tG2\tG3\n", encoding="utf-8")
    cfg = mini_world.config(project / "a.gmt")
    cfg["sources"][0]["path"] = "a.gmt"
    cfg["project"]["source_dir"] = "../sources"
    cfg["project"]["outdir"] = "out"
    cfg["targets"][0]["annotation"]["path"] = "../genes.gtf"
    (project / "c.yml").write_text(yaml.safe_dump(cfg), encoding="utf-8")
    elsewhere = mini_world.root / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)

    result = build(config=Path("../proj/c.yml"))

    assert result.outdir.resolve() == (project / "out").resolve()
    assert "LOCAL__T1" in read_gmt(project / "out" / "gmt" / "mini.gmt")
