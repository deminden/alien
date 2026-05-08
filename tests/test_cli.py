from pathlib import Path
from types import SimpleNamespace

from alien import cli


def test_cli_uses_config_outdir_when_outdir_override_is_omitted(tmp_path, monkeypatch):
    config = tmp_path / "config.yml"
    config.write_text("project:\n  outdir: configured-out\n", encoding="utf-8")
    call = {}

    def fake_build(**kwargs):
        call.update(kwargs)
        return SimpleNamespace(outdir=Path("configured-out"), namespaces=("symbols",))

    monkeypatch.setattr(cli, "build", fake_build)

    assert cli.main(["build", "--config", str(config), "--dry-run"]) == 0
    assert call["config"] == config
    assert call["outdir"] is None


def test_cli_outdir_argument_overrides_config_outdir(tmp_path, monkeypatch):
    config = tmp_path / "config.yml"
    override = tmp_path / "override-out"
    config.write_text("project:\n  outdir: configured-out\n", encoding="utf-8")
    call = {}

    def fake_build(**kwargs):
        call.update(kwargs)
        return SimpleNamespace(outdir=override, namespaces=("symbols",))

    monkeypatch.setattr(cli, "build", fake_build)

    assert cli.main(["build", "--config", str(config), "--outdir", str(override), "--dry-run"]) == 0
    assert call["outdir"] == override
