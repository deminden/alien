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


def test_cli_accepts_positional_preset_with_overrides(tmp_path, monkeypatch):
    overlay = tmp_path / "overlay.yml"
    source_dir = tmp_path / "sources"
    outdir = tmp_path / "out"
    overlay.write_text("filtering:\n  min_size_default: 5\n", encoding="utf-8")
    call = {}

    def fake_build(**kwargs):
        call.update(kwargs)
        return SimpleNamespace(outdir=outdir, namespaces=("tcga_recount3_gencode29",))

    monkeypatch.setattr(cli, "build", fake_build)

    assert (
        cli.main(
            [
                "build",
                "cancer",
                "--override",
                str(overlay),
                "--source-dir",
                str(source_dir),
                "--outdir",
                str(outdir),
                "--output-mode",
                "minimal",
                "--dry-run",
            ]
        )
        == 0
    )
    assert call["config"] == "cancer"
    assert call["overrides"] == [overlay]
    assert call["source_dir"] == source_dir
    assert call["outdir"] == outdir
    assert call["output_mode"] == "minimal"


def test_cli_accepts_positional_config_path(tmp_path, monkeypatch):
    config = tmp_path / "pathways.yml"
    config.write_text("project:\n  outdir: configured-out\n", encoding="utf-8")
    call = {}

    def fake_build(**kwargs):
        call.update(kwargs)
        return SimpleNamespace(outdir=Path("configured-out"), namespaces=("symbols",))

    monkeypatch.setattr(cli, "build", fake_build)

    assert cli.main(["build", str(config), "--dry-run"]) == 0
    assert call["config"] == config


def test_cli_config_option_treats_value_as_path(monkeypatch):
    call = {}

    def fake_build(**kwargs):
        call.update(kwargs)
        return SimpleNamespace(outdir=Path("configured-out"), namespaces=("symbols",))

    monkeypatch.setattr(cli, "build", fake_build)

    assert cli.main(["build", "--config", "cancer", "--dry-run"]) == 0
    assert call["config"] == Path("cancer")


def test_cli_presets_lists_bundled_presets(capsys):
    assert cli.main(["presets"]) == 0
    output = capsys.readouterr().out
    assert "cancer_dependency" in output
    assert "cancer" in output
    assert "function_location" in output
