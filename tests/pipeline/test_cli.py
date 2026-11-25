from pathlib import Path
from types import SimpleNamespace

from alien import cli


def test_cli_uses_config_outdir_when_outdir_override_is_omitted(tmp_path, monkeypatch):
    config = tmp_path / "config.yml"
    config.write_text("project:\n  outdir: configured-out\n", encoding="utf-8")
    call = {}

    def fake_build(config=None, **kwargs):
        call.update(kwargs, config=config)
        return SimpleNamespace(outdir=Path("configured-out"), namespaces=("symbols",), plan=())

    monkeypatch.setattr(cli, "build", fake_build)

    assert cli.main(["build", "--config", str(config), "--dry-run"]) == 0
    assert call["config"] == config
    assert call["outdir"] is None


def test_cli_outdir_argument_overrides_config_outdir(tmp_path, monkeypatch):
    config = tmp_path / "config.yml"
    override = tmp_path / "override-out"
    config.write_text("project:\n  outdir: configured-out\n", encoding="utf-8")
    call = {}

    def fake_build(config=None, **kwargs):
        call.update(kwargs, config=config)
        return SimpleNamespace(outdir=override, namespaces=("symbols",), plan=())

    monkeypatch.setattr(cli, "build", fake_build)

    assert cli.main(["build", "--config", str(config), "--outdir", str(override), "--dry-run"]) == 0
    assert call["outdir"] == override


def test_cli_accepts_positional_preset_with_overrides(tmp_path, monkeypatch):
    overlay = tmp_path / "overlay.yml"
    source_dir = tmp_path / "sources"
    outdir = tmp_path / "out"
    overlay.write_text("filtering:\n  min_size_default: 5\n", encoding="utf-8")
    call = {}

    def fake_build(config=None, **kwargs):
        call.update(kwargs, config=config)
        return SimpleNamespace(outdir=outdir, namespaces=("tcga_recount3_gencode29",), plan=())

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

    def fake_build(config=None, **kwargs):
        call.update(kwargs, config=config)
        return SimpleNamespace(outdir=Path("configured-out"), namespaces=("symbols",), plan=())

    monkeypatch.setattr(cli, "build", fake_build)

    assert cli.main(["build", str(config), "--dry-run"]) == 0
    assert call["config"] == config


def test_cli_config_option_treats_value_as_path(monkeypatch):
    call = {}

    def fake_build(config=None, **kwargs):
        call.update(kwargs, config=config)
        return SimpleNamespace(outdir=Path("configured-out"), namespaces=("symbols",), plan=())

    monkeypatch.setattr(cli, "build", fake_build)

    assert cli.main(["build", "--config", "cancer", "--dry-run"]) == 0
    assert call["config"] == Path("cancer")


def test_cli_presets_lists_bundled_presets(capsys):
    assert cli.main(["presets"]) == 0
    output = capsys.readouterr().out
    assert "cancer_dependency" in output
    assert "cancer" in output
    assert "function_location" in output


def test_cli_maps_errors_to_exit_codes(tmp_path, capsys):
    bad = tmp_path / "bad.yml"
    bad.write_text("filterring: {}\nsources: [{type: symbol_gmt, path: x.gmt}]\n", encoding="utf-8")
    assert cli.main(["validate", str(bad)]) == 2
    assert "filterring: unknown option" in capsys.readouterr().err
    missing = tmp_path / "missing.yml"
    missing.write_text(
        f"sources: [{{type: symbol_gmt, path: {tmp_path / 'nope.gmt'}}}]\ntargets: []\n", encoding="utf-8"
    )
    assert cli.main(["build", str(missing), "--outdir", str(tmp_path / "out")]) == 3
    assert "nope.gmt" in capsys.readouterr().err


def test_cli_validate_version_and_module_entry(tmp_path, capsys):
    import subprocess
    import sys

    good = tmp_path / "good.yml"
    good.write_text("sources: [{type: symbol_gmt, path: x.gmt}]\ntargets: []\n", encoding="utf-8")
    assert cli.main(["validate", str(good)]) == 0
    assert "Configuration is valid" in capsys.readouterr().out
    result = subprocess.run([sys.executable, "-m", "alien", "--version"], capture_output=True, text=True, check=False)
    assert result.returncode == 0 and result.stdout.startswith("alien ")


def test_cli_requires_config(capsys):
    import pytest

    with pytest.raises(SystemExit):
        cli.main(["build"])
    assert "requires a preset name" in capsys.readouterr().err


def test_cli_unexpected_error_is_one_line(monkeypatch, capsys):
    def boom(config=None, **kwargs):
        raise KeyError("x")

    monkeypatch.setattr(cli, "build", boom)
    assert cli.main(["build", "pathways"]) == 1
    err = capsys.readouterr().err
    assert "unexpected error" in err and "Traceback" not in err
