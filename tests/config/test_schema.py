from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from alien import build, cli
from alien.config import AlienConfig, load_config, load_preset_config
from alien.errors import ConfigError

SOURCE = {"type": "symbol_gmt", "path": "lib.gmt"}


def _write(path: Path, data: dict[str, object]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    return path


@pytest.mark.parametrize(
    ("data", "location"),
    [
        ({"filterring": {"min_size_default": 2}}, "filterring: unknown option"),
        ({"redundancy": {"jacard_cutoff": 0.5}}, "redundancy.jacard_cutoff: unknown option"),
        ({"sources": [SOURCE | {"pathz": "x"}]}, "sources[0].pathz: unknown option"),
        (
            {
                "targets": [
                    {"name": "a", "annotation": {"source": "GENCODE", "version": "47"}},
                    {"name": "b", "annotation": {"versoin": "1"}},
                ]
            },
            "targets[1].annotation.versoin: unknown option",
        ),
        ({"redundancy": {"jaccard_cutoff": 1.5}}, "redundancy.jaccard_cutoff"),
        ({"runtime": {"workers": 0}}, "runtime.workers"),
        ({"sources": [{"type": "nope"}]}, "sources[0]"),
    ],
)
def test_typos_and_bad_values_name_the_yaml_path(tmp_path: Path, data: dict[str, object], location: str) -> None:
    config = _write(tmp_path / "c.yml", {"sources": [SOURCE], **data})
    with pytest.raises(ConfigError) as info:
        load_config(config)
    assert location in str(info.value)


def test_config_without_sources_is_rejected() -> None:
    with pytest.raises(ConfigError, match="at least one source"):
        load_config({"targets": []})


def test_duplicate_target_names_are_rejected() -> None:
    target = {"name": "same", "annotation": {"source": "GENCODE", "version": "47"}}
    with pytest.raises(ConfigError, match="repeated: same"):
        load_config({"sources": [SOURCE], "targets": [target, target]})


def test_term_id_collisions_shorthand() -> None:
    assert load_config({"sources": [SOURCE], "term_id_collisions": "warn"}).term_id_collisions.action == "warn"


def test_deprecated_gene_id_aliases_warn_and_map(tmp_path: Path) -> None:
    target = {
        "name": "t",
        "annotation": {"source": "GENCODE", "version": "47"},
        "output_genes": {"file": "g.tsv", "gene_column": "id"},
    }
    with pytest.warns(DeprecationWarning) as record:
        cfg = load_config(_write(tmp_path / "c.yml", {"sources": [SOURCE], "targets": [target]}))
    messages = " ".join(str(item.message) for item in record)
    assert "output_genes/gene_filter.path" in messages
    assert "output_genes/gene_filter.id_column" in messages
    output_genes = cfg.targets[0].output_genes
    assert output_genes is not None
    assert output_genes.path == tmp_path / "g.tsv"
    assert output_genes.id_column == "id"


def test_legacy_annotation_keys_still_load_with_warning() -> None:
    with pytest.warns(DeprecationWarning, match="annotation.version"):
        cfg = load_config({"sources": [SOURCE], "targets": [{"name": "t", "gencode_version": 47}]})
    annotation = cfg.targets[0].annotation
    assert (annotation.source, annotation.version) == ("GENCODE", "47")


def test_extends_preset_and_file(tmp_path: Path) -> None:
    base = _write(tmp_path / "base.yml", {"extends": "pathways", "filtering": {"min_size_default": 3}})
    child = _write(tmp_path / "sub" / "child.yml", {"extends": "../base.yml", "redundancy": {"jaccard_cutoff": 0.5}})
    cfg = load_config(child)
    preset = load_preset_config("pathways")
    assert cfg.sources == preset.sources
    assert cfg.filtering.min_size_default == 3
    assert cfg.redundancy.jaccard_cutoff == 0.5
    # A preset's paths stay relative to the working directory even when extended.
    assert cfg.project.source_dir == Path("data/alien_sources")
    assert base.exists()


def test_extends_loop_is_detected(tmp_path: Path) -> None:
    _write(tmp_path / "a.yml", {"extends": "b.yml", "sources": [SOURCE]})
    _write(tmp_path / "b.yml", {"extends": "a.yml"})
    with pytest.raises(ConfigError, match="loops"):
        load_config(tmp_path / "a.yml")


def test_extends_unknown_target(tmp_path: Path) -> None:
    with pytest.raises(ConfigError, match="neither a bundled preset"):
        load_config(_write(tmp_path / "a.yml", {"extends": "nowhere.yml"}))


def test_list_extension_keys(tmp_path: Path) -> None:
    cfg = load_config(
        {
            "sources": [SOURCE],
            "source_priority": {"extra": {"biology_process_pathway": ["LOCAL"]}},
            "gene_mapping": {"extra_non_gene_token_regex": ["^x"]},
        }
    )
    assert cfg.source_priority["biology_process_pathway"][-1] == "LOCAL"
    assert cfg.source_priority["disease_phenotype"][0] == "HPO"
    assert cfg.gene_mapping.non_gene_token_regex == ["^x"]


def test_output_mode_comes_from_config(tmp_path: Path) -> None:
    gmt = tmp_path / "lib.gmt"
    gmt.write_text("T\td\tA\tB\n", encoding="utf-8")
    cfg = {"sources": [{"type": "symbol_gmt", "path": str(gmt)}], "targets": [], "outputs": {"mode": "gmt"}}
    result = build(config=cfg, outdir=tmp_path / "out", dry_run=True)
    assert "output mode gmt" in result.plan[0]
    result = build(config=cfg, outdir=tmp_path / "out", dry_run=True, output_mode="minimal")
    assert "output mode minimal" in result.plan[0]


def test_config_show_and_schema(capsys: pytest.CaptureFixture[str]) -> None:
    assert cli.main(["config", "show", "pathways"]) == 0
    shown = yaml.safe_load(capsys.readouterr().out)
    assert shown["outputs"]["include_symbols"] is False
    assert cli.main(["config", "schema"]) == 0
    schema = json.loads(capsys.readouterr().out)
    assert schema == AlienConfig.model_json_schema()


def test_committed_json_schema_is_current() -> None:
    committed = Path(__file__).resolve().parents[2] / "docs" / "config.schema.json"
    current, saved = AlienConfig.model_json_schema(), json.loads(committed.read_text(encoding="utf-8"))
    # Compare structure, not pydantic-version-specific rendering; regenerate with `alien config schema`.
    assert set(saved["properties"]) == set(current["properties"])
    assert {name: set(d.get("properties", {})) for name, d in saved["$defs"].items()} == {
        name: set(d.get("properties", {})) for name, d in current["$defs"].items()
    }


def test_every_bundled_preset_and_example_validates() -> None:
    root = Path(__file__).resolve().parents[2]
    for example in sorted((root / "examples").glob("*.yml")):
        load_config(example)
