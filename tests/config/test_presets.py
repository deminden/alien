from pathlib import Path

from alien.config import list_presets, load_preset_config, preset_path


def test_packaged_preset_files_exist_for_all_examples():
    root = Path(__file__).resolve().parents[2]
    for preset in list_presets():
        assert preset_path(preset.name).exists()
        assert (root / "examples" / preset.filename).exists()


def test_bundled_presets_are_gtex_gencode47_only_without_external_output_gene_files():
    for name in ["pathways", "function_location", "disease_phenotype", "cancer_programs", "cell_line_signatures"]:
        cfg = load_preset_config(name)
        assert [target.name for target in cfg.targets] == ["gtex_v11_gencode47"]
        gtex_target = cfg.targets[0]

        assert (gtex_target.annotation.source, gtex_target.annotation.version) == ("GENCODE", "47")
        assert gtex_target.annotation.path is None
        assert gtex_target.output_genes is None
        assert gtex_target.annotation.metadata_fallbacks == []


def test_examples_keep_study_specific_gtex_output_gene_paths():
    root = Path(__file__).resolve().parents[2]
    for name in ["pathways", "function_location", "disease_phenotype"]:
        text = (root / "examples" / f"{name}.yml").read_text(encoding="utf-8")

        assert "data/gtex/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_reads.parquet" in text


def test_preset_alias_loads_canonical_config():
    cfg = load_preset_config("cancer")

    assert cfg.project.outdir == Path("results/alien_outputs/cancer_programs")
    assert cfg.targets[0].name == "gtex_v11_gencode47"
    assert load_preset_config("cancer_dependency").project.outdir == cfg.project.outdir


def _tags(cfg) -> set[str]:
    tags: set[str] = set()
    for source in cfg.sources:
        libraries = getattr(source, "libraries", None)
        if libraries:
            tags.update(library.source_tag for library in libraries)
        else:
            tags.add(f"{source.collection}:{source.subcollection or ''}")
    return tags


def test_cell_line_libraries_are_opt_in_only():
    """Per-cell-line libraries and C4:CGN live only in the cell_line_signatures preset."""
    cell_line_tags = {"CCLE", "NCI60", "DEPMAP"}
    for name in ["pathways", "function_location", "disease_phenotype", "cancer_programs"]:
        tags = _tags(load_preset_config(name))
        assert not tags & cell_line_tags, name
        assert "C4:CGN" not in tags, name
    assert _tags(load_preset_config("cell_lines")) == cell_line_tags
    assert _tags(load_preset_config("cancer")) == {"C6:", "C4:3CA", "C9:"}
