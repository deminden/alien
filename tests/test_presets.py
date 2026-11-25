from pathlib import Path

from alien.config import list_presets, load_preset_config, preset_path


def test_packaged_preset_files_exist_for_all_examples():
    root = Path(__file__).resolve().parents[1]
    for preset in list_presets():
        assert preset_path(preset.name).exists()
        assert (root / "examples" / preset.filename).exists()


def test_bundled_presets_are_gtex_gencode47_only_without_external_output_gene_files():
    for name in ["pathways", "function_location", "disease_phenotype", "cancer_dependency"]:
        cfg = load_preset_config(name)
        assert [target.name for target in cfg.targets] == ["gtex_v11_gencode47"]
        gtex_target = cfg.targets[0]

        assert (gtex_target.annotation.source, gtex_target.annotation.version) == ("GENCODE", "47")
        assert gtex_target.annotation.path is None
        assert gtex_target.output_genes is None
        assert gtex_target.annotation.metadata_fallbacks == []


def test_examples_keep_study_specific_gtex_output_gene_paths():
    root = Path(__file__).resolve().parents[1]
    for name in ["pathways", "function_location", "disease_phenotype"]:
        text = (root / "examples" / f"{name}.yml").read_text(encoding="utf-8")

        assert "data/gtex/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_reads.parquet" in text


def test_preset_alias_loads_canonical_config():
    cfg = load_preset_config("cancer")

    assert cfg.project.outdir == Path("results/alien_outputs/cancer_dependency")
    assert cfg.targets[0].name == "gtex_v11_gencode47"
