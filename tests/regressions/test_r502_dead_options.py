"""5.2: documented options either work or are rejected."""

from __future__ import annotations

import pytest

from alien import build

from conftest import MiniWorld, read_gmt, read_tsv


def test_exact_duplicate_removal_can_be_disabled(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"T1": ["G1", "G2", "G3"], "T2": ["G1", "G2", "G3"]})
    cfg = mini_world.config(gmt, redundancy={"exact_duplicate_removal": False, "jaccard_cutoff": 1.0})
    result = build(config=cfg)
    assert set(read_gmt(result.outdir / "gmt" / "mini.gmt")) == {"LOCAL__T1", "LOCAL__T2"}


def test_apply_per_family_false_compares_across_families(mini_world: MiniWorld) -> None:
    gmt_a = mini_world.gmt("a.gmt", {"T1": ["G1", "G2", "G3"]})
    gmt_b = mini_world.gmt("b.gmt", {"T2": ["G1", "G2", "G3"]})
    cfg = mini_world.config(gmt_a)
    cfg["sources"].append(
        {"type": "symbol_gmt", "path": str(gmt_b), "source_tag": "DIS", "family": "disease_phenotype"}
    )
    per_family = build(config=cfg, outdir=mini_world.root / "per_family")
    assert set(read_gmt(per_family.outdir / "gmt" / "mini.gmt")) == {"LOCAL__T1", "DIS__T2"}

    cfg["redundancy"] = {"apply_per_family": False}
    across = build(config=cfg, outdir=mini_world.root / "across")
    assert len(read_gmt(across.outdir / "gmt" / "mini.gmt")) == 1
    across_parallel = build(config=cfg, outdir=mini_world.root / "across_parallel", workers=2)
    assert read_gmt(across_parallel.outdir / "gmt" / "mini.gmt") == read_gmt(across.outdir / "gmt" / "mini.gmt")


def test_use_hgnc_alias_mapping_false_skips_aliases(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"T1": ["ALIASG1", "G2", "G3"]})
    cfg = mini_world.config(gmt, gene_mapping={"use_hgnc_alias_mapping": False})
    result = build(config=cfg)
    genes = read_gmt(result.outdir / "gmt" / "mini.gmt")["LOCAL__T1"]
    assert mini_world.ensembl("G1") not in genes
    unmapped = read_tsv(result.outdir / "metadata" / "unmapped_genes.tsv.gz")
    assert "ALIASG1" in set(unmapped["input_gene"])


@pytest.mark.parametrize(
    ("section", "key"),
    [
        ("redundancy", "apply_after_output_gene_projection"),
        ("redundancy", "nested_overlap_filter"),
        ("filtering", "remove_empty_after_output_gene_projection"),
        ("gene_mapping", "strip_ensembl_version"),
        ("gene_mapping", "allow_ambiguous_alias_mapping"),
        ("gene_mapping", "prefer_gencode_exact_gene_name"),
        ("gene_mapping", "prefer_current_hgnc_symbol"),
        ("gene_mapping", "drop_unmapped_from_ensg_gmt"),
        ("ensembl_archive", "unresolved_example_threshold"),
    ],
)
def test_removed_options_are_rejected(mini_world: MiniWorld, section: str, key: str) -> None:
    gmt = mini_world.gmt("lib.gmt", {"T1": ["G1", "G2"]})
    cfg = mini_world.config(gmt, **{section: {key: True}})
    with pytest.raises(ValueError, match=key):
        build(config=cfg)
