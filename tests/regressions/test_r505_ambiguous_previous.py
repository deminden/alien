"""5.5: a previous symbol shared by several genes is reported as ambiguous."""

from __future__ import annotations

from alien import build

from conftest import MiniWorld, read_tsv


def test_ambiguous_previous_symbol_is_audited(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"T1": ["PREVDUP", "G4", "G5"]})
    result = build(config=mini_world.config(gmt))
    ambiguous = read_tsv(result.outdir / "metadata" / "ambiguous_gene_mappings.tsv.gz")
    hits = ambiguous[ambiguous["input_gene"].eq("PREVDUP")]
    assert set(hits["target_namespace"]) == {"symbols", "mini"}
    assert set(hits["mapping_stage"]) == {"hgnc_previous_symbol"}
    assert set(hits["possible_symbols"]) == {"G1,G2"}
