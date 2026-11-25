"""5.1: the term manifest records removal per (term, namespace)."""

from __future__ import annotations

from conftest import MiniWorld, read_tsv

from alien import build


def test_term_manifest_is_long_format_per_namespace(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"T1": ["G1", "G2", "G3"], "T2": ["G1", "G2", "G3"]})
    targets = [mini_world.target("nsA"), mini_world.target("nsB", output_genes=["G1", "G4"])]
    cfg = mini_world.config(gmt, targets=targets, outputs={"include_symbols": False})
    result = build(config=cfg)

    manifest = read_tsv(result.outdir / "metadata" / "term_manifest.tsv.gz")
    rows = {(row.term_id, row.target_namespace): row for row in manifest.itertuples()}
    assert set(rows) == {("LOCAL__T1", "nsA"), ("LOCAL__T2", "nsA"), ("LOCAL__T1", "nsB"), ("LOCAL__T2", "nsB")}

    t1_a = rows[("LOCAL__T1", "nsA")]
    assert (t1_a.kept, t1_a.removed_reason, t1_a.size_after_mapping) == ("True", "", "3")
    t2_a = rows[("LOCAL__T2", "nsA")]
    assert (t2_a.kept, t2_a.removed_reason, t2_a.representative_term_id) == ("False", "exact_duplicate", "LOCAL__T1")
    t1_b = rows[("LOCAL__T1", "nsB")]
    assert (t1_b.kept, t1_b.removed_reason, t1_b.size_after_mapping) == ("False", "too_small", "1")
    assert t1_b.representative_term_id == "LOCAL__T1"
