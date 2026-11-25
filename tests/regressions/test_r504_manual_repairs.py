"""5.4: manual symbol repairs apply to every namespace and win over alias/previous lookups."""

from __future__ import annotations

from alien import build

from conftest import MiniWorld, read_gmt


def test_manual_repair_applies_to_symbols_namespace(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"T1": ["BROKEN-G4", "G2", "G3"]})
    cfg = mini_world.config(gmt, gene_mapping={"manual_symbol_repairs": {"BROKEN-G4": "G4"}})
    result = build(config=cfg)
    symbols = read_gmt(result.outdir / "gmt" / "symbols.gmt")["LOCAL__T1"]
    assert "G4" in symbols
    assert "BROKEN-G4" not in symbols


def test_manual_repair_overrides_alias(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"T1": ["WRONGALIAS", "G2", "G5"]})
    cfg = mini_world.config(gmt, gene_mapping={"manual_symbol_repairs": {"WRONGALIAS": "G4"}})
    result = build(config=cfg)
    assert "G4" in read_gmt(result.outdir / "gmt" / "symbols.gmt")["LOCAL__T1"]
    assert mini_world.ensembl("G4") in read_gmt(result.outdir / "gmt" / "mini.gmt")["LOCAL__T1"]


def test_manual_repair_does_not_override_exact_current_symbol(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"T1": ["G1", "G2", "G5"]})
    cfg = mini_world.config(gmt, gene_mapping={"manual_symbol_repairs": {"G1": "G4"}})
    result = build(config=cfg)
    assert "G1" in read_gmt(result.outdir / "gmt" / "symbols.gmt")["LOCAL__T1"]
