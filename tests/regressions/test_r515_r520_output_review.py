"""Fixes from the 0.2.0 output review: term names, member weights, letter case and the output-directory error."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from alien import build, cli
from alien.errors import OutputDirectoryError, TermIdCollisionError
from alien.io.gmt import description_from_meta, member_symbol, parse_gmt
from alien.mapping import resolve_symbol

from conftest import MiniWorld, hgnc_maps, read_gmt, read_tsv


def test_original_name_keeps_the_source_spelling(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"Weighted set (HP:0001)": ["G1", "G2", "G3"]})
    result = build(config=mini_world.config(gmt))
    manifest = read_tsv(result.outdir / "metadata" / "term_manifest.tsv.gz")
    assert set(manifest["term_id"]) == {"LOCAL__Weighted_set_HP_0001"}
    assert set(manifest["original_name"]) == {"Weighted set (HP:0001)"}
    line = (result.outdir / "gmt" / "mini.gmt").read_text(encoding="utf-8")
    assert "original_name=Weighted set (HP:0001)\t" in line


@pytest.mark.parametrize(
    "names",
    [("TGF beta/SMAD", "TGF beta SMAD"), ("\u03b1-thalassemia", "\u03b2-thalassemia")],
    ids=["punct", "greek"],
)
def test_names_that_sanitize_to_one_term_id_are_a_collision(mini_world: MiniWorld, names: tuple[str, str]) -> None:
    gmt = mini_world.gmt("lib.gmt", {names[0]: ["G1", "G2"], names[1]: ["G3", "G4", "G5"]})
    with pytest.raises(TermIdCollisionError):
        build(config=mini_world.config(gmt, term_id_collisions={"action": "error"}))
    audit = read_tsv(mini_world.root / "out" / "failed" / "metadata" / "term_id_collisions.tsv")
    assert set(audit["original_name"]) == set(names)


def test_colliding_names_merge_only_when_asked_and_are_audited(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"TGF beta/SMAD": ["G1", "G2"], "TGF beta SMAD": ["G3", "G4", "G5"]})
    result = build(config=mini_world.config(gmt, term_id_collisions={"action": "merge"}))
    assert len(read_gmt(result.outdir / "gmt" / "mini.gmt")["LOCAL__TGF_beta_SMAD"]) == 5
    assert (result.outdir / "metadata" / "term_id_collisions.tsv").exists()
    assert any("conflicting term/source identities" in warning for warning in result.warnings)


def test_multi_word_broad_disease_patterns_match_gmt_terms(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"Malignant neoplasm": ["G1", "G2", "G3"], "Rare syndrome": ["G3", "G4", "G5"]})
    cfg = mini_world.config(gmt, filtering={"broad_disease_regex": ["^malignant neoplasm$"]})
    cfg["sources"][0]["family"] = "disease_phenotype"
    result = build(config=cfg)
    assert set(read_gmt(result.outdir / "gmt" / "mini.gmt")) == {"LOCAL__Rare_syndrome"}
    removed = read_tsv(result.outdir / "metadata" / "removed_terms_size_filter.tsv.gz")
    assert set(removed["reason"]) == {"broad_disease_term"}


def test_repeated_term_name_in_one_gmt_is_logged(mini_world: MiniWorld, caplog: pytest.LogCaptureFixture) -> None:
    path = mini_world.root / "lib.gmt"
    path.write_text("SAME\td\tG1\tG2\nSAME\td\tG3\n", encoding="utf-8")
    with caplog.at_level("WARNING", logger="alien"):
        result = build(config=mini_world.config(path))
    assert len(read_gmt(result.outdir / "gmt" / "mini.gmt")["LOCAL__SAME"]) == 3
    assert "more than once" in caplog.text


@pytest.mark.parametrize(
    ("token", "symbol"),
    [
        ("TP53,1.0", "TP53"),
        (" TP53 , 0.25 ", "TP53"),
        ("TP53,1e-3", "TP53"),
        ("TP53,5", "TP53"),
        ("HLA-DRB1", "HLA-DRB1"),
    ],
)
def test_member_weights_are_stripped(token: str, symbol: str) -> None:
    assert member_symbol(token) == symbol


def test_weighted_gmt_members_map(mini_world: MiniWorld, tmp_path: Path) -> None:
    gmt = mini_world.gmt("lib.gmt", {"FUZZY": ["G1,1.0", "G2,0.5", "G3,0.2"]})
    assert parse_gmt(gmt)[0].genes == ("G1", "G2", "G3")
    result = build(config=mini_world.config(gmt))
    assert len(read_gmt(result.outdir / "gmt" / "mini.gmt")["LOCAL__FUZZY"]) == 3
    assert read_tsv(result.outdir / "metadata" / "unmapped_genes.tsv.gz").empty


def test_description_values_cannot_break_the_key_value_format() -> None:
    assert description_from_meta({"source_tag": "S", "original_name": "a; b=c"}) == "source=S;original_name=a, b=c"


HGNC = [  # symbol, previous symbols, alias symbols, Ensembl ID
    ("C1orf112", "", "", "ENSG01"),
    ("FIRRM", "C1orf999", "Apolo1", "ENSG02"),
    ("TP53", "", "p53", "ENSG03"),
    ("AAA", "", "Shared", "ENSG04"),
    ("BBB", "", "SHARED", "ENSG05"),
    ("Abc1", "", "", "ENSG06"),
    ("ABC1", "", "", "ENSG07"),
]


def test_case_insensitive_fallback_order_and_status() -> None:
    hgnc = hgnc_maps(HGNC)
    assert resolve_symbol("C1ORF112", hgnc).status == "unmapped"  # off unless requested
    current = resolve_symbol("C1ORF112", hgnc, case_insensitive=True)
    assert (current.symbol, current.status) == ("C1orf112", "hgnc_case_insensitive_symbol")
    assert resolve_symbol("C1ORF999", hgnc, case_insensitive=True).symbol == "FIRRM"  # previous symbol
    assert resolve_symbol("APOLO1", hgnc, case_insensitive=True).symbol == "FIRRM"  # alias
    assert resolve_symbol("APOLO1", hgnc, case_insensitive=True, use_alias=False).status == "unmapped"
    assert resolve_symbol("tp53", hgnc, case_insensitive=True).symbol == "TP53"


def test_exact_matches_win_and_folded_ties_are_ambiguous() -> None:
    hgnc = hgnc_maps(HGNC)
    assert resolve_symbol("p53", hgnc, case_insensitive=True).status == "hgnc_alias_symbol"
    assert resolve_symbol("Abc1", hgnc, case_insensitive=True).status == "exact_current_symbol"
    tie = resolve_symbol("aBc1", hgnc, case_insensitive=True)
    assert (tie.symbol, tie.status, tie.candidates) == ("", "ambiguous_case_insensitive", ("ABC1", "Abc1"))
    assert resolve_symbol("shared", hgnc, case_insensitive=True).candidates == ("AAA", "BBB")


def test_case_variants_map_in_a_build_and_are_audited(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"T1": ["g1", "G2", "G3"]})
    result = build(config=mini_world.config(gmt))
    assert mini_world.ensembl("G1") in read_gmt(result.outdir / "gmt" / "mini.gmt")["LOCAL__T1"]
    assert "G1" in read_gmt(result.outdir / "gmt" / "symbols.gmt")["LOCAL__T1"]
    status = read_tsv(result.outdir / "qc" / "mapping_status_summary.tsv")
    assert "hgnc_case_insensitive_symbol" in set(status["mapping_status"])
    strict = build(config=mini_world.config(gmt, gene_mapping={"case_insensitive_symbols": False}))
    assert mini_world.ensembl("G1") not in read_gmt(strict.outdir / "gmt" / "mini.gmt")["LOCAL__T1"]


def test_refusing_a_foreign_output_directory_is_a_usage_error(
    mini_world: MiniWorld, capsys: pytest.CaptureFixture[str]
) -> None:
    gmt = mini_world.gmt("lib.gmt", {"T1": ["G1", "G2", "G3"]})
    config = mini_world.root / "config.yml"
    config.write_text(yaml.safe_dump(mini_world.config(gmt)), encoding="utf-8")
    foreign = mini_world.root / "notes"
    foreign.mkdir()
    (foreign / "notes.txt").write_text("mine", encoding="utf-8")
    with pytest.raises(OutputDirectoryError):
        build(config=config, outdir=foreign)
    assert cli.main(["build", str(config), "--outdir", str(foreign)]) == 2
    assert "unexpected error" not in capsys.readouterr().err
    assert (foreign / "notes.txt").read_text(encoding="utf-8") == "mine"
