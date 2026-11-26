"""List-quality changes from the 2026-11-26 real-data review.

Collisions are renamed instead of failing, curated exclusions and broad categories are removed with an audit row,
short text-mined aliases, legacy mitochondrial names, spreadsheet dates and ``LOC<GeneID>`` tokens are handled, and
every non-exact mapping is listed once in ``metadata/symbol_rescues``.
"""

from __future__ import annotations

from typing import Any

import pandas as pd
import pytest

from alien import build
from alien.config import load_config
from alien.mapping import map_namespaces, resolve_symbol
from alien.mapping.symbols import MITOCHONDRIAL_LEGACY, SymbolOptions
from alien.resources.ncbi import NcbiMaps, NcbiRecord

from conftest import MiniWorld, hgnc_maps, make_target, membership, read_gmt, read_tsv

# ------------------------------------------------------------------------------------------- collisions


def test_colliding_names_are_kept_apart_by_default(mini_world: MiniWorld) -> None:
    """DisGeNET has 'beta Thalassemia' and 'beta^+^ Thalassemia': two diseases, one sanitized ID."""
    gmt = mini_world.gmt("lib.gmt", {"beta Thalassemia": ["G1", "G2"], "beta^+^ Thalassemia": ["G3", "G4", "G5"]})
    result = build(config=mini_world.config(gmt))
    terms = read_gmt(result.outdir / "gmt" / "symbols.gmt")
    renamed = [t for t in terms if t.startswith("LOCAL__beta_Thalassemia__")]
    assert sorted(terms["LOCAL__beta_Thalassemia"]) == ["G1", "G2"]
    assert len(renamed) == 1 and sorted(terms[renamed[0]]) == ["G3", "G4", "G5"]
    audit = read_tsv(result.outdir / "metadata" / "term_id_collisions.tsv")
    assert set(audit["action"]) == {"rename"}
    assert set(audit["resolved_term_id"]) == {"LOCAL__beta_Thalassemia", renamed[0]}
    manifest = read_tsv(result.outdir / "metadata" / "term_manifest.tsv.gz")
    assert set(manifest.loc[manifest["term_id"] == renamed[0], "original_name"]) == {"beta^+^ Thalassemia"}
    assert not result.warnings


def test_renamed_ids_are_stable_across_builds(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"TGF beta/SMAD": ["G1", "G2"], "TGF beta SMAD": ["G3", "G4", "G5"]})
    first = set(read_gmt(build(config=mini_world.config(gmt)).outdir / "gmt" / "mini.gmt"))
    second = set(read_gmt(build(config=mini_world.config(gmt)).outdir / "gmt" / "mini.gmt"))
    assert first == second and len(first) == 2


# ------------------------------------------------------------------------------------------- term filters


def _removed(result: Any) -> pd.DataFrame:
    return read_tsv(result.outdir / "metadata" / "removed_terms_size_filter.tsv.gz")


def test_exclusion_rules_remove_terms_with_the_rule_name(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt(
        "lib.gmt",
        {
            "KEGG_MEDICUS_REFERENCE_EGF_EGFR_RAS_ERK": ["G1", "G2"],
            "KEGG_MEDICUS_VARIANT_MUTATION_ACTIVATED_KRAS": ["G1", "G2", "G3"],
            "KEGG_MEDICUS_PATHOGEN_HPV_E5_TO_EGFR": ["G1", "G4"],
        },
    )
    cfg = mini_world.config(gmt)
    cfg["sources"][0]["source_tag"] = "KEGG_MEDICUS"
    result = build(config=cfg)
    kept = set(read_gmt(result.outdir / "gmt" / "mini.gmt"))
    assert kept == {"KEGG_MEDICUS__KEGG_MEDICUS_REFERENCE_EGF_EGFR_RAS_ERK"}
    removed = _removed(result)
    removed = removed[removed["target_namespace"] == "mini"]
    assert set(removed["reason"]) == {"excluded_term"}
    assert set(removed["exclusion_rule"]) == {"kegg_medicus_variant_and_pathogen"}
    manifest = read_tsv(result.outdir / "metadata" / "term_manifest.tsv.gz")
    assert set(manifest.loc[manifest["kept"] == "False", "removed_reason"]) == {"excluded_term"}


def test_exclusion_rules_only_apply_to_their_source_tags_and_can_be_switched_off(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"KEGG_MEDICUS_VARIANT_X": ["G1", "G2"], "SARS-CoV-2 infection": ["G3", "G4"]})
    assert len(read_gmt(build(config=mini_world.config(gmt)).outdir / "gmt" / "mini.gmt")) == 2  # tag LOCAL
    cfg = mini_world.config(gmt, filtering={"exclude_terms": {"pathogen_specific_pathways": {"enabled": False}}})
    cfg["sources"][0]["source_tag"] = "REACTOME"
    assert len(read_gmt(build(config=cfg).outdir / "gmt" / "mini.gmt")) == 2
    cfg["filtering"]["exclude_terms"] = {"pathogen_specific_pathways": {"enabled": True}}
    assert set(read_gmt(build(config=cfg).outdir / "gmt" / "mini.gmt")) == {"REACTOME__KEGG_MEDICUS_VARIANT_X"}


def test_custom_exclusion_rule_merges_with_the_defaults() -> None:
    cfg = load_config(
        {
            "sources": [{"type": "symbol_gmt", "path": "x.gmt"}],
            "filtering": {"exclude_terms": {"mine": {"patterns": ["^OLD "]}}},
        }
    )
    assert {"mine", "hpo_modifiers", "kegg_medicus_variant_and_pathogen"} <= set(cfg.filtering.exclude_terms)


@pytest.mark.parametrize(
    ("name", "excluded"),
    [
        ("HP_X_LINKED_RECESSIVE_INHERITANCE", True),
        ("HP_ADULT_ONSET", True),
        ("HP_AGE_OF_DEATH", True),
        ("HP_TYPIFIED_BY_INCOMPLETE_PENETRANCE", True),
        ("HP_BILATERAL_TONIC_CLONIC_SEIZURE_WITH_FOCAL_ONSET", False),
        ("HP_SUDDEN_CARDIAC_DEATH", False),
        ("HP_ATRIAL_FIBRILLATION", False),
    ],
)
def test_hpo_modifier_terms(mini_world: MiniWorld, name: str, excluded: bool) -> None:
    gmt = mini_world.gmt("lib.gmt", {name: ["G1", "G2"]})
    cfg = mini_world.config(gmt)
    cfg["sources"][0].update({"source_tag": "HPO", "family": "disease_phenotype"})
    terms = read_gmt(build(config=cfg).outdir / "gmt" / "mini.gmt")
    assert (not terms) is excluded


@pytest.mark.parametrize(
    "name",
    [
        "IMMUNE SYSTEM DISEASE",
        "Immune System Diseases",
        "DISEASE OF CELLULAR PROLIFERATION",
        "SYNDROME",
        "Neoplasm Metastasis",
        "HP_ABNORMALITY_OF_THE_EYE",
        "ABNORMALITY OF THE NERVOUS SYSTEM",
    ],
)
def test_broad_disease_categories(mini_world: MiniWorld, name: str) -> None:
    gmt = mini_world.gmt("lib.gmt", {name: ["G1", "G2"], "Breast cancer": ["G3", "G4"]})
    cfg = mini_world.config(gmt)
    cfg["sources"][0]["family"] = "disease_phenotype"
    result = build(config=cfg)
    assert set(read_gmt(result.outdir / "gmt" / "mini.gmt")) == {"LOCAL__Breast_cancer"}
    assert "broad_disease_term" in set(_removed(result)["reason"])


# ------------------------------------------------------------------------------------------- symbols

HGNC = [
    ("MT-ND1", "", "ND1", "ENSG_MTND1"),
    ("IVNS1ABP", "", "ND1", "ENSG_IVNS"),
    ("MT-TP", "", "", "ENSG_MTTP"),
    ("TRNP1", "", "TRNP", "ENSG_TRNP1"),
    ("HYCC1", "", "HCC", "ENSG_HYCC1"),
    ("KMT2A", "", "HRX|ALL-1", "ENSG_KMT2A"),
    ("SEPTIN1", "SEPT1", "", "ENSG_SEPT1"),
    ("MARCHF6", "", "", "ENSG_MARCH6"),
    ("MARCHF1", "", "", "ENSG_MARCH1"),
    ("MTARC1", "", "", "ENSG_MTARC1"),
    ("PTGS2", "", "COX2", "ENSG_PTGS2"),
    ("MT-CO2", "", "COX2", "ENSG_MTCO2"),
]
DEFAULTS = SymbolOptions(min_alias_length=4, mitochondrial_legacy=True, excel_dates=True, case_insensitive=True)


@pytest.mark.parametrize(
    ("token", "symbol", "status"),
    [
        ("ND1", "MT-ND1", "mitochondrial_legacy_symbol"),  # alias of IVNS1ABP too; the MT name wins
        ("TRNP", "MT-TP", "mitochondrial_legacy_symbol"),  # HGNC alias of TRNP1
        ("COX2", "", "ambiguous_alias"),  # PTGS2 or MT-CO2: left ambiguous
        ("HCC", "", "short_alias_rejected"),  # hepatocellular carcinoma, not HYCC1
        ("hcc", "", "short_alias_rejected"),
        ("ALL-1", "KMT2A", "hgnc_alias_symbol"),
        ("1-Sep", "SEPTIN1", "excel_date_symbol"),
        ("SEP-1", "SEPTIN1", "excel_date_symbol"),
        ("6-MAR", "MARCHF6", "excel_date_symbol"),
        ("1-Mar", "", "ambiguous_excel_date"),  # MARCHF1 or MTARC1
        ("13-Sep", "", "unmapped"),  # no such gene in this table
    ],
)
def test_symbol_rules(token: str, symbol: str, status: str) -> None:
    resolution = resolve_symbol(token, hgnc_maps(HGNC), options=DEFAULTS)
    assert (resolution.symbol, resolution.status) == (symbol, status)


def test_symbol_rules_are_off_with_legacy_options() -> None:
    hgnc = hgnc_maps(HGNC)
    assert resolve_symbol("HCC", hgnc).symbol == "HYCC1"
    assert resolve_symbol("TRNP", hgnc).symbol == "TRNP1"
    assert resolve_symbol("1-Sep", hgnc).status == "unmapped"


def test_mitochondrial_table_targets_mt_genes_only() -> None:
    assert all(value.startswith("MT-") for value in MITOCHONDRIAL_LEGACY.values())
    assert "COX1" not in MITOCHONDRIAL_LEGACY and "COX2" not in MITOCHONDRIAL_LEGACY


def test_short_alias_is_reported_as_unmapped_with_its_candidate(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"T1": ["G4", "G5", "AG1"]})
    cfg = mini_world.config(gmt)
    (mini_world.source_dir / "hgnc" / "hgnc_complete_set.txt").write_text(
        "symbol\tprev_symbol\talias_symbol\tensembl_gene_id\n"
        + "".join(
            f"{s}\t\t{'AG1' if s == 'G1' else ''}\t{e}\n"
            for s, e in [
                ("G1", mini_world.ensembl("G1")),
                ("G4", mini_world.ensembl("G4")),
                ("G5", mini_world.ensembl("G5")),
            ]
        ),
        encoding="utf-8",
    )
    # LOCAL is not a text-mined source: short aliases are accepted.
    assert mini_world.ensembl("G1") in read_gmt(build(config=cfg).outdir / "gmt" / "mini.gmt")["LOCAL__T1"]
    cfg["gene_mapping"] = {"short_alias_source_tags": ["LOCAL"]}
    result = build(config=cfg)
    assert mini_world.ensembl("G1") not in read_gmt(result.outdir / "gmt" / "mini.gmt")["LOCAL__T1"]
    assert "G1" not in read_gmt(result.outdir / "gmt" / "symbols.gmt")["LOCAL__T1"]
    unmapped = read_tsv(result.outdir / "metadata" / "unmapped_genes.tsv.gz")
    row = unmapped[(unmapped["target_namespace"] == "mini") & (unmapped["input_gene"] == "AG1")].iloc[0]
    assert (row["reason"], row["possible_matches"]) == ("short_alias_rejected", "G1")
    cfg["gene_mapping"] = {"short_alias_source_tags": ["*"], "min_alias_length": 3}
    assert mini_world.ensembl("G1") in read_gmt(build(config=cfg).outdir / "gmt" / "mini.gmt")["LOCAL__T1"]


# ------------------------------------------------------------------------------------------- LOC<GeneID>


def _loc_world() -> tuple[Any, Any, NcbiMaps]:
    hgnc = hgnc_maps([("GGTLC2", "", "", "ENSG00000237766"), ("KEEP", "", "", "ENSG00000000001")])
    target = make_target(
        "t",
        [
            {"gene_symbol": "GGTLC2", "ensembl_gene_id": "ENSG00000237766"},
            {"gene_symbol": "KEEP", "ensembl_gene_id": "ENSG00000000001"},
            {"gene_symbol": "LINC9", "ensembl_gene_id": "ENSG00000000009"},
        ],
    )
    ncbi = NcbiMaps(
        gene_ids={
            "92249": NcbiRecord("GGTLC2", ("92249",), ("ENSG00000237766",)),
            "999": NcbiRecord("LINC9", ("999",), ("ENSG00000000009",)),
        },
        replaced_ids={"555": "999", "777": "-"},
    )
    return hgnc, target, ncbi


def test_loc_tokens_map_through_their_gene_id_for_every_source() -> None:
    hgnc, target, ncbi = _loc_world()
    rows = [
        membership("T1", gene, source_tag="DISGENET") for gene in ["KEEP", "LOC92249", "LOC555", "LOC777", "LOC123"]
    ]
    config = load_config({"sources": [{"type": "symbol_gmt", "path": "x.gmt"}], "outputs": {"include_symbols": False}})
    result = map_namespaces(pd.DataFrame(rows), config, [target], hgnc=hgnc, ncbi=ncbi, resolver=None)
    assert result.terms["t"]["T1"].genes == {"ENSG00000000001", "ENSG00000237766", "ENSG00000000009"}
    statuses = dict(zip(result.status_summary["mapping_status"], result.status_summary["n_gene_members"], strict=True))
    assert statuses["ncbi_gene_id"] == 2 and statuses["unmapped_drop"] == 2
    audit = result.ncbi_audit.set_index("input_gene")
    assert audit.loc["LOC555", "ncbi_current_symbol"] == "LINC9"
    assert audit.loc["LOC777", "ncbi_status"] == "ncbi_gene_id_withdrawn"
    rescues = result.symbol_rescues.set_index("input_gene")
    assert set(rescues.index) == {"LOC92249", "LOC555"}
    assert rescues.loc["LOC92249", "mapped_symbol"] == "GGTLC2"


# ------------------------------------------------------------------------------------------- rescue audit


def test_symbol_rescues_lists_each_non_exact_token_once(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"T1": ["OLDG1", "G2", "g3"], "T2": ["OLDG1", "G4", "G5"]})
    result = build(config=mini_world.config(gmt))
    rescues = read_tsv(result.outdir / "metadata" / "symbol_rescues.tsv.gz")
    rescues = rescues[rescues["target_namespace"] == "mini"].set_index("input_gene")
    assert set(rescues.index) == {"OLDG1", "g3"}
    assert rescues.loc["OLDG1", ["mapping_status", "mapped_symbol", "n_memberships", "n_terms"]].tolist() == [
        "hgnc_previous_symbol",
        "G1",
        "2",
        "2",
    ]
    assert rescues.loc["g3", "mapping_status"] == "hgnc_case_insensitive_symbol"
