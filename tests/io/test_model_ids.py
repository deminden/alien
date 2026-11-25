from __future__ import annotations

import math

from alien.model import MappingStatus, Target, clean_field, ensembl_stable_id, is_ensembl_gene_id, sanitize_term_id


def test_ensembl_stable_id() -> None:
    assert ensembl_stable_id("ENSG00000141510.18") == "ENSG00000141510"
    assert ensembl_stable_id("ENSG00000141510.3_PAR_Y") == "ENSG00000141510"
    assert ensembl_stable_id("ENSG00000141510") == "ENSG00000141510"
    assert ensembl_stable_id(" TP53 ") == "TP53"
    assert ensembl_stable_id("NA") == ""
    assert ensembl_stable_id(math.nan) == ""
    assert ensembl_stable_id(None) == ""
    assert is_ensembl_gene_id("ENSG1.2") and not is_ensembl_gene_id("TP53")


def test_sanitize_and_clean() -> None:
    assert sanitize_term_id("  A b/c  ") == "A_b_c"
    assert sanitize_term_id("***") == "UNNAMED"
    assert clean_field("a\tb\r\nc") == "a b c"
    assert clean_field(None) == ""


def test_mapping_status_properties() -> None:
    assert MappingStatus.ARCHIVE_ENSEMBL_RESCUE.is_mapped and MappingStatus.ARCHIVE_ENSEMBL_RESCUE.is_rescue
    assert MappingStatus.HGNC_CURRENT_SYMBOL.is_mapped and not MappingStatus.HGNC_CURRENT_SYMBOL.is_rescue
    assert not MappingStatus.AMBIGUOUS_DROP.is_mapped


def test_target_effective_ids() -> None:
    from conftest import make_target

    target = make_target(
        "t", [{"gene_symbol": "A", "ensembl_gene_id": "E1"}, {"gene_symbol": "B", "ensembl_gene_id": "E2"}]
    )
    assert target.effective_ids == {"E1", "E2"} and target.output_gene_source == "annotation"
    restricted = make_target("t", [{"gene_symbol": "A", "ensembl_gene_id": "E1"}], {"E1", "E9"})
    assert restricted.effective_ids == {"E1", "E9"} and restricted.output_gene_source == "output_genes"
    assert isinstance(restricted, Target)
