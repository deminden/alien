from alien.utils import strip_ensembl_version


def test_strip_ensembl_version():
    assert strip_ensembl_version("ENSG00000141510.18") == "ENSG00000141510"
    assert strip_ensembl_version("ENSG00000141510") == "ENSG00000141510"
    assert strip_ensembl_version("TP53") == "TP53"
    assert strip_ensembl_version("") == ""
    assert strip_ensembl_version(None) == ""
