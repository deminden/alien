from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from alien.cache import ResourceCache
from alien.config.schema import AnnotationSpec, TargetSpec
from alien.errors import ConfigError, ResourceError
from alien.resources.annotation import (
    gencode_annotation_url,
    load_annotation_file,
    normalize_gencode_version,
    prepare_target,
)

from conftest import FakeServer

GTF_LINE = 'chr1\tT\tgene\t1\t10\t.\t+\t.\tgene_id "{id}"; gene_name "{name}"; gene_type "protein_coding";\n'


def _gz(text: str) -> bytes:
    return gzip.compress(text.encode("utf-8"))


def test_gencode_version_downloads_without_explicit_path(tmp_path: Path, fake_http: FakeServer) -> None:
    fake_http.routes[gencode_annotation_url("49")] = _gz(GTF_LINE.format(id="ENSG00000141510.18", name="TP53"))
    spec = AnnotationSpec(source="GENCODE", version="49")
    frame, label, path = load_annotation_file("human", spec, "annotation", ResourceCache(tmp_path))
    assert label == "GENCODE v49"
    assert path == tmp_path / "gencode" / "gencode.v49.annotation.gtf.gz"
    assert frame["ensembl_gene_id"].tolist() == ["ENSG00000141510"]
    # The parsed gene table is cached by content hash; a second parse does not need the network.
    fake_http.routes.clear()
    again, _, _ = load_annotation_file("human", spec, "annotation", ResourceCache(tmp_path))
    assert again.equals(frame)
    assert list((tmp_path / "gtf_cache").glob("genes__*.parquet"))


def test_gencode_accepts_v_prefixed_version(tmp_path: Path, fake_http: FakeServer) -> None:
    url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz"
    fake_http.routes[url] = _gz(GTF_LINE.format(id="ENSG000002.1", name="GENE2"))
    _, label, path = load_annotation_file(
        "t", AnnotationSpec(source="GENCODE", version="v50"), "annotation", ResourceCache(tmp_path)
    )
    assert label == "GENCODE v50" and path.name == "gencode.v50.annotation.gtf.gz"
    assert normalize_gencode_version("V7") == "7"
    with pytest.raises(ConfigError, match="numeric"):
        normalize_gencode_version("latest")


def test_annotation_url_and_custom_attributes(tmp_path: Path, fake_http: FakeServer) -> None:
    fake_http.routes["https://example.test/custom.gtf?x=1"] = (
        'chr1\tT\tgene\t1\t10\t.\t+\t.\tID "ENSG00000141510.18"; Name "TP53"; biotype "protein_coding";\n'
        "chr1\tT\tgene\t20\t30\t.\t+\t.\tID=ENSG000002.1;Name=GENE2;biotype=lncRNA;\n"
    )
    spec = AnnotationSpec.model_validate(
        {
            "source": "Custom",
            "url": "https://example.test/custom.gtf?x=1",
            "attributes": {"gene_id": "ID", "gene_name": "Name", "gene_biotype": "biotype"},
        }
    )
    frame, label, path = load_annotation_file("custom", spec, "annotation", ResourceCache(tmp_path))
    assert label == "Custom" and path == tmp_path / "targets" / "custom" / "annotation" / "custom.gtf"
    assert frame["ensembl_gene_id"].tolist() == ["ENSG00000141510", "ENSG000002"]
    assert frame["gene_biotype"].tolist() == ["protein_coding", "lncRNA"]


def test_missing_local_annotation_is_an_error(tmp_path: Path) -> None:
    spec = AnnotationSpec(source="Local", path=tmp_path / "nope.gtf")
    with pytest.raises(ResourceError, match="does not exist"):
        load_annotation_file("t", spec, "annotation", ResourceCache(tmp_path))


def test_prepare_target_with_gene_filter(tmp_path: Path) -> None:
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(GTF_LINE.format(id="ENSG1.1", name="A") + GTF_LINE.format(id="ENSG2.1", name="B"), encoding="utf-8")
    spec = TargetSpec.model_validate(
        {"name": "t", "annotation": {"source": "X", "path": str(gtf)}, "gene_filter": ["ENSG2", "ENSG9"]}
    )
    target = prepare_target(spec, ResourceCache(tmp_path), [])
    assert target.effective_ids == {"ENSG2"}
    assert target.annotation.set_index("ensembl_gene_id")["is_in_output_genes"].to_dict() == {
        "ENSG1": False,
        "ENSG2": True,
    }
