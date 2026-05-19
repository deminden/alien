import importlib
from pathlib import Path

from alien import BuildResult, build as run_build


MINIMAL_OUTPUTS = {
    "gmt/symbols.gmt",
    "metadata/effective_config.yml",
    "metadata/source_manifest.tsv",
    "metadata/source_provenance.json",
    "metadata/term_manifest.tsv.gz",
    "metadata/target_metadata_fallbacks.tsv",
    "metadata/target_output_genes.tsv.gz",
    "metadata/target_gene_filter.tsv.gz",
    "qc/target_namespace_summary.tsv",
    "qc/mapping_summary.tsv",
    "qc/redundancy_summary.tsv",
    "qc/warnings.txt",
}


def test_overlay_precedence_before_build(tmp_path, monkeypatch):
    build_module = importlib.import_module("alien.build")
    base = tmp_path / "base.yml"
    overlay1 = tmp_path / "overlay1.yml"
    overlay2 = tmp_path / "overlay2.yml"
    outdir = tmp_path / "cli-out"
    source_dir = tmp_path / "cli-source"
    captured = {}
    base.write_text(
        """
project:
  outdir: base-out
  source_dir: base-source
runtime:
  workers: 1
sources:
  - type: symbol_gmt
    path: base.gmt
targets:
  - name: human_test
    type: ensembl_gtf
    annotation:
      source: GENCODE
      version: "47"
""",
        encoding="utf-8",
    )
    overlay1.write_text(
        """
filtering:
  min_size_default: 2
sources:
  - type: symbol_gmt
    path: overlay1.gmt
""",
        encoding="utf-8",
    )
    overlay2.write_text(
        """
filtering:
  max_size_default: 77
sources:
  - type: symbol_gmt
    path: overlay2.gmt
""",
        encoding="utf-8",
    )

    def fake_build_to_output_dir(cfg, output_dir, source_dir, worker_count, dry_run, force_download):
        captured.update({"cfg": cfg, "output_dir": output_dir, "source_dir": source_dir, "worker_count": worker_count})
        return BuildResult(outdir=output_dir, namespaces=("symbols",), n_source_memberships=1, n_source_terms=1, warnings=())

    monkeypatch.setattr(build_module, "_build_to_output_dir", fake_build_to_output_dir)

    run_build(config=base, overrides=[overlay1, overlay2], outdir=outdir, source_dir=source_dir, workers=8)

    assert captured["cfg"]["project"]["outdir"] == str(outdir)
    assert captured["cfg"]["project"]["source_dir"] == str(source_dir)
    assert captured["cfg"]["runtime"]["workers"] == 8
    assert captured["cfg"]["filtering"]["min_size_default"] == 2
    assert captured["cfg"]["filtering"]["max_size_default"] == 77
    assert captured["cfg"]["sources"] == [{"type": "symbol_gmt", "path": "overlay2.gmt"}]
    assert captured["output_dir"] == outdir
    assert captured["source_dir"] == source_dir
    assert captured["worker_count"] == 8


def test_output_mode_gmt_copies_only_gmt_files(tmp_path, monkeypatch):
    build_module = importlib.import_module("alien.build")
    outdir = tmp_path / "final"
    monkeypatch.setattr(build_module, "_build_to_output_dir", _fake_staged_build)

    result = run_build(outdir=outdir, output_mode="gmt")

    assert result.outdir == outdir
    assert _files_under(outdir) == {"gmt/symbols.gmt"}


def test_output_mode_minimal_copies_reproducibility_pack(tmp_path, monkeypatch):
    build_module = importlib.import_module("alien.build")
    outdir = tmp_path / "final"
    monkeypatch.setattr(build_module, "_build_to_output_dir", _fake_staged_build)

    result = run_build(outdir=outdir, output_mode="minimal")

    assert result.outdir == outdir
    assert _files_under(outdir) == MINIMAL_OUTPUTS


def _fake_staged_build(cfg, output_dir, source_dir, worker_count, dry_run, force_download):
    assert not dry_run
    for relative_path in MINIMAL_OUTPUTS | {"metadata/gene_symbol_master.tsv.gz", "qc/extra.tsv"}:
        path = output_dir / relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(relative_path, encoding="utf-8")
    return BuildResult(outdir=output_dir, namespaces=("symbols",), n_source_memberships=1, n_source_terms=1, warnings=())


def _files_under(path: Path) -> set[str]:
    return {file.relative_to(path).as_posix() for file in path.rglob("*") if file.is_file()}
