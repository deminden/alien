"""Phase 6: lock file, fetch, --frozen and --offline builds, verify."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from alien import build, cli, fetch
from alien.cache import LockFile, ResourceCache, derived_cache_path, file_sha256
from alien.errors import ResourceError
from alien.outputs import verify_output_dir

FIXTURE = Path(__file__).resolve().parents[1] / "golden" / "fixture"


@pytest.fixture
def project(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    root = tmp_path / "project"
    shutil.copytree(FIXTURE, root)
    monkeypatch.chdir(tmp_path)
    return root


def test_fetch_writes_lock_and_frozen_build_succeeds_offline(project: Path) -> None:
    config = project / "config.yml"
    lock_path = fetch(config)
    assert lock_path == project / "alien.lock.json"
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    assert set(lock["resources"]) == {
        "hgnc",
        "ncbi_gene:gene_info",
        "ncbi_gene:gene_history",
        "ensembl_archive:lookups:fixture_annotation",
        "ensembl_archive:lookups:fixture_measured",
    }
    assert lock["resources"]["hgnc"]["sha256"] == file_sha256(project / "sources" / "hgnc" / "hgnc_complete_set.txt")
    result = build(config, frozen=True)
    provenance = json.loads((result.outdir / "metadata" / "source_provenance.json").read_text(encoding="utf-8"))
    assert provenance["cache_mode"] == "frozen"
    assert provenance["resources"]["hgnc"]["sha256"] == lock["resources"]["hgnc"]["sha256"]
    assert verify_output_dir(result.outdir) == []


def test_frozen_build_fails_when_a_cached_file_changes(project: Path) -> None:
    config = project / "config.yml"
    fetch(config)
    hgnc = project / "sources" / "hgnc" / "hgnc_complete_set.txt"
    hgnc.write_text(hgnc.read_text(encoding="utf-8") + "ZZZ1\t\t\tENSG00000000999\n", encoding="utf-8")
    with pytest.raises(ResourceError, match=r"differs from .*alien\.lock\.json"):
        build(config, frozen=True)
    assert cli.main(["build", str(config), "--frozen"]) == 4


def test_frozen_needs_a_lock_and_offline_needs_a_cache(project: Path) -> None:
    config = project / "config.yml"
    with pytest.raises(ResourceError, match="needs a lock file"):
        build(config, frozen=True)
    shutil.rmtree(project / "sources" / "hgnc")
    with pytest.raises(ResourceError, match=r"not cached and network access is disabled \(--offline\)"):
        build(config, offline=True)


def test_verify_detects_changes(project: Path) -> None:
    result = build(project / "config.yml", output_mode="minimal")
    assert verify_output_dir(result.outdir) == []
    gmt = result.outdir / "gmt" / "symbols.gmt"
    gmt.write_text(gmt.read_text(encoding="utf-8") + "X\tx\tA\n", encoding="utf-8")
    (result.outdir / "extra.txt").write_text("?", encoding="utf-8")
    (result.outdir / "qc" / "warnings.txt").unlink()
    assert verify_output_dir(result.outdir) == [
        "changed: gmt/symbols.gmt",
        "missing: qc/warnings.txt",
        "unexpected: extra.txt",
    ]
    assert cli.main(["verify", str(result.outdir)]) == 3
    assert verify_output_dir(project / "nowhere") == [f"{project / 'nowhere' / 'MANIFEST.json'} does not exist"]


def test_cli_fetch_and_lock_option(project: Path, tmp_path: Path) -> None:
    lock = tmp_path / "custom.lock.json"
    assert cli.main(["fetch", str(project / "config.yml"), "--lock", str(lock)]) == 0
    assert LockFile.load(lock).resources
    assert (
        cli.main(["build", str(project / "config.yml"), "--frozen", "--lock", str(lock), "--output-mode", "gmt"]) == 0
    )


def test_lock_file_format_errors(tmp_path: Path) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text('{"lock_version": 99}', encoding="utf-8")
    with pytest.raises(ResourceError, match="unsupported format"):
        LockFile.load(bad)
    bad.write_text("{", encoding="utf-8")
    with pytest.raises(ResourceError, match="unreadable"):
        LockFile.load(bad)
    with pytest.raises(ResourceError, match="No lock file path"):
        LockFile().save()


def test_derived_cache_path_depends_on_inputs(tmp_path: Path) -> None:
    a = derived_cache_path(tmp_path, "x", {"sha": "1", "opt": True})
    assert a == derived_cache_path(tmp_path, "x", {"opt": True, "sha": "1"})
    assert a != derived_cache_path(tmp_path, "x", {"sha": "1", "opt": False})


def test_file_sha256_memo_tracks_modification(tmp_path: Path) -> None:
    path = tmp_path / "f.txt"
    path.write_text("one", encoding="utf-8")
    first = file_sha256(path)
    assert file_sha256(path) == first
    path.write_text("two!", encoding="utf-8")
    assert file_sha256(path) != first


def test_concurrent_fetches_share_one_download(tmp_path: Path, fake_http: object) -> None:
    from concurrent.futures import ThreadPoolExecutor

    fake_http.routes["https://example.test/r.txt"] = b"payload"  # type: ignore[attr-defined]
    cache = ResourceCache(tmp_path)
    with ThreadPoolExecutor(4) as pool:
        paths = list(pool.map(lambda _: cache.fetch("r", "https://example.test/r.txt", tmp_path / "r.txt"), range(4)))
    assert {p.read_bytes() for p in paths} == {b"payload"}
    assert fake_http.urls() == ["https://example.test/r.txt"]  # type: ignore[attr-defined]
    assert not list(tmp_path.glob(".r.txt.*.part"))


def test_frozen_build_survives_unrelated_additions_to_the_shared_archive_cache(project: Path) -> None:
    """The archive lookup cache is shared between projects; only the lookups a build reads are pinned."""
    config = project / "config.yml"
    fetch(config)
    cache_file = project / "sources" / "ensembl_archive" / "archive_id_cache.json"
    lookups = json.loads(cache_file.read_text(encoding="utf-8"))
    lookups["ENSG00000999999"] = {"query_id": "ENSG00000999999", "found": False, "checked_release": 115}
    cache_file.write_text(json.dumps(lookups), encoding="utf-8")
    assert build(config, frozen=True).namespaces


def test_frozen_build_fails_when_a_lookup_it_reads_changes(project: Path) -> None:
    config = project / "config.yml"
    fetch(config)
    cache_file = project / "sources" / "ensembl_archive" / "archive_id_cache.json"
    lookups = json.loads(cache_file.read_text(encoding="utf-8"))
    lookups["ENSG00000000099"]["possible_replacement"] = ["ENSG00000000001"]
    cache_file.write_text(json.dumps(lookups), encoding="utf-8")
    with pytest.raises(ResourceError, match="re-run `alien fetch`"):
        build(config, frozen=True)
