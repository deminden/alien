from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from alien.cache import LockFile, ResourceCache
from alien.config.schema import EnsemblArchiveConfig
from alien.errors import ResourceError
from alien.resources.ensembl_archive import EnsemblArchiveResolver, archive_target_matches, unresolved_status

from conftest import FakeServer

ENDPOINT = "https://example.test/archive/id"
INFO = "https://example.test/info/data/?content-type=application/json"


def _resolver(tmp_path: Path, **cache: Any) -> EnsemblArchiveResolver:
    config = EnsemblArchiveConfig(endpoint=ENDPOINT, batch_size=50, retry_attempts=1, retry_backoff_seconds=0)
    return EnsemblArchiveResolver(config, ResourceCache(tmp_path, **cache))


def _ok(payload: Any) -> list[dict[str, Any]]:
    return [{"id": stable_id, "is_current": "1"} for stable_id in payload["id"]]


def test_retries_cached_errors_and_saves_atomically(tmp_path: Path, fake_http: FakeServer) -> None:
    fake_http.routes[ENDPOINT] = _ok
    resolver = _resolver(tmp_path)
    resolver.cache = {
        "ENSG_SUCCESS": {"query_id": "ENSG_SUCCESS", "found": True},
        "ENSG_ERROR": {"query_id": "ENSG_ERROR", "found": False, "error": "timeout"},
    }
    resolver.prefetch(["ENSG_SUCCESS", "ENSG_ERROR", "ENSG_NEW"])
    assert [payload["id"] for _, _, payload in fake_http.calls if payload] == [["ENSG_ERROR", "ENSG_NEW"]]
    on_disk = json.loads(resolver.cache_path.read_text(encoding="utf-8"))
    assert on_disk["ENSG_ERROR"]["found"] is True and on_disk["ENSG_NEW"]["found"] is True
    assert "ensembl_archive:lookups" in resolver.resources.used


def test_splits_failed_batches(tmp_path: Path, fake_http: FakeServer) -> None:
    def flaky(payload: Any) -> Any:
        return (504, b"timeout") if len(payload["id"]) > 1 else json.dumps(_ok(payload)).encode()

    fake_http.routes[ENDPOINT] = flaky
    resolver = _resolver(tmp_path)
    resolver.prefetch(["ENSG_A", "ENSG_B"])
    assert [payload["id"] for _, _, payload in fake_http.calls if payload] == [
        ["ENSG_A", "ENSG_B"],
        ["ENSG_A"],
        ["ENSG_B"],
    ]
    assert resolver.lookup("ENSG_A")["found"] and resolver.lookup("ENSG_B")["found"]


def test_single_id_failure_is_recorded_as_error(tmp_path: Path, fake_http: FakeServer) -> None:
    fake_http.routes[ENDPOINT] = lambda payload: (500, b"boom")
    resolver = _resolver(tmp_path)
    resolver.prefetch(["ENSG_X"])
    assert resolver.lookup("ENSG_X")["error"]


def test_not_found_answers_are_rechecked_after_a_new_release(tmp_path: Path, fake_http: FakeServer) -> None:
    fake_http.routes[INFO] = json.dumps({"releases": [113]})
    fake_http.routes[ENDPOINT] = lambda payload: []
    resolver = _resolver(tmp_path)
    resolver.prefetch(["ENSG_GONE"])
    assert resolver.lookup("ENSG_GONE") == {"query_id": "ENSG_GONE", "found": False, "checked_release": 113}
    fake_http.calls.clear()
    _resolver(tmp_path).prefetch(["ENSG_GONE"])
    assert [m for m, *_ in fake_http.calls] == ["GET"]  # same release: no new lookup
    fake_http.routes[INFO] = json.dumps({"releases": [114]})
    fake_http.calls.clear()
    _resolver(tmp_path).prefetch(["ENSG_GONE"])
    assert [m for m, *_ in fake_http.calls] == ["GET", "POST"]


def test_offline_and_frozen(tmp_path: Path) -> None:
    offline = _resolver(tmp_path, offline=True)
    offline.prefetch(["ENSG_Q"])
    assert offline.lookup("ENSG_Q")["error"] == "offline"
    assert not offline.cache_path.exists()
    frozen = _resolver(tmp_path, frozen=True, lock=LockFile(path=tmp_path / "alien.lock.json"))
    with pytest.raises(ResourceError, match="--frozen"):
        frozen.prefetch(["ENSG_Q"])


def test_target_matches_and_statuses() -> None:
    info = {
        "found": True,
        "archive_id": "ENSG1",
        "is_current": "1",
        "latest": "ENSG2.3",
        "possible_replacement": ["ENSG3"],
    }
    assert archive_target_matches(info, {"ENSG1", "ENSG2", "ENSG3", "ENSG4"}) == ["ENSG1", "ENSG2", "ENSG3"]
    assert archive_target_matches({"found": False}, {"E1"}) == []
    assert unresolved_status(info, ["E1", "E2"], "x") == "archive_ambiguous_target_replacement__x"
    assert unresolved_status({"error": "t"}, [], "x") == "archive_query_failed__x"
    assert unresolved_status({"found": False}, [], "x") == "archive_not_found__x"
    assert unresolved_status({"found": True}, [], "x") == "archive_no_target_replacement__x"
    assert unresolved_status({}, [], "x") == "archive_not_checked__x"
