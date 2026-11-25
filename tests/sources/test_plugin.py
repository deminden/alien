from __future__ import annotations

from importlib.metadata import EntryPoint
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from alien.cache import ResourceCache
from alien.config import load_config
from alien.errors import ConfigError, SourceError
from alien.sources import load_source_memberships


def plugin_loader(options: dict[str, Any], cache: ResourceCache) -> pd.DataFrame:
    return pd.DataFrame({"term_id": [options["term"]] * 2, "gene_symbol": ["A", "B"]})


def bad_loader(options: dict[str, Any], cache: ResourceCache) -> Any:
    return {"not": "a frame"}


def _install(monkeypatch: pytest.MonkeyPatch, name: str, target: str) -> None:
    point = EntryPoint(name=name, value=f"{__name__}:{target}", group="alien.sources")
    monkeypatch.setattr("alien.sources.entry_points", lambda group: [point] if group == "alien.sources" else [])


def _config(name: str) -> Any:
    return load_config({"sources": [{"type": "plugin", "plugin": name, "options": {"term": "P1"}}], "targets": []})


def test_plugin_source_through_entry_point(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    _install(monkeypatch, "internal_db", "plugin_loader")
    frame = load_source_memberships(_config("internal_db"), ResourceCache(tmp_path))
    assert frame["term_id"].tolist() == ["P1", "P1"]
    assert "metadata_json" in frame.columns


def test_plugin_errors(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    _install(monkeypatch, "bad", "bad_loader")
    with pytest.raises(SourceError, match="did not return a DataFrame"):
        load_source_memberships(_config("bad"), ResourceCache(tmp_path))
    with pytest.raises(ConfigError, match="No source plugin named 'missing'"):
        load_source_memberships(_config("missing"), ResourceCache(tmp_path))
