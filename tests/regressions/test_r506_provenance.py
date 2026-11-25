"""5.6: provenance reflects what the build actually did."""

from __future__ import annotations

import json

from alien import build

from conftest import MiniWorld


def test_provenance_is_built_from_the_run(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"T1": ["G1", "G2", "G3"]})
    result = build(config=mini_world.config(gmt))
    text = (result.outdir / "metadata" / "source_provenance.json").read_text(encoding="utf-8")
    provenance = json.loads(text)
    assert "2026.1.Hs" not in text
    assert "official_v0_1_scope" not in text
    assert provenance["alien_version"]
    assert {"pandas", "pyarrow", "PyYAML", "requests", "rdata", "urllib3"} <= set(provenance["packages"])
    assert "NCBI_Gene" not in provenance["resources"]
    assert "MSigDB" not in provenance["resources"]
    inputs = {entry["path"]: entry for entry in provenance["input_files"]}
    assert str(gmt) in inputs
    assert len(inputs[str(gmt)]["sha256"]) == 64
    assert len(provenance["config_sha256"]) == 64
