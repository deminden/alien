from __future__ import annotations

import socket
from pathlib import Path
from typing import Any

import pandas as pd
import pytest


@pytest.fixture(autouse=True)
def _block_network(request: pytest.FixtureRequest, monkeypatch: pytest.MonkeyPatch) -> None:
    """Fail any test that tries to open a real network connection unless it is marked ``network``."""
    if request.node.get_closest_marker("network"):
        return

    original_connect = socket.socket.connect

    def guard(sock: socket.socket, address: object) -> None:
        if sock.family == getattr(socket, "AF_UNIX", None):
            original_connect(sock, address)  # process pools talk over local sockets
            return
        raise RuntimeError("Network access is disabled in the test suite; mark the test with @pytest.mark.network.")

    monkeypatch.setattr(socket.socket, "connect", guard)


MINI_GENES: tuple[tuple[str, str, str, str], ...] = (
    # symbol, stable Ensembl ID, previous symbols, alias symbols
    ("G1", "ENSG00000000101", "OLDG1|PREVDUP", "ALIASG1"),
    ("G2", "ENSG00000000102", "PREVDUP", "SHAREDALIAS"),
    ("G3", "ENSG00000000103", "", "SHAREDALIAS|WRONGALIAS"),
    ("G4", "ENSG00000000104", "", ""),
    ("G5", "ENSG00000000105", "", ""),
    ("G6", "ENSG00000000106", "", ""),
    ("G7", "ENSG00000000107", "", ""),
    ("G8", "ENSG00000000108", "", ""),
)


class MiniWorld:
    """A tiny offline gene world: HGNC table, GTF and helpers to write sources and configs."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self.source_dir = root / "sources"
        hgnc = self.source_dir / "hgnc"
        hgnc.mkdir(parents=True)
        rows = ["symbol\tprev_symbol\talias_symbol\tensembl_gene_id"]
        rows += [f"{s}\t{p}\t{a}\t{e}" for s, e, p, a in MINI_GENES]
        (hgnc / "hgnc_complete_set.txt").write_text("\n".join(rows) + "\n", encoding="utf-8")
        self.gtf = root / "genes.gtf"
        self.gtf.write_text(
            "".join(
                f'chr1\tT\tgene\t{i * 10 + 1}\t{i * 10 + 5}\t.\t+\t.\tgene_id "{e}.1"; gene_name "{s}"; '
                'gene_type "protein_coding";\n'
                for i, (s, e, _, _) in enumerate(MINI_GENES)
            ),
            encoding="utf-8",
        )

    def ensembl(self, symbol: str) -> str:
        return next(e for s, e, _, _ in MINI_GENES if s == symbol)

    def gmt(self, name: str, terms: dict[str, list[str]]) -> Path:
        path = self.root / name
        path.write_text("".join(f"{t}\tdesc\t" + "\t".join(g) + "\n" for t, g in terms.items()), encoding="utf-8")
        return path

    def target(self, name: str, output_genes: list[str] | None = None) -> dict[str, Any]:
        target: dict[str, Any] = {"name": name, "annotation": {"source": "Mini", "version": "1", "path": str(self.gtf)}}
        if output_genes is not None:
            target["output_genes"] = {"ids": [self.ensembl(symbol) for symbol in output_genes]}
        return target

    def config(self, gmt: Path, targets: list[dict[str, Any]] | None = None, **sections: Any) -> dict[str, Any]:
        cfg: dict[str, Any] = {
            "project": {"source_dir": str(self.source_dir), "outdir": str(self.root / "out")},
            "outputs": {"include_symbols": True},
            "sources": [
                {"type": "symbol_gmt", "path": str(gmt), "source_tag": "LOCAL", "family": "biology_process_pathway"}
            ],
            "targets": targets if targets is not None else [self.target("mini")],
            "filtering": {"min_size_default": 2, "max_size_default": 50},
            "ncbi_gene": {"enabled": False},
            "ensembl_archive": {"enabled": False},
        }
        for key, value in sections.items():
            if isinstance(value, dict) and isinstance(cfg.get(key), dict):
                cfg[key] = {**cfg[key], **value}
            else:
                cfg[key] = value
        return cfg


@pytest.fixture
def mini_world(tmp_path: Path) -> MiniWorld:
    return MiniWorld(tmp_path)


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def read_gmt(path: Path) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        result[parts[0]] = parts[2:]
    return result
