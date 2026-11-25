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


# ----------------------------------------------------------------------------------------------- fake HTTP


class FakeResponse:
    def __init__(self, url: str, status: int, body: bytes, headers: dict[str, str] | None = None) -> None:
        self.url, self.status_code, self.content = url, status, body
        self.headers = headers or {"ETag": '"fixture"'}

    def __enter__(self) -> FakeResponse:
        return self

    def __exit__(self, *args: object) -> None:
        return None

    def raise_for_status(self) -> None:
        import requests

        if self.status_code >= 400:
            raise requests.HTTPError(f"{self.status_code} for {self.url}")

    def iter_content(self, chunk_size: int = 1024) -> Any:
        for start in range(0, len(self.content), chunk_size):
            yield self.content[start : start + chunk_size]

    def json(self) -> Any:
        import json

        return json.loads(self.content.decode("utf-8"))

    @property
    def text(self) -> str:
        return self.content.decode("utf-8")


class FakeServer:
    """Routes URL -> bytes, (status, bytes), an exception, or a callable(json_payload) for POST."""

    def __init__(self) -> None:
        self.routes: dict[str, Any] = {}
        self.calls: list[tuple[str, str, Any]] = []

    def respond(self, method: str, url: str, payload: Any = None) -> FakeResponse:
        self.calls.append((method, url, payload))
        route = self.routes.get(url)
        if route is None:
            return FakeResponse(url, 404, b"not found")
        if isinstance(route, Exception):
            raise route
        if callable(route):
            route = route(payload)
        if isinstance(route, tuple):
            return FakeResponse(url, route[0], route[1])
        if isinstance(route, list | dict):
            import json

            route = json.dumps(route)
        if isinstance(route, str):
            route = route.encode("utf-8")
        return FakeResponse(url, 200, route)

    def urls(self, method: str = "GET") -> list[str]:
        return [url for m, url, _ in self.calls if m == method]


class FakeSession:
    def __init__(self, server: FakeServer) -> None:
        self.server = server
        self.headers: dict[str, str] = {}

    def __enter__(self) -> FakeSession:
        return self

    def __exit__(self, *args: object) -> None:
        return None

    def get(self, url: str, **kwargs: Any) -> FakeResponse:
        return self.server.respond("GET", url)

    def post(self, url: str, json: Any = None, **kwargs: Any) -> FakeResponse:
        return self.server.respond("POST", url, json)


@pytest.fixture
def fake_http(monkeypatch: pytest.MonkeyPatch) -> FakeServer:
    from alien.io import http

    server = FakeServer()
    monkeypatch.setattr(http, "SESSION_FACTORY", lambda **kwargs: FakeSession(server))
    return server


# ----------------------------------------------------------------------------------------------- targets


def annotation_frame(rows: list[dict[str, Any]]) -> pd.DataFrame:
    from alien.schema import ANNOTATION_COLUMNS

    defaults = {
        "gencode_version": "v1",
        "ensembl_gene_id_versioned": "",
        "gene_biotype": "protein_coding",
        "seqname": "chr1",
        "start": "1",
        "end": "2",
        "strand": "+",
        "source_gtf": "fixture.gtf",
        "target_id_metadata_source": "annotation_gtf",
        "is_in_output_genes": True,
    }
    return pd.DataFrame([{**defaults, **row} for row in rows], columns=list(ANNOTATION_COLUMNS))


def make_target(name: str, rows: list[dict[str, Any]], output_ids: set[str] | None = None) -> Any:
    from alien.model import GeneIdList, Target
    from alien.resources.output_genes import empty_gene_id_list

    annotation = annotation_frame(rows)
    ids = frozenset(annotation["ensembl_gene_id"])
    output = empty_gene_id_list()
    if output_ids is not None:
        table = pd.DataFrame(
            {
                "input_gene_id": sorted(output_ids),
                "ensembl_gene_id": sorted(output_ids),
                "gene_symbol": "",
                "id_type": "ensembl_stable",
            }
        )
        output = GeneIdList(frozenset(output_ids), table, "ensembl_stable")
    return Target(name, annotation, "Fixture v1", Path("fixture.gtf"), output, empty_gene_id_list(), ids, ids)


def membership(
    term_id: str, gene: str, ensembl: str = "", source_tag: str = "REACTOME", family: str = "biology_process_pathway"
) -> dict[str, str]:
    return {
        "term_id": term_id,
        "original_name": term_id,
        "display_name": term_id,
        "description": "",
        "source": "test",
        "source_tag": source_tag,
        "collection": "C2",
        "subcollection": "",
        "family": family,
        "aspect": "pathway",
        "db_version": "",
        "source_url": "",
        "source_license_note": "",
        "gene_symbol": gene,
        "gene_ensembl_from_source": ensembl,
        "metadata_json": "{}",
    }


def hgnc_maps(rows: list[tuple[str, str, str, str]]) -> Any:
    """rows: (symbol, previous symbols, alias symbols, Ensembl ID)."""
    from alien.resources.hgnc import build_hgnc_maps

    return build_hgnc_maps(
        pd.DataFrame([{"symbol": s, "prev_symbol": p, "alias_symbol": a, "ensembl_gene_id": e} for s, p, a, e in rows])
    )
