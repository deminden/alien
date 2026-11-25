"""Read output_genes / gene_filter Ensembl ID lists from inline IDs, matrices, tables or plain lists."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Final

import pandas as pd

from ..config.schema import GeneIdSpec
from ..errors import SourceError
from ..io.gtf import open_text
from ..model import GeneIdList, ensembl_stable_id, is_ensembl_gene_id
from ..schema import GENE_ID_LIST_COLUMNS

LOGGER = logging.getLogger("alien")
GENE_COLUMNS: Final = ("Ensembl_gene_ID", "ensembl_gene_id", "gene_id", "gene", "Gene", "id")
INDEX_COLUMNS: Final = frozenset({"index", "__index__"})


def empty_gene_id_list() -> GeneIdList:
    return GeneIdList(ids=None, table=pd.DataFrame(columns=list(GENE_ID_LIST_COLUMNS)), id_type="none")


def read_gene_id_list(spec: GeneIdSpec | None, *, with_symbols: bool = True) -> GeneIdList:
    if spec is None:
        return empty_gene_id_list()
    symbol_column = spec.symbol_column if with_symbols else None
    if spec.ids is not None:
        genes, symbols = [str(gene) for gene in spec.ids], [""] * len(spec.ids)
    else:
        assert spec.path is not None  # guaranteed by GeneIdSpec validation
        genes, symbols = _read_values(spec.path, spec.id_column, symbol_column)
    ids, table, id_type = _rows(genes, symbols)
    return GeneIdList(ids, table, id_type, spec.path, spec.id_column, symbol_column)


def _read_values(path: Path, id_column: str | None, symbol_column: str | None) -> tuple[list[str], list[str]]:
    if not path.exists():
        raise SourceError(f"Gene ID file {path} does not exist.")
    if path.suffix == ".parquet":
        wanted = [column for column in (id_column, symbol_column) if column and column not in INDEX_COLUMNS]
        try:
            frame = pd.read_parquet(path, columns=list(dict.fromkeys(wanted)))
        except (KeyError, ValueError) as exc:
            raise SourceError(f"Output gene ID column or symbol column was not found in {path}.") from exc
        genes = _id_values(frame, id_column, path, allow_index=True)
        return genes, _symbol_values(frame, symbol_column, path, len(genes))
    if id_column or symbol_column or _looks_tabular(path):
        frame = _read_table(path)
        genes = _id_values(frame, id_column, path, allow_index=False)
        return genes, _symbol_values(frame, symbol_column, path, len(genes))
    with open_text(path) as handle:
        genes = [line.strip().split()[0] for line in handle if line.strip()]
    return genes, [""] * len(genes)


def _row(gene: str, symbol: str) -> dict[str, str]:
    stripped = ensembl_stable_id(gene)
    if not is_ensembl_gene_id(gene):
        return {"input_gene_id": gene, "ensembl_gene_id": "", "gene_symbol": str(symbol).strip(), "id_type": "symbol"}
    id_type = "ensembl_versioned" if stripped != gene else "ensembl_stable"
    return {"input_gene_id": gene, "ensembl_gene_id": stripped, "gene_symbol": str(symbol).strip(), "id_type": id_type}


def _check_ids(n_rows: int, n_symbol_rows: int, stable: set[str]) -> None:
    if n_rows and not stable:
        raise SourceError(
            "output_genes and gene_filter must provide Ensembl IDs in id_column/ids; "
            "symbol_column is optional metadata."
        )
    if n_symbol_rows:
        LOGGER.warning(
            "%d output gene row(s) do not hold Ensembl gene IDs; they are kept in the audit table only.", n_symbol_rows
        )


def _rows(genes: list[str], symbols: list[str]) -> tuple[frozenset[str] | None, pd.DataFrame, str]:
    rows = [_row(gene, symbol) for gene, symbol in zip(genes, symbols, strict=True)]
    stable = {row["ensembl_gene_id"] for row in rows if row["id_type"] != "symbol"}
    _check_ids(len(rows), len(rows) - sum(row["id_type"] != "symbol" for row in rows), stable)
    types = {row["id_type"] for row in rows}
    id_type = "none" if not types else next(iter(types)) if len(types) == 1 else "mixed"
    table = pd.DataFrame(rows, columns=list(GENE_ID_LIST_COLUMNS)).drop_duplicates()
    return (frozenset(stable) if stable else None), table, id_type


def _id_values(frame: pd.DataFrame, column: str | None, path: Path, *, allow_index: bool) -> list[str]:
    if column:
        if column in INDEX_COLUMNS and allow_index:
            return frame.index.astype(str).tolist()
        if column not in frame.columns:
            raise SourceError(f"Output gene ID column {column!r} was not found in {path}.")
        return frame[column].astype(str).tolist()
    lowered = {name.lower(): name for name in frame.columns}
    for candidate in GENE_COLUMNS:
        match = candidate if candidate in frame.columns else lowered.get(candidate.lower())
        if match is not None:
            return frame[match].astype(str).tolist()
    if allow_index:
        return frame.index.astype(str).tolist()
    return frame[frame.columns[0]].astype(str).tolist()


def _symbol_values(frame: pd.DataFrame, column: str | None, path: Path, n_rows: int) -> list[str]:
    if not column:
        return [""] * n_rows
    if column not in frame:
        raise SourceError(f"Output gene symbol column {column!r} was not found in {path}.")
    return frame[column].fillna("").astype(str).tolist()


def _read_table(path: Path) -> pd.DataFrame:
    name = path.name.lower()
    if name.endswith((".tsv", ".tsv.gz", ".tab", ".tab.gz")):
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if name.endswith((".csv", ".csv.gz")):
        return pd.read_csv(path, sep=",", dtype=str, keep_default_na=False)
    return pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)


def _looks_tabular(path: Path) -> bool:
    with open_text(path) as handle:
        first = handle.readline()
    return "\t" in first or "," in first
