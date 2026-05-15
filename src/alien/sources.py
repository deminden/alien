from __future__ import annotations

import hashlib
import json
import re
import shutil
import zipfile
from functools import lru_cache
from pathlib import Path
from typing import Any
from urllib.parse import quote

import pandas as pd
import requests

from .utils import (
    LOGGER,
    concat_frames,
    download_file,
    empty_term_gene_frame,
    ensure_dirs,
    first_existing_column,
    open_text,
    parse_gmt,
    sanitize_id,
    strip_ensembl_version,
    write_json,
    write_tsv_gz,
)


MSIGDB_SOURCES = {
    "REACTOME": ("C2", "CP:REACTOME", "biology_process_pathway", "pathway"),
    "WIKIPATHWAYS": ("C2", "CP:WIKIPATHWAYS", "biology_process_pathway", "pathway"),
    "KEGG_MEDICUS": ("C2", "CP:KEGG_MEDICUS", "biology_process_pathway", "pathway"),
    "GOBP": ("C5", "GO:BP", "biology_process_pathway", "biological_process"),
    "GOMF": ("C5", "GO:MF", "biology_function_location", "molecular_function"),
    "GOCC": ("C5", "GO:CC", "biology_function_location", "cellular_component"),
    "HPO": ("C5", "HPO", "disease_phenotype", "phenotype"),
    "C6": ("C6", "", "cancer_dependency_state", "oncogenic_signature"),
    "C9": ("C9", "", "cancer_dependency_state", "dependency_signature"),
    "C4_3CA": ("C4", "3CA", "cancer_dependency_state", "cancer_signature"),
    "C4_CGN": ("C4", "CGN", "cancer_dependency_state", "cancer_signature"),
    "C4_CM": ("C4", "CM", "cancer_dependency_state", "cancer_signature"),
}

ENRICHR_BASE_URL = "https://maayanlab.cloud/Enrichr"
ENRICHR_DATASET_STATISTICS_ENDPOINT = f"{ENRICHR_BASE_URL}/datasetStatistics"
ENRICHR_LIBRARY_DOWNLOAD_ENDPOINT = f"{ENRICHR_BASE_URL}/geneSetLibrary?mode=text&libraryName={{library}}"

HGNC_URL = "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt"
HGNC_FALLBACK_URLS = [
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
]

GENCODE_HUMAN_BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"

NCBI_GENE_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
NCBI_GENE_HISTORY_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_history.gz"
MSIGDB_RELEASES = {
    "2026.1": {
        "zip_url": "https://zenodo.org/records/18968178/files/msigdb.2026.1.zip?download=1",
        "zip_md5": "512ba99c6827141a9d471972b812d4ac",
        "zip_name": "msigdb.2026.1.zip",
        "summary_rds": "msigdb.2026.1.summary.rds",
    }
}

GENE_COLUMNS = ["Ensembl_gene_ID", "ensembl_gene_id", "gene_id", "gene", "Gene", "id"]
CANONICAL_COLUMNS = list(empty_term_gene_frame().columns)
DEFAULT_GTF_ATTRIBUTES = {
    "gene_id": ["gene_id"],
    "gene_name": ["gene_name"],
    "gene_biotype": ["gene_type", "gene_biotype"],
}


def load_source_memberships(cfg: dict[str, Any], force_download: bool = False) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    source_dir = Path(cfg.get("project", {}).get("source_dir", "data/alien_sources"))
    for spec in cfg.get("sources", []):
        if not spec.get("enabled", True):
            continue
        frames.append(_load_source_spec(spec, source_dir, force_download=force_download))
    return concat_frames(frames)


def _load_source_spec(
    spec: dict[str, Any],
    source_dir: Path,
    force_download: bool = False,
) -> pd.DataFrame:
    _reject_optional_source_key(spec)
    source_type = str(spec.get("type", "")).strip()
    if source_type == "canonical_tsv":
        return concat_frames([read_canonical_memberships(path, spec) for path in _source_paths(spec, source_dir)])
    if source_type == "msigdb_cache":
        path = Path(spec.get("path") or source_dir / "msigdb")
        return read_msigdb_cache(path, bool(spec.get("include_c4_cm", False)), spec)
    if source_type == "msigdb_remote":
        return read_msigdb_remote(spec, source_dir, force=force_download)
    if source_type == "enrichr_remote":
        return read_enrichr_remote(spec, source_dir, force=force_download)
    if source_type == "msigdb_tsv":
        return concat_frames([read_msigdb_like_table(path, spec) for path in _source_paths(spec, source_dir)])
    if source_type in {"symbol_gmt", "enrichr_gmt"}:
        return concat_frames([read_symbol_gmt_source(path, spec) for path in _source_paths(spec, source_dir)])
    raise ValueError(f"Unsupported source type: {source_type or '<empty>'}")


def _reject_optional_source_key(spec: dict[str, Any]) -> None:
    if "optional" in spec:
        raise ValueError("Source option 'optional' is not supported; configured sources are required. Use enabled: false to exclude a source.")


def _source_paths(spec: dict[str, Any], source_dir: Path) -> list[Path]:
    values = spec.get("paths")
    if values is None:
        values = [spec.get("path")]
    paths = [Path(value) for value in values if value]
    resolved = [path if path.is_absolute() else Path(path) for path in paths]
    if not resolved and spec.get("type") in {"symbol_gmt", "enrichr_gmt"}:
        default_dir = source_dir / str(spec.get("name", "symbol_gmt"))
        resolved = sorted(default_dir.glob("*.gmt"))
    return resolved


def read_canonical_memberships(path: Path, spec: dict[str, Any] | None = None) -> pd.DataFrame:
    spec = spec or {}
    raw = pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)
    missing = [column for column in ["term_id", "gene_symbol"] if column not in raw]
    if missing:
        raise ValueError(f"{path} is missing canonical columns: {', '.join(missing)}")
    frame = pd.DataFrame({column: raw[column] if column in raw else "" for column in CANONICAL_COLUMNS})
    frame["term_id"] = frame["term_id"].map(sanitize_id)
    frame["source"] = frame["source"].replace("", spec.get("source", spec.get("name", "canonical")))
    frame["source_tag"] = frame["source_tag"].replace("", spec.get("source_tag", sanitize_id(frame["source"].iloc[0]).upper()))
    frame["family"] = frame["family"].replace("", spec.get("family", "biology_process_pathway"))
    frame["aspect"] = frame["aspect"].replace("", spec.get("aspect", "library"))
    frame["collection"] = frame["collection"].replace("", spec.get("collection", Path(path).stem))
    frame["subcollection"] = frame["subcollection"].replace("", spec.get("subcollection", ""))
    frame["original_name"] = frame["original_name"].where(frame["original_name"].ne(""), frame["term_id"])
    frame["display_name"] = frame["display_name"].where(frame["display_name"].ne(""), frame["original_name"])
    frame["description"] = frame["description"].where(frame["description"].ne(""), frame["display_name"])
    if "gene_id" in raw and "gene_id_namespace" in raw:
        ensembl_mask = raw["gene_id_namespace"].str.lower().str.contains("ensembl|ensg", regex=True, na=False)
        frame.loc[ensembl_mask, "gene_ensembl_from_source"] = raw.loc[ensembl_mask, "gene_id"]
    elif "gene_id" in raw:
        frame["gene_ensembl_from_source"] = raw["gene_id"].map(lambda value: value if strip_ensembl_version(value).startswith("ENSG") else "")
    frame["metadata_json"] = frame["metadata_json"].replace("", "{}")
    return frame[CANONICAL_COLUMNS]


def read_msigdb_cache(
    path: Path,
    include_c4_cm: bool,
    spec: dict[str, Any] | None = None,
) -> pd.DataFrame:
    spec = spec or {}
    _reject_optional_source_key(spec)
    frames: list[pd.DataFrame] = []
    for tag, (collection, subcollection, family, aspect) in MSIGDB_SOURCES.items():
        if tag == "C4_CM" and not include_c4_cm:
            continue
        file_path = path / f"msigdbr_{tag}.tsv.gz"
        if not file_path.exists():
            raise FileNotFoundError(f"Missing MSigDB cache file: {file_path}")
        frames.append(read_msigdb_like_table(file_path, {"source_tag": tag, "collection": collection, "subcollection": subcollection, "family": family, "aspect": aspect}))
    return concat_frames(frames)


def read_msigdb_remote(spec: dict[str, Any], source_dir: Path, force: bool = False) -> pd.DataFrame:
    cache_dir = Path(spec.get("cache_dir") or source_dir / "msigdb_remote")
    force_refresh = force or bool(spec.get("force", False))
    release = _msigdb_release_info(spec)
    release_dir = ensure_msigdb_remote_cache(
        cache_dir=cache_dir,
        release=release,
        force=force_refresh,
        timeout=int(spec.get("timeout_seconds", 600)),
    )
    alien_cache_path = _msigdb_remote_alien_cache_path(release_dir, release, spec)
    if alien_cache_path.exists() and not force_refresh:
        LOGGER.info("Loading cached normalized MSigDB memberships from %s", alien_cache_path)
        return _read_canonical_parquet(alien_cache_path)

    summary = _read_msigdb_summary(release_dir, release, force_refresh)
    db_species = str(spec.get("db_species", "HS")).upper()
    collection = spec.get("collection")
    subcollection = spec.get("subcollection")
    summary = summary[summary["db_target_species"].astype(str).str.upper().eq(db_species)]
    if collection:
        summary = summary[summary["gs_collection"].astype(str).eq(str(collection))]
    if subcollection:
        sub = str(subcollection)
        summary = summary[
            summary["gs_subcollection"].astype(str).eq(sub)
            | summary["gs_subcollection"].astype(str).str.replace(r".*:", "", regex=True).eq(sub)
        ]
    if summary.empty:
        raise ValueError("No MSigDB remote collection files matched the requested filters.")

    frames = []
    for rds_name in sorted(set(summary["df_rds"].astype(str))):
        frame = _read_rds_dataframe(release_dir / rds_name)
        frames.append(_normalize_msigdb_remote_frame(frame, spec))
    memberships = concat_frames(frames)
    _write_canonical_parquet(memberships, alien_cache_path)
    return memberships


def _read_msigdb_summary(release_dir: Path, release: dict[str, str], force: bool = False) -> pd.DataFrame:
    summary_rds = release_dir / str(release["summary_rds"])
    summary_cache = release_dir / "alien_cache" / f"{summary_rds.stem}.parquet"
    if summary_cache.exists() and not force:
        return pd.read_parquet(summary_cache).fillna("").astype(str)
    summary = _read_rds_dataframe(summary_rds)
    ensure_dirs(summary_cache.parent)
    summary.to_parquet(summary_cache, index=False)
    return summary


def _msigdb_remote_alien_cache_path(release_dir: Path, release: dict[str, str], spec: dict[str, Any]) -> Path:
    relevant = {
        "zip_md5": release.get("zip_md5", ""),
        "summary_rds": release.get("summary_rds", ""),
        "db_species": str(spec.get("db_species", "HS")).upper(),
        "collection": spec.get("collection", ""),
        "subcollection": spec.get("subcollection", ""),
        "family": spec.get("family", ""),
        "aspect": spec.get("aspect", ""),
    }
    digest = hashlib.md5(json.dumps(relevant, sort_keys=True).encode("utf-8")).hexdigest()[:16]
    label_bits = [
        relevant["db_species"],
        str(relevant["collection"] or "all"),
        str(relevant["subcollection"] or "all"),
        str(relevant["family"] or "default"),
        str(relevant["aspect"] or "default"),
    ]
    label = sanitize_id("__".join(label_bits))
    return release_dir / "alien_cache" / f"memberships__{label}__{digest}.parquet"


def _read_canonical_parquet(path: Path) -> pd.DataFrame:
    frame = pd.read_parquet(path).fillna("").astype(str)
    for column in CANONICAL_COLUMNS:
        if column not in frame:
            frame[column] = ""
    return frame[CANONICAL_COLUMNS]


def _write_canonical_parquet(frame: pd.DataFrame, path: Path) -> None:
    ensure_dirs(path.parent)
    frame[CANONICAL_COLUMNS].to_parquet(path, index=False)


def ensure_msigdb_remote_cache(cache_dir: Path, release: dict[str, str], force: bool = False, timeout: int = 600) -> Path:
    ensure_dirs(cache_dir)
    zip_name = str(release["zip_name"])
    zip_path = cache_dir / zip_name
    release_dir = cache_dir / zip_path.stem
    summary_path = release_dir / str(release["summary_rds"])
    if summary_path.exists() and not force:
        return release_dir
    if force and release_dir.exists():
        shutil.rmtree(release_dir)
    if not zip_path.exists():
        _download_checked_file(str(release["zip_url"]), zip_path, str(release["zip_md5"]), timeout)
    else:
        _verify_md5(zip_path, str(release["zip_md5"]))
    ensure_dirs(release_dir)
    with zipfile.ZipFile(zip_path) as archive:
        archive.extractall(release_dir)
    if not summary_path.exists():
        nested = list(release_dir.rglob(str(release["summary_rds"])))
        if nested:
            summary_path = nested[0]
        else:
            raise FileNotFoundError(f"Expected MSigDB summary RDS was not found after extraction: {release['summary_rds']}")
    return summary_path.parent


def _msigdb_release_info(spec: dict[str, Any]) -> dict[str, str]:
    version = str(spec.get("version", "2026.1"))
    release = dict(MSIGDB_RELEASES.get(version, {}))
    release.update({key: str(value) for key, value in spec.get("release", {}).items() if value is not None})
    if not release and {"zip_url", "zip_md5"}.issubset(spec):
        release = {key: str(spec[key]) for key in ["zip_url", "zip_md5"]}
    if "zip_url" not in release or "zip_md5" not in release:
        raise ValueError(f"Unknown MSigDB release {version!r}; provide release.zip_url and release.zip_md5.")
    zip_name = str(release.get("zip_name") or re.sub(r"\?.*$", "", Path(str(release["zip_url"])).name))
    release["zip_name"] = zip_name
    release.setdefault("summary_rds", zip_name.replace(".zip", ".summary.rds"))
    return release


def _download_checked_file(url: str, path: Path, md5: str, timeout: int) -> None:
    ensure_dirs(path.parent)
    tmp = path.with_suffix(path.suffix + ".tmp")
    LOGGER.info("Downloading MSigDB release archive from %s", url)
    with requests.get(url, stream=True, timeout=(15, timeout)) as response:
        response.raise_for_status()
        digest = hashlib.md5()
        with tmp.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    digest.update(chunk)
                    handle.write(chunk)
    if digest.hexdigest() != md5:
        tmp.unlink(missing_ok=True)
        raise RuntimeError("Downloaded MSigDB archive does not match the expected MD5 checksum.")
    shutil.move(tmp, path)


def _verify_md5(path: Path, expected: str) -> None:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    if digest.hexdigest() != expected:
        raise RuntimeError(f"Cached MSigDB archive failed MD5 verification: {path}")


@lru_cache(maxsize=16)
def _read_rds_dataframe(path: Path) -> pd.DataFrame:
    try:
        import rdata
    except ImportError as exc:
        raise RuntimeError("Reading MSigDB remote releases requires the Python package 'rdata'.") from exc
    obj = rdata.read_rds(path)
    if isinstance(obj, pd.DataFrame):
        return obj.fillna("").astype(str)
    raise TypeError(f"Expected RDS data frame in {path}, got {type(obj).__name__}.")


def _normalize_msigdb_remote_frame(raw: pd.DataFrame, spec: dict[str, Any]) -> pd.DataFrame:
    if "db_gene_symbol" not in raw or "gs_name" not in raw:
        raise ValueError("MSigDB remote RDS must contain db_gene_symbol and gs_name columns.")
    raw = _filter_msigdb_remote_rows(raw, spec).copy()
    if raw.empty:
        raise ValueError("No MSigDB remote rows matched the requested filters.")
    raw["gene_symbol"] = raw["db_gene_symbol"]
    raw["ensembl_gene"] = raw.get("db_ensembl_gene", "")
    subcollection = raw.get("gs_subcollection", pd.Series([""] * len(raw), index=raw.index)).fillna("").astype(str).str.strip()
    collection = raw.get("gs_collection", pd.Series([""] * len(raw), index=raw.index)).fillna("").astype(str).str.strip()
    raw["_alien_source_tag"] = subcollection.where(subcollection.ne(""), collection).replace("", "MSIGDB").map(_msigdb_remote_source_tag)
    frames = []
    for tag, tag_frame in raw.groupby("_alien_source_tag", sort=True):
        family = str(spec.get("family") or _family_for_msigdb_tag(tag))
        aspect = str(spec.get("aspect") or _aspect_for_msigdb_tag(tag))
        collection = str(tag_frame["gs_collection"].iloc[0]) if "gs_collection" in tag_frame else ""
        subcollection = str(tag_frame["gs_subcollection"].iloc[0]) if "gs_subcollection" in tag_frame else ""
        frames.append(_normalize_msigdb(tag_frame, tag, collection, subcollection, family, aspect))
    return concat_frames(frames)


def _filter_msigdb_remote_rows(raw: pd.DataFrame, spec: dict[str, Any]) -> pd.DataFrame:
    frame = raw
    db_species = str(spec.get("db_species", "HS")).upper()
    if "db_target_species" in frame:
        frame = frame[frame["db_target_species"].astype(str).str.upper().eq(db_species)]
    collection = spec.get("collection")
    if collection and "gs_collection" in frame:
        frame = frame[frame["gs_collection"].astype(str).eq(str(collection))]
    subcollection = spec.get("subcollection")
    if subcollection and "gs_subcollection" in frame:
        sub = str(subcollection)
        frame = frame[
            frame["gs_subcollection"].astype(str).eq(sub)
            | frame["gs_subcollection"].astype(str).str.replace(r".*:", "", regex=True).eq(sub)
        ]
    return frame


def _msigdb_remote_source_tag(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return "MSIGDB"
    tail = text.split(":")[-1]
    if tail in {"BP", "MF", "CC"}:
        return "GO" + tail
    return sanitize_id(tail).upper()


def _msigdb_remote_row_source_tag(row: pd.Series) -> str:
    subcollection = str(row.get("gs_subcollection", "") or "").strip()
    collection = str(row.get("gs_collection", "") or "").strip()
    return _msigdb_remote_source_tag(subcollection or collection or "MSIGDB")


def _family_for_msigdb_tag(tag: str) -> str:
    return MSIGDB_SOURCES.get(tag, ("", "", "biology_process_pathway", ""))[2]


def _aspect_for_msigdb_tag(tag: str) -> str:
    return MSIGDB_SOURCES.get(tag, ("", "", "", "library"))[3]


def read_msigdb_like_table(path: Path, spec: dict[str, Any] | None = None) -> pd.DataFrame:
    spec = spec or {}
    raw = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    name_col = str(spec.get("term_name_column", "gs_name"))
    symbol_col = str(spec.get("symbol_column", "gene_symbol"))
    if name_col not in raw or symbol_col not in raw:
        raise ValueError(f"{path} must contain {name_col} and {symbol_col} columns.")
    tag = str(spec.get("source_tag") or spec.get("name") or "MSIGDB")
    collection = str(spec.get("collection", ""))
    subcollection = str(spec.get("subcollection", ""))
    family = str(spec.get("family", "biology_process_pathway"))
    aspect = str(spec.get("aspect", "library"))
    renamed = raw.rename(columns={name_col: "gs_name", symbol_col: "gene_symbol"}).copy()
    if "ensembl_gene" not in renamed:
        gene_id_col = spec.get("gene_id_column")
        renamed["ensembl_gene"] = renamed.get(gene_id_col, "") if gene_id_col else ""
    return _normalize_msigdb(renamed, tag, collection, subcollection, family, aspect)


def read_symbol_gmt_source(path: Path, spec: dict[str, Any] | None = None) -> pd.DataFrame:
    spec = spec or {}
    rows: list[dict[str, str]] = []
    source_name = str(spec.get("source", spec.get("name", "symbol_gmt")))
    source_tag = str(spec.get("source_tag", sanitize_id(source_name).upper()))
    family = str(spec.get("family", "biology_process_pathway"))
    aspect = str(spec.get("aspect", "library"))
    collection = str(spec.get("collection", path.stem))
    prefix = sanitize_id(str(spec.get("term_prefix", source_tag)))
    for record in parse_gmt(path):
        original_name = str(record["term_id"])
        term_id = f"{prefix}__{sanitize_id(original_name)}"
        for gene in record["genes"]:
            rows.append(
                {
                    "term_id": term_id,
                    "original_name": original_name,
                    "display_name": original_name,
                    "description": str(record["description"]),
                    "source": source_name,
                    "source_tag": source_tag,
                    "collection": collection,
                    "subcollection": str(spec.get("subcollection", "")),
                    "family": family,
                    "aspect": aspect,
                    "db_version": str(spec.get("db_version", "")),
                    "source_url": str(spec.get("source_url", "")),
                    "source_license_note": str(spec.get("source_license_note", "")),
                    "gene_symbol": gene,
                    "gene_ensembl_from_source": "",
                    "metadata_json": "{}",
                }
            )
    return pd.DataFrame(rows, columns=CANONICAL_COLUMNS)


def read_enrichr_remote(
    spec: dict[str, Any],
    source_dir: Path,
    force: bool = False,
) -> pd.DataFrame:
    spec = spec or {}
    _reject_optional_source_key(spec)
    cache_dir = Path(spec.get("cache_dir") or source_dir / "enrichr")
    force_refresh = force or bool(spec.get("force", False))
    library_specs = _enrichr_library_specs(spec)
    metadata = ensure_enrichr_metadata(cache_dir, spec, force=force_refresh)
    names = _extract_library_names(metadata)
    if not names:
        raise ValueError("Enrichr metadata did not contain any library names.")

    frames: list[pd.DataFrame] = []
    for library_spec in library_specs:
        selected, candidates, match_method = resolve_enrichr_library(names, library_spec)
        gmt_path = ensure_enrichr_library(cache_dir, selected, spec, force=force_refresh)
        frames.append(
            _normalize_enrichr_gmt(
                gmt_path,
                selected,
                source_tag=str(library_spec.get("source_tag") or spec.get("source_tag") or sanitize_id(selected).upper()),
                family=str(library_spec.get("family") or spec.get("family") or "biology_process_pathway"),
                aspect=str(library_spec.get("aspect") or spec.get("aspect") or "library"),
                metadata={
                    "selected_library": selected,
                    "configured_name": str(library_spec.get("name", "")),
                    "match": str(library_spec.get("match", "")),
                    "match_method": match_method,
                    "candidate_libraries": candidates,
                },
            )
        )
    return concat_frames(frames)


def ensure_enrichr_metadata(cache_dir: Path, spec: dict[str, Any], force: bool = False) -> object:
    ensure_dirs(cache_dir)
    path = cache_dir / "datasetStatistics.json"
    if force or not path.exists():
        endpoint = str(spec.get("dataset_statistics_endpoint", ENRICHR_DATASET_STATISTICS_ENDPOINT))
        timeout = int(spec.get("timeout_seconds", 120))
        LOGGER.info("Downloading Enrichr library metadata from %s", endpoint)
        response = requests.get(endpoint, timeout=(15, timeout))
        response.raise_for_status()
        path.write_text(response.text, encoding="utf-8")
    return json.loads(path.read_text(encoding="utf-8"))


def ensure_enrichr_library(cache_dir: Path, library: str, spec: dict[str, Any], force: bool = False) -> Path:
    ensure_dirs(cache_dir)
    path = cache_dir / f"{sanitize_id(library)}.gmt"
    if path.exists() and not force:
        return path
    endpoint = str(spec.get("library_download_endpoint", ENRICHR_LIBRARY_DOWNLOAD_ENDPOINT))
    timeout = int(spec.get("timeout_seconds", 120))
    url = endpoint.format(library=quote(library))
    LOGGER.info("Downloading Enrichr library %s", library)
    response = requests.get(url, timeout=(15, timeout))
    response.raise_for_status()
    if not response.text.strip():
        raise ValueError(f"Downloaded Enrichr library {library!r} was empty.")
    path.write_text(response.text, encoding="utf-8")
    return path


def resolve_enrichr_library(names: list[str], library_spec: dict[str, Any]) -> tuple[str, list[str], str]:
    name = str(library_spec.get("name", "")).strip()
    if name and name in names:
        return name, [name], "exact"

    match = str(library_spec.get("match", "")).strip()
    if match:
        pattern = re.compile(match, re.I)
        candidates = sorted({library for library in names if pattern.search(library)})
        if candidates:
            selected = sorted(candidates, key=lambda value: (-_latest_year(value), len(value), value))[0]
            return selected, candidates, "regex"
        label = name or match
        raise ValueError(f"No Enrichr libraries matched {label!r}.")

    if name:
        raise ValueError(f"Enrichr library {name!r} was not found in datasetStatistics metadata.")
    raise ValueError("Every Enrichr remote library must define name or match.")


def _enrichr_library_specs(spec: dict[str, Any]) -> list[dict[str, Any]]:
    libraries = spec.get("libraries")
    if libraries is None:
        if spec.get("library"):
            libraries = [{"name": spec["library"]}]
        elif spec.get("name") or spec.get("match"):
            libraries = [spec]
        else:
            raise ValueError("enrichr_remote source must define libraries.")
    result: list[dict[str, Any]] = []
    for library in libraries:
        if isinstance(library, str):
            result.append({"name": library})
        elif isinstance(library, dict):
            _reject_optional_source_key(library)
            result.append(dict(library))
        else:
            raise TypeError("Each enrichr_remote library must be a string or mapping.")
    return result


def load_hgnc(source_dir: Path, force: bool = False) -> dict[str, object]:
    # Download HGNC symbol history
    hgnc_dir = source_dir / "hgnc"
    path = hgnc_dir / "hgnc_complete_set.txt"
    last_error: Exception | None = None
    for url in [HGNC_URL] + HGNC_FALLBACK_URLS:
        try:
            download_file(url, path, force=force)
            last_error = None
            break
        except Exception as exc:
            last_error = exc
            LOGGER.warning("HGNC download failed from %s: %s", url, exc)
    if last_error is not None:
        raise RuntimeError("Could not download HGNC complete set from primary or fallback URLs.") from last_error
    raw = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return build_hgnc_maps(raw)


def load_ncbi_gene_maps(cfg: dict[str, Any], source_dir: Path, force: bool = False) -> dict[str, object]:
    # Build lower-confidence rescue maps from NCBI Gene history and gene info
    if not cfg.get("enabled", False):
        return {}
    ncbi_dir = source_dir / "ncbi_gene"
    ensure_dirs(ncbi_dir)
    cache_path = ncbi_dir / "ncbi_gene_symbol_rescue.tsv.gz"
    map_cache_path = ncbi_dir / "ncbi_gene_symbol_rescue_maps.json"
    if map_cache_path.exists() and not force:
        LOGGER.info("Loading cached NCBI Gene rescue maps from %s", map_cache_path)
        maps = _read_ncbi_gene_map_json(map_cache_path)
        if maps is not None:
            return maps
    if cache_path.exists() and not force:
        LOGGER.info("Loading cached NCBI Gene rescue map from %s", cache_path)
        maps = _read_ncbi_gene_map_cache(cache_path)
        write_json(maps, map_cache_path)
        return maps

    info_path = ncbi_dir / "Homo_sapiens.gene_info.gz"
    history_path = ncbi_dir / "gene_history.gz"
    download_file(cfg.get("gene_info_url", NCBI_GENE_INFO_URL), info_path, force=force)
    download_file(cfg.get("gene_history_url", NCBI_GENE_HISTORY_URL), history_path, force=force)

    LOGGER.info("Parsing NCBI Gene info and history rescue tables")
    gene_info = _read_ncbi_gene_info(info_path)
    rescue_rows = _build_ncbi_rescue_rows(gene_info, history_path, bool(cfg.get("use_gene_info_synonyms", True)))
    LOGGER.info("Writing %d NCBI Gene rescue records to %s", len(rescue_rows), cache_path)
    rescue_rows.to_csv(cache_path, sep="\t", index=False, compression="gzip")
    maps = _ncbi_maps_from_rows(rescue_rows)
    write_json(maps, map_cache_path)
    return maps


def build_hgnc_maps(raw: pd.DataFrame) -> dict[str, object]:
    # Build current previous and alias symbol maps
    current: dict[str, str] = {}
    previous: dict[str, set[str]] = {}
    alias: dict[str, set[str]] = {}
    ensembl: dict[str, str] = {}

    for _, row in raw.iterrows():
        symbol = str(row.get("symbol", "")).strip()
        if not symbol:
            continue
        current[symbol] = symbol
        if row.get("ensembl_gene_id", ""):
            ensembl[symbol] = str(row["ensembl_gene_id"]).strip()
        for prev in _split_hgnc_list(row.get("prev_symbol", "")):
            previous.setdefault(prev, set()).add(symbol)
        for alias_symbol in _split_hgnc_list(row.get("alias_symbol", "")):
            alias.setdefault(alias_symbol, set()).add(symbol)

    previous_unambiguous = {key: next(iter(vals)) for key, vals in previous.items() if len(vals) == 1}
    alias_unambiguous = {key: next(iter(vals)) for key, vals in alias.items() if len(vals) == 1}
    ambiguous = {
        "previous": {key: sorted(vals) for key, vals in previous.items() if len(vals) > 1},
        "alias": {key: sorted(vals) for key, vals in alias.items() if len(vals) > 1},
    }
    return {
        "current": current,
        "previous": previous_unambiguous,
        "alias": alias_unambiguous,
        "ambiguous": ambiguous,
        "current_to_ensembl": ensembl,
        "raw": raw,
    }


def _read_ncbi_gene_map_cache(path: Path) -> dict[str, object]:
    rows = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return _ncbi_maps_from_rows(rows)


def _read_ncbi_gene_map_json(path: Path) -> dict[str, object] | None:
    try:
        with path.open("r", encoding="utf-8") as handle:
            maps = json.load(handle)
    except Exception as exc:
        LOGGER.warning("Could not read cached NCBI Gene rescue maps %s: %s", path, exc)
        return None
    if not isinstance(maps, dict) or not isinstance(maps.get("stages"), dict) or not isinstance(maps.get("ambiguous"), dict):
        LOGGER.warning("Cached NCBI Gene rescue maps have an unexpected shape: %s", path)
        return None
    return maps


def _read_ncbi_gene_info(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    tax_col = "#tax_id" if "#tax_id" in raw else "tax_id"
    raw = raw[raw[tax_col].eq("9606")].copy()
    raw["ensembl_gene_ids"] = raw.get("dbXrefs", pd.Series([""] * len(raw), index=raw.index)).map(_ensembl_ids_from_dbxrefs)
    return raw


def _build_ncbi_rescue_rows(gene_info: pd.DataFrame, history_path: Path, use_synonyms: bool) -> pd.DataFrame:
    current_by_gene_id: dict[str, dict[str, str]] = {}
    for _, row in gene_info.iterrows():
        gene_id = str(row.get("GeneID", "")).strip()
        symbol = str(row.get("Symbol", "")).strip()
        if not gene_id or not symbol or symbol == "-":
            continue
        current_by_gene_id[gene_id] = {
            "current_symbol": symbol,
            "current_gene_id": gene_id,
            "ensembl_gene_ids": str(row.get("ensembl_gene_ids", "")),
        }

    rows: list[dict[str, str]] = []
    history_iter = pd.read_csv(history_path, sep="\t", dtype=str, keep_default_na=False, chunksize=500_000)
    for chunk in history_iter:
        tax_col = "#tax_id" if "#tax_id" in chunk else "tax_id"
        chunk = chunk[chunk[tax_col].eq("9606")]
        for _, row in chunk.iterrows():
            discontinued_symbol = str(row.get("Discontinued_Symbol", "")).strip()
            gene_id = str(row.get("GeneID", "")).strip()
            if not discontinued_symbol or discontinued_symbol == "-" or not gene_id or gene_id == "-":
                continue
            current = current_by_gene_id.get(gene_id)
            if not current:
                continue
            rows.append(
                {
                    "input_symbol": discontinued_symbol,
                    "stage": "ncbi_gene_history",
                    **current,
                }
            )

    if use_synonyms:
        for _, row in gene_info.iterrows():
            current_symbol = str(row.get("Symbol", "")).strip()
            if not current_symbol or current_symbol == "-":
                continue
            current = {
                "current_symbol": current_symbol,
                "current_gene_id": str(row.get("GeneID", "")).strip(),
                "ensembl_gene_ids": str(row.get("ensembl_gene_ids", "")),
            }
            for synonym in _split_ncbi_synonyms(row.get("Synonyms", "")):
                if synonym != current_symbol:
                    rows.append({"input_symbol": synonym, "stage": "ncbi_gene_info_synonym", **current})

    columns = ["input_symbol", "stage", "current_symbol", "current_gene_id", "ensembl_gene_ids"]
    return pd.DataFrame(rows, columns=columns).drop_duplicates()


def _ncbi_maps_from_rows(rows: pd.DataFrame) -> dict[str, object]:
    stages: dict[str, dict[str, dict[str, object]]] = {}
    ambiguous: dict[str, dict[str, list[str]]] = {}
    if rows.empty:
        return {"stages": stages, "ambiguous": ambiguous}
    for stage, stage_rows in rows.groupby("stage"):
        stage_map: dict[str, dict[str, object]] = {}
        stage_ambiguous: dict[str, list[str]] = {}
        for input_symbol, symbol_rows in stage_rows.groupby("input_symbol"):
            targets = symbol_rows[["current_symbol", "current_gene_id", "ensembl_gene_ids"]].drop_duplicates()
            current_symbols = sorted(set(targets["current_symbol"]))
            if len(current_symbols) != 1:
                stage_ambiguous[str(input_symbol)] = current_symbols
                continue
            ensembl_ids = sorted(
                {
                    ensembl
                    for value in targets["ensembl_gene_ids"]
                    for ensembl in str(value).split(",")
                    if ensembl
                }
            )
            stage_map[str(input_symbol)] = {
                "current_symbol": current_symbols[0],
                "current_gene_ids": sorted(set(targets["current_gene_id"])),
                "ensembl_gene_ids": ensembl_ids,
            }
        stages[str(stage)] = stage_map
        ambiguous[str(stage)] = stage_ambiguous
    return {"stages": stages, "ambiguous": ambiguous}


def _ensembl_ids_from_dbxrefs(value: object) -> str:
    text = "" if value is None else str(value)
    matches = re.findall(r"(?:^|\|)Ensembl:(ENSG[0-9]+)", text)
    return ",".join(sorted(set(matches)))


def _split_ncbi_synonyms(value: object) -> list[str]:
    text = "" if value is None else str(value).strip()
    if not text or text == "-":
        return []
    return [item.strip() for item in text.split("|") if item.strip() and item.strip() != "-"]


def write_gene_symbol_master(hgnc_maps: dict[str, object], out_path: Path) -> None:
    raw = hgnc_maps.get("raw")
    if isinstance(raw, pd.DataFrame):
        write_tsv_gz(raw, out_path)


def load_gencode(
    source_dir: Path,
    version: str,
    supplied_gtf: Path | None = None,
    force: bool = False,
) -> pd.DataFrame:
    # Backward-compatible helper for human GENCODE annotation releases.
    version = normalize_gencode_version(version)
    if supplied_gtf:
        path = supplied_gtf
    else:
        path = source_dir / "gencode" / f"gencode.v{version}.annotation.gtf.gz"
        url = gencode_annotation_url(version)
        try:
            download_file(url, path, force=force)
        except Exception as exc:
            raise RuntimeError(
                f"Could not download GENCODE v{version} from {url}. "
                "Supply an annotation_gtf path in the target config if the release URL changed."
            ) from exc
    LOGGER.info("Parsing GENCODE v%s genes.", version)
    return parse_gtf_genes(path, version)


def load_ensembl_gtf_target(
    source_dir: Path,
    target: dict[str, Any],
    force: bool = False,
) -> tuple[pd.DataFrame, str, Path]:
    # Read an Ensembl-style target namespace from any GTF origin.
    return _load_ensembl_gtf_annotation(source_dir, target, target.get("annotation", {}), force, role="annotation")


def _load_ensembl_gtf_annotation(
    source_dir: Path,
    target: dict[str, Any],
    annotation_config: object,
    force: bool = False,
    role: str = "annotation",
) -> tuple[pd.DataFrame, str, Path]:
    annotation = target.get("annotation", {}) if isinstance(target.get("annotation", {}), dict) else {}
    if isinstance(annotation_config, dict):
        annotation = annotation_config
    legacy_version = target.get("gencode_version") if role == "annotation" else ""
    version = str(annotation.get("version") or legacy_version or "")
    source_name = str(annotation.get("source") or ("GENCODE" if legacy_version else "GTF"))
    if source_name.upper() == "GENCODE" and version:
        version = normalize_gencode_version(version)
    legacy_path = target.get("annotation_gtf") if role == "annotation" else None
    supplied_path = _optional_local_path(annotation.get("path") or legacy_path)
    if supplied_path:
        path = supplied_path
    elif annotation.get("url"):
        file_name = str(annotation.get("file_name") or Path(str(annotation["url"]).split("?", 1)[0]).name)
        if not file_name:
            file_name = f"{sanitize_id(str(target['name']))}.gtf.gz"
        path = source_dir / "targets" / sanitize_id(str(target["name"])) / sanitize_id(role) / file_name
        download_file(str(annotation["url"]), path, force=force)
    else:
        if not version:
            raise ValueError(f"Target {target.get('name', '<unnamed>')} must define an annotation path, URL, or version.")
        if source_name.upper() != "GENCODE":
            raise ValueError(
                f"Target {target.get('name', '<unnamed>')} uses annotation source {source_name!r}; "
                "provide annotation.path or annotation.url for non-GENCODE target annotations."
        )
        path = source_dir / "gencode" / f"gencode.v{version}.annotation.gtf.gz"
        url = gencode_annotation_url(version)
        try:
            download_file(url, path, force=force)
        except Exception as exc:
            raise RuntimeError(
                f"Could not download GENCODE v{version} from {url}. "
                "Supply annotation.path or annotation.url in the target config if the release URL changed."
            ) from exc
        source_name = "GENCODE"
    label = source_name
    if version:
        label = f"{source_name} v{version}"
    LOGGER.info("Parsing %s target annotation genes from %s.", label, path)
    return parse_gtf_genes(path, version, _gtf_attribute_names(annotation.get("attributes", {}))), label, path


def load_metadata_fallbacks(
    source_dir: Path,
    target: dict[str, Any],
    output_genes: set[str] | None,
    primary_annotation_ids: set[str],
    force: bool = False,
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    # Fill measured output IDs missing from the primary annotation helper.
    annotation = target.get("annotation", {}) if isinstance(target.get("annotation", {}), dict) else {}
    if "supplements" in annotation:
        raise ValueError(
            f"Target {target.get('name', '<unnamed>')} uses unsupported key annotation.supplements; "
            "use annotation.metadata_fallbacks instead."
        )
    if "annotation_supplements" in target:
        raise ValueError(
            f"Target {target.get('name', '<unnamed>')} uses unsupported key annotation_supplements; "
            "use annotation.metadata_fallbacks instead."
        )
    fallbacks = annotation.get("metadata_fallbacks", [])
    if fallbacks is None:
        fallbacks = []
    if isinstance(fallbacks, dict):
        fallbacks = [fallbacks]
    if not isinstance(fallbacks, list):
        raise TypeError(f"Target {target.get('name', '<unnamed>')} annotation.metadata_fallbacks must be a list.")

    if not fallbacks:
        return pd.DataFrame(), []

    records: list[dict[str, object]] = []
    frames: list[pd.DataFrame] = []
    output_gene_ids = {strip_ensembl_version(gene) for gene in output_genes} if output_genes else set()
    covered_ids = set(primary_annotation_ids)
    for index, spec in enumerate(fallbacks, start=1):
        if not isinstance(spec, dict):
            raise TypeError(f"Target {target.get('name', '<unnamed>')} annotation metadata fallback {index} must be a mapping.")
        mode = str(spec.get("mode", "fill_missing_output_metadata"))
        if mode != "fill_missing_output_metadata":
            raise ValueError(f"Unsupported annotation metadata fallback mode {mode!r}; use 'fill_missing_output_metadata'.")
        frame, label, path = _load_ensembl_gtf_annotation(source_dir, target, spec, force, role=f"metadata_fallback_{index}")
        fallback_ids = {strip_ensembl_version(gene) for gene in frame.get("ensembl_gene_id", pd.Series(dtype=str))}
        missing_before = output_gene_ids - covered_ids if output_gene_ids else set()
        add_ids = missing_before & fallback_ids
        selected = frame[frame["ensembl_gene_id"].isin(add_ids)].copy() if add_ids else frame.iloc[0:0].copy()
        if not selected.empty:
            selected["target_id_metadata_source"] = "metadata_fallback"
            frames.append(selected)
            covered_ids.update(add_ids)
        records.append(
            {
                "target_namespace": str(target.get("name", "")),
                "fallback_index": index,
                "mode": mode,
                "annotation": label,
                "annotation_path": str(path),
                "n_fallback_genes": len(fallback_ids),
                "n_missing_before": len(missing_before),
                "n_rows_added": len(selected),
                "n_missing_after": len(output_gene_ids - covered_ids) if output_gene_ids else 0,
            }
        )
    if frames:
        fallback_frame = pd.concat(frames, ignore_index=True).drop_duplicates()
    else:
        fallback_frame = pd.DataFrame()
    return fallback_frame, records


def normalize_gencode_version(version: object) -> str:
    text = str(version).strip()
    if text.lower().startswith("v"):
        text = text[1:]
    if not re.fullmatch(r"[0-9]+", text):
        raise ValueError(f"GENCODE version must be numeric, got {version!r}.")
    return text


def gencode_annotation_url(version: object) -> str:
    version = normalize_gencode_version(version)
    return f"{GENCODE_HUMAN_BASE_URL}/release_{version}/gencode.v{version}.annotation.gtf.gz"


def _optional_local_path(value: object) -> Path | None:
    if value in {None, "", "null"}:
        return None
    return Path(str(value))


def parse_gtf_genes(path: Path, version: str, attributes: dict[str, list[str]] | None = None) -> pd.DataFrame:
    # Extract gene records from GTF annotation
    attributes = attributes or DEFAULT_GTF_ATTRIBUTES
    rows: list[dict[str, object]] = []
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = _parse_attrs(parts[8])
            gene_id_versioned = _first_attr(attrs, attributes["gene_id"])
            rows.append(
                {
                    "gencode_version": f"v{version}",
                    "ensembl_gene_id_versioned": gene_id_versioned,
                    "ensembl_gene_id": strip_ensembl_version(gene_id_versioned),
                    "gene_symbol": _first_attr(attrs, attributes["gene_name"]),
                    "gene_biotype": _first_attr(attrs, attributes["gene_biotype"]),
                    "seqname": parts[0],
                    "start": parts[3],
                    "end": parts[4],
                    "strand": parts[6],
                    "source_gtf": str(path),
                    "target_id_metadata_source": "annotation_gtf",
                    "is_in_output_genes": False,
                }
            )
    return pd.DataFrame(rows).drop_duplicates()


def _gtf_attribute_names(config: object) -> dict[str, list[str]]:
    result = {key: list(values) for key, values in DEFAULT_GTF_ATTRIBUTES.items()}
    if not isinstance(config, dict):
        return result
    for canonical in result:
        value = config.get(canonical)
        if value is None:
            continue
        if isinstance(value, str):
            result[canonical] = [value]
        else:
            result[canonical] = [str(item) for item in value if str(item)]
    return result


def _first_attr(attrs: dict[str, str], names: list[str]) -> str:
    for name in names:
        value = attrs.get(name, "")
        if value:
            return value
    return ""


def augment_annotation_with_output_gene_symbols(
    annotation: pd.DataFrame,
    output_gene_table: pd.DataFrame,
    output_genes: set[str] | None,
    source_path: Path | None = None,
) -> pd.DataFrame:
    # Add dataset-provided symbol metadata for output genes absent from the annotation GTF.
    if not output_genes or output_gene_table.empty or "gene_symbol" not in output_gene_table:
        return annotation
    annotation = annotation.copy()
    if "target_id_metadata_source" not in annotation:
        annotation["target_id_metadata_source"] = "annotation_gtf"
    existing_ids = {strip_ensembl_version(gene) for gene in annotation.get("ensembl_gene_id", pd.Series(dtype=str))}
    rows: list[dict[str, object]] = []
    seen: set[str] = set()
    for _, row in output_gene_table.iterrows():
        stable = strip_ensembl_version(row.get("ensembl_gene_id", ""))
        symbol = str(row.get("gene_symbol", "")).strip()
        if not stable or not symbol or stable in existing_ids or stable in seen:
            continue
        seen.add(stable)
        record = {column: "" for column in annotation.columns}
        record.update(
            {
                "gencode_version": "",
                "ensembl_gene_id_versioned": str(row.get("input_gene_id", "") or stable),
                "ensembl_gene_id": stable,
                "gene_symbol": symbol,
                "gene_biotype": "",
                "source_gtf": f"output_genes:{source_path}" if source_path else "output_genes",
                "target_id_metadata_source": "output_genes",
                "is_in_output_genes": True,
            }
        )
        rows.append(record)
    if not rows:
        return annotation
    return pd.concat([annotation, pd.DataFrame(rows, columns=annotation.columns)], ignore_index=True)


def mark_output_genes(gencode: pd.DataFrame, output_genes: set[str] | None) -> pd.DataFrame:
    # Mark genes present in the configured output gene set.
    gencode = gencode.copy()
    if output_genes:
        gencode["is_in_output_genes"] = gencode["ensembl_gene_id"].isin(output_genes)
    else:
        gencode["is_in_output_genes"] = True
    return gencode


def write_mapping(gencode: pd.DataFrame, out_path: Path) -> None:
    write_tsv_gz(gencode, out_path)


def read_output_genes(
    path: Path | None,
    id_column: str | None = None,
    ids: list[str] | tuple[str, ...] | set[str] | None = None,
    symbol_column: str | None = None,
) -> tuple[set[str] | None, pd.DataFrame, str]:
    # Read Ensembl IDs from a count matrix or gene list.
    if ids is not None:
        genes = [str(gene) for gene in ids]
        symbols = [""] * len(genes)
    elif path is None:
        return None, pd.DataFrame(columns=["input_gene_id", "ensembl_gene_id", "gene_symbol", "id_type"]), "none"
    elif path.suffix == ".parquet":
        columns = []
        if id_column and id_column not in {"index", "__index__"}:
            columns.append(id_column)
        if symbol_column:
            columns.append(symbol_column)
        if columns:
            try:
                df = pd.read_parquet(path, columns=list(dict.fromkeys(columns)))
            except Exception as exc:
                raise ValueError(f"Output gene ID column or symbol column was not found in {path}.") from exc
        else:
            df = pd.read_parquet(path, columns=[])
        genes = _output_id_values(df, id_column, path, allow_index=True)
        symbols = _output_symbol_values(df, symbol_column, path, len(genes))
    elif id_column or symbol_column or _looks_tabular(path):
        df = _read_output_id_table(path)
        genes = _output_id_values(df, id_column, path, allow_index=False)
        symbols = _output_symbol_values(df, symbol_column, path, len(genes))
    else:
        with open_text(path) as handle:
            genes = [line.strip().split()[0] for line in handle if line.strip()]
        symbols = [""] * len(genes)

    rows = []
    stable = set()
    id_types: list[str] = []
    for gene, symbol in zip(genes, symbols):
        stripped = strip_ensembl_version(gene)
        row_id_type = "symbol"
        if stripped.startswith("ENSG"):
            row_id_type = "ensembl_versioned" if stripped != gene else "ensembl_stable"
            stable.add(stripped)
        id_types.append(row_id_type)
        rows.append({"input_gene_id": gene, "ensembl_gene_id": stripped, "gene_symbol": str(symbol).strip(), "id_type": row_id_type})
    if genes and not stable:
        raise ValueError("output_genes and gene_filter must provide Ensembl IDs in id_column/ids; symbol_column is optional metadata.")
    unique_id_types = set(id_types)
    id_type = "none" if not unique_id_types else next(iter(unique_id_types)) if len(unique_id_types) == 1 else "mixed"
    return stable if stable else None, pd.DataFrame(rows).drop_duplicates(), id_type


def _output_id_values(df: pd.DataFrame, column: str | None, path: Path, allow_index: bool) -> list[str]:
    if column:
        if column in {"index", "__index__"} and allow_index:
            return df.index.astype(str).tolist()
        if column not in df.columns:
            raise ValueError(f"Output gene ID column {column!r} was not found in {path}.")
        return df[column].astype(str).tolist()
    col = first_existing_column(list(df.columns), GENE_COLUMNS)
    if col:
        return df[col].astype(str).tolist()
    if allow_index:
        return df.index.astype(str).tolist()
    return df[df.columns[0]].astype(str).tolist()


def _output_symbol_values(df: pd.DataFrame, symbol_column: str | None, path: Path, n_rows: int) -> list[str]:
    if not symbol_column:
        return [""] * n_rows
    if symbol_column not in df:
        raise ValueError(f"Output gene symbol column {symbol_column!r} was not found in {path}.")
    return df[symbol_column].fillna("").astype(str).tolist()


def _read_output_id_table(path: Path) -> pd.DataFrame:
    name = path.name.lower()
    if name.endswith((".tsv", ".tsv.gz", ".tab", ".tab.gz")):
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if name.endswith((".csv", ".csv.gz")):
        return pd.read_csv(path, sep=",", dtype=str, keep_default_na=False)
    return pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)


def _normalize_msigdb(
    raw: pd.DataFrame,
    tag: str,
    collection: str,
    subcollection: str,
    family: str,
    aspect: str,
) -> pd.DataFrame:
    # Convert MSigDB rows to the common source schema
    if raw.empty:
        return empty_term_gene_frame()
    term_ids = "MSIGDB_" + tag + "__" + raw["gs_name"].map(sanitize_id)
    return pd.DataFrame(
        {
            "term_id": term_ids,
            "original_name": raw["gs_name"],
            "display_name": raw["gs_name"],
            "description": _col(raw, "gs_description"),
            "source": "MSigDB",
            "source_tag": tag,
            "collection": _col(raw, "gs_collection", collection).replace("", collection),
            "subcollection": _col(raw, "gs_subcollection", subcollection).replace("", subcollection),
            "family": family,
            "aspect": aspect,
            "db_version": _col(raw, "db_version"),
            "source_url": "https://www.gsea-msigdb.org/gsea/msigdb",
            "source_license_note": "MSigDB license and attribution terms apply.",
            "gene_symbol": raw["gene_symbol"],
            "gene_ensembl_from_source": _col(raw, "ensembl_gene"),
            "metadata_json": "{}",
        }
    )


def _normalize_enrichr_gmt(
    path: Path,
    library: str,
    source_tag: str,
    family: str,
    aspect: str,
    metadata: dict[str, object] | None = None,
) -> pd.DataFrame:
    # Convert Enrichr GMT rows to the common source schema
    rows: list[dict[str, str]] = []
    safe_library = sanitize_id(library)
    metadata_json = json.dumps(metadata or {}, sort_keys=True)
    for record in parse_gmt(path):
        term_name = str(record["term_id"])
        term_id = f"ENRICHR_{safe_library}__{sanitize_id(term_name)}"
        for gene in record["genes"]:
            rows.append(
                {
                    "term_id": term_id,
                    "original_name": term_name,
                    "display_name": term_name,
                    "description": str(record["description"]),
                    "source": "Enrichr",
                    "source_tag": source_tag,
                    "collection": library,
                    "subcollection": "",
                    "family": family,
                    "aspect": aspect,
                    "db_version": "",
                    "source_url": "https://maayanlab.cloud/Enrichr",
                    "source_license_note": "Enrichr library provenance and upstream licenses apply.",
                    "gene_symbol": gene,
                    "gene_ensembl_from_source": "",
                    "metadata_json": metadata_json,
                }
            )
    return pd.DataFrame(rows, columns=CANONICAL_COLUMNS)


def _col(raw: pd.DataFrame, name: str, default: str = "") -> pd.Series:
    if name in raw:
        return raw[name].fillna("")
    return pd.Series([default] * len(raw), index=raw.index)


def _extract_library_names(stats: object) -> list[str]:
    if isinstance(stats, dict) and "statistics" in stats:
        stats = stats["statistics"]
    names: list[str] = []
    if isinstance(stats, list):
        for item in stats:
            if isinstance(item, dict):
                name = item.get("libraryName") or item.get("library_name") or item.get("name")
                if name:
                    names.append(str(name))
    return sorted(set(names))


def _latest_year(name: str) -> int:
    years = [int(year) for year in re.findall(r"(20[0-9]{2})", name)]
    return max(years) if years else 0


def _split_hgnc_list(value: object) -> list[str]:
    if value is None:
        return []
    text = str(value).strip()
    if not text:
        return []
    return [item.strip() for item in text.split("|") if item.strip()]


def _parse_attrs(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        else:
            parts = item.split(None, 1)
            if len(parts) != 2:
                continue
            key, value = parts
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def _looks_tabular(path: Path) -> bool:
    with open_text(path) as handle:
        first = handle.readline()
    return "\t" in first or "," in first
