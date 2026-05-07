from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
import zipfile
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
    utc_now,
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

ENRICHR_LOGICAL = {
    "disgenet": ("DISGENET", "disease_phenotype", "disease"),
    "clinvar": ("CLINVAR", "disease_phenotype", "disease"),
    "gwas_catalog": ("GWAS_CATALOG", "disease_phenotype", "disease"),
    "jensen_diseases": ("JENSEN_DISEASES_CURATED", "disease_phenotype", "disease"),
    "depmap": ("DEPMAP", "cancer_dependency_state", "dependency_signature"),
    "ccle": ("CCLE", "cancer_dependency_state", "cancer_signature"),
    "nci60": ("NCI60", "cancer_dependency_state", "cancer_signature"),
}

HGNC_URL = "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt"
HGNC_FALLBACK_URLS = [
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
]

GENCODE_URLS = {
    "47": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz",
    "29": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz",
}

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


def load_source_memberships(cfg: dict[str, Any]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    source_dir = Path(cfg.get("project", {}).get("source_dir", "data/alien_sources"))
    for spec in cfg.get("sources", []):
        if not spec.get("enabled", True):
            continue
        frames.append(_load_source_spec(spec, source_dir))
    return concat_frames(frames)


def _load_source_spec(spec: dict[str, Any], source_dir: Path) -> pd.DataFrame:
    source_type = str(spec.get("type", "")).strip()
    if source_type == "canonical_tsv":
        return concat_frames([read_canonical_memberships(path, spec) for path in _source_paths(spec, source_dir)])
    if source_type == "msigdb_cache":
        path = Path(spec.get("path") or source_dir / "msigdb")
        return read_msigdb_cache(path, bool(spec.get("include_c4_cm", False)), spec)
    if source_type == "msigdb_remote":
        return read_msigdb_remote(spec, source_dir)
    if source_type == "msigdb_tsv":
        return concat_frames([read_msigdb_like_table(path, spec) for path in _source_paths(spec, source_dir)])
    if source_type in {"symbol_gmt", "enrichr_gmt"}:
        return concat_frames([read_symbol_gmt_source(path, spec) for path in _source_paths(spec, source_dir)])
    raise ValueError(f"Unsupported source type: {source_type or '<empty>'}")


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


def read_msigdb_cache(path: Path, include_c4_cm: bool, spec: dict[str, Any] | None = None) -> pd.DataFrame:
    spec = spec or {}
    frames: list[pd.DataFrame] = []
    for tag, (collection, subcollection, family, aspect) in MSIGDB_SOURCES.items():
        if tag == "C4_CM" and not include_c4_cm:
            continue
        file_path = path / f"msigdbr_{tag}.tsv.gz"
        if not file_path.exists():
            if spec.get("optional", False):
                LOGGER.warning("Optional MSigDB cache file is missing: %s", file_path)
                continue
            raise FileNotFoundError(f"Missing MSigDB cache file: {file_path}")
        frames.append(read_msigdb_like_table(file_path, {"source_tag": tag, "collection": collection, "subcollection": subcollection, "family": family, "aspect": aspect}))
    return concat_frames(frames)


def read_msigdb_remote(spec: dict[str, Any], source_dir: Path) -> pd.DataFrame:
    cache_dir = Path(spec.get("cache_dir") or source_dir / "msigdb_remote")
    release = _msigdb_release_info(spec)
    release_dir = ensure_msigdb_remote_cache(
        cache_dir=cache_dir,
        release=release,
        force=bool(spec.get("force", False)),
        timeout=int(spec.get("timeout_seconds", 600)),
    )
    summary = _read_rds_dataframe(release_dir / str(release["summary_rds"]))
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
    return concat_frames(frames)


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
    raw = raw.copy()
    raw["gene_symbol"] = raw["db_gene_symbol"]
    raw["ensembl_gene"] = raw.get("db_ensembl_gene", "")
    raw["_alien_source_tag"] = raw.apply(_msigdb_remote_row_source_tag, axis=1)
    frames = []
    for tag, tag_frame in raw.groupby("_alien_source_tag", sort=True):
        family = str(spec.get("family") or _family_for_msigdb_tag(tag))
        aspect = str(spec.get("aspect") or _aspect_for_msigdb_tag(tag))
        collection = str(tag_frame["gs_collection"].iloc[0]) if "gs_collection" in tag_frame else ""
        subcollection = str(tag_frame["gs_subcollection"].iloc[0]) if "gs_subcollection" in tag_frame else ""
        frames.append(_normalize_msigdb(tag_frame, tag, collection, subcollection, family, aspect))
    return concat_frames(frames)


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


def ensure_msigdb_sources(
    source_dir: Path,
    include_c4_cm: bool,
    min_version: str,
    force: bool = False,
    r_script: Path | None = None,
) -> None:
    # Fetch MSigDB cache through the R helper when needed
    msigdb_dir = source_dir / "msigdb"
    wanted = [tag for tag in MSIGDB_SOURCES if include_c4_cm or tag != "C4_CM"]
    missing = [tag for tag in wanted if not (msigdb_dir / f"msigdbr_{tag}.tsv.gz").exists()]
    if not missing and not force:
        return

    r_script = r_script or Path("scripts/fetch_msigdb.R")
    cmd = ["Rscript", str(r_script), str(msigdb_dir), str(include_c4_cm).upper(), min_version]
    LOGGER.info("Fetching MSigDB through msigdbr.")
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "Rscript is unavailable. Install R and msigdbr >= "
            f"{min_version}, then rerun the GMT builder."
        ) from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            "MSigDB fetch failed. Install/update msigdbr in R, for example: "
            "install.packages('msigdbr'), and ensure msigdbr meets the configured version."
        ) from exc


def read_msigdb(source_dir: Path, include_c4_cm: bool) -> pd.DataFrame:
    # Read cached MSigDB collections
    frames: list[pd.DataFrame] = []
    msigdb_dir = source_dir / "msigdb"
    for tag, (collection, subcollection, family, aspect) in MSIGDB_SOURCES.items():
        if tag == "C4_CM" and not include_c4_cm:
            continue
        path = msigdb_dir / f"msigdbr_{tag}.tsv.gz"
        if not path.exists():
            raise FileNotFoundError(f"Missing MSigDB cache file: {path}")
        raw = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        if "gs_name" not in raw or "gene_symbol" not in raw:
            raise ValueError(f"{path} must contain gs_name and gene_symbol columns.")
        frames.append(_normalize_msigdb(raw, tag, collection, subcollection, family, aspect))
    return concat_frames(frames)


def load_enrichr(cfg: dict, source_dir: Path, force: bool = False) -> tuple[pd.DataFrame, dict, list[str]]:
    # Select and fetch optional Enrichr libraries
    enrichr_dir = source_dir / "enrichr"
    ensure_dirs(enrichr_dir)
    warnings: list[str] = []
    provenance: dict = {"selected_libraries": {}, "candidate_libraries": {}, "warnings": []}
    stats_path = enrichr_dir / "datasetStatistics.json"

    try:
        if force or not stats_path.exists():
            LOGGER.info("Downloading Enrichr library metadata.")
            response = requests.get(cfg["dataset_statistics_endpoint"], timeout=(15, 120))
            response.raise_for_status()
            stats_path.write_text(response.text, encoding="utf-8")
        stats = json.loads(stats_path.read_text(encoding="utf-8"))
        names = _extract_library_names(stats)
    except Exception as exc:
        message = f"Enrichr metadata unavailable: {exc}"
        if cfg.get("optional", True):
            warnings.append(message)
            provenance["warnings"].append(message)
            return pd.DataFrame(), provenance, warnings
        raise

    frames: list[pd.DataFrame] = []
    for group_name, group_patterns in cfg.get("library_patterns", {}).items():
        for logical, patterns in group_patterns.items():
            source_tag, family, aspect = ENRICHR_LOGICAL.get(logical, (sanitize_id(logical).upper(), group_name, group_name))
            selected, candidates = _select_library(names, patterns)
            provenance["candidate_libraries"][logical] = candidates
            provenance["selected_libraries"][logical] = selected
            if not selected:
                warnings.append(f"Optional Enrichr library missing for {logical}.")
                continue
            try:
                gmt_path = _download_library(cfg, enrichr_dir, selected, force)
                frames.append(_normalize_enrichr_gmt(gmt_path, selected, source_tag, family, aspect))
            except Exception as exc:
                message = f"Optional Enrichr library {selected} failed: {exc}"
                if cfg.get("optional", True):
                    warnings.append(message)
                    continue
                raise

    provenance["source_url"] = cfg.get("base_url", "")
    provenance["download_timestamp_utc"] = utc_now()
    provenance["local_metadata_path"] = str(stats_path)
    return concat_frames(frames), provenance, warnings


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
    if cache_path.exists() and not force:
        LOGGER.info("Loading cached NCBI Gene rescue map from %s", cache_path)
        return _read_ncbi_gene_map_cache(cache_path)

    info_path = ncbi_dir / "Homo_sapiens.gene_info.gz"
    history_path = ncbi_dir / "gene_history.gz"
    download_file(cfg.get("gene_info_url", NCBI_GENE_INFO_URL), info_path, force=force)
    download_file(cfg.get("gene_history_url", NCBI_GENE_HISTORY_URL), history_path, force=force)

    LOGGER.info("Parsing NCBI Gene info and history rescue tables")
    gene_info = _read_ncbi_gene_info(info_path)
    rescue_rows = _build_ncbi_rescue_rows(gene_info, history_path, bool(cfg.get("use_gene_info_synonyms", True)))
    LOGGER.info("Writing %d NCBI Gene rescue records to %s", len(rescue_rows), cache_path)
    rescue_rows.to_csv(cache_path, sep="\t", index=False, compression="gzip")
    return _ncbi_maps_from_rows(rescue_rows)


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
    # Backward-compatible helper for the built-in GENCODE URLs.
    if supplied_gtf:
        path = supplied_gtf
    else:
        path = source_dir / "gencode" / f"gencode.v{version}.annotation.gtf.gz"
        try:
            download_file(GENCODE_URLS[version], path, force=force)
        except Exception as exc:
            raise RuntimeError(
                f"Could not download GENCODE v{version} from {GENCODE_URLS[version]}. "
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
    annotation = target.get("annotation", {}) if isinstance(target.get("annotation", {}), dict) else {}
    version = str(annotation.get("version") or target.get("gencode_version") or "")
    source_name = str(annotation.get("source") or ("GENCODE" if target.get("gencode_version") else "GTF"))
    supplied_path = _optional_local_path(annotation.get("path") or target.get("annotation_gtf"))
    if supplied_path:
        path = supplied_path
    elif annotation.get("url"):
        file_name = str(annotation.get("file_name") or Path(str(annotation["url"]).split("?", 1)[0]).name)
        if not file_name:
            file_name = f"{sanitize_id(str(target['name']))}.gtf.gz"
        path = source_dir / "targets" / sanitize_id(str(target["name"])) / file_name
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
        try:
            download_file(GENCODE_URLS[version], path, force=force)
        except Exception as exc:
            raise RuntimeError(
                f"Could not download GENCODE v{version} from {GENCODE_URLS[version]}. "
                "Supply annotation.path or annotation.url in the target config if the release URL changed."
            ) from exc
        source_name = "GENCODE"
    label = source_name
    if version:
        label = f"{source_name} v{version}"
    LOGGER.info("Parsing %s target annotation genes from %s.", label, path)
    return parse_gtf_genes(path, version), label, path


def _optional_local_path(value: object) -> Path | None:
    if value in {None, "", "null"}:
        return None
    return Path(str(value))


def parse_gtf_genes(path: Path, version: str) -> pd.DataFrame:
    # Extract gene records from GTF annotation
    rows: list[dict[str, object]] = []
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = _parse_attrs(parts[8])
            gene_id_versioned = attrs.get("gene_id", "")
            rows.append(
                {
                    "gencode_version": f"v{version}",
                    "ensembl_gene_id_versioned": gene_id_versioned,
                    "ensembl_gene_id": strip_ensembl_version(gene_id_versioned),
                    "gene_symbol": attrs.get("gene_name", ""),
                    "gene_biotype": attrs.get("gene_type") or attrs.get("gene_biotype", ""),
                    "seqname": parts[0],
                    "start": parts[3],
                    "end": parts[4],
                    "strand": parts[6],
                    "source_gtf": str(path),
                    "is_in_expression_universe": False,
                }
            )
    return pd.DataFrame(rows).drop_duplicates()


def mark_universe(gencode: pd.DataFrame, universe: set[str] | None) -> pd.DataFrame:
    # Mark genes present in the expression matrix
    gencode = gencode.copy()
    if universe:
        gencode["is_in_expression_universe"] = gencode["ensembl_gene_id"].isin(universe)
    else:
        gencode["is_in_expression_universe"] = True
    return gencode


def write_mapping(gencode: pd.DataFrame, out_path: Path) -> None:
    write_tsv_gz(gencode, out_path)


def read_gene_universe(path: Path | None) -> tuple[set[str] | None, pd.DataFrame, str]:
    # Read Ensembl IDs from a count matrix or gene list
    if path is None:
        return None, pd.DataFrame(columns=["input_gene_id", "ensembl_gene_id", "id_type"]), "none"

    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
        col = first_existing_column(list(df.columns), GENE_COLUMNS)
        genes = df[col].astype(str).tolist() if col else df.index.astype(str).tolist()
    elif _looks_tabular(path):
        df = pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)
        col = first_existing_column(list(df.columns), GENE_COLUMNS) or df.columns[0]
        genes = df[col].astype(str).tolist()
    else:
        with open_text(path) as handle:
            genes = [line.strip().split()[0] for line in handle if line.strip()]

    rows = []
    stable = set()
    id_type = "symbol"
    for gene in genes:
        stripped = strip_ensembl_version(gene)
        if stripped.startswith("ENSG"):
            id_type = "ensembl_versioned" if stripped != gene else "ensembl_stable"
            stable.add(stripped)
        rows.append({"input_gene_id": gene, "ensembl_gene_id": stripped, "id_type": id_type})
    return stable if stable else None, pd.DataFrame(rows).drop_duplicates(), id_type


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


def _normalize_enrichr_gmt(path: Path, library: str, source_tag: str, family: str, aspect: str) -> pd.DataFrame:
    # Convert Enrichr GMT rows to the common source schema
    rows: list[dict[str, str]] = []
    safe_library = sanitize_id(library)
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
                    "metadata_json": "{}",
                }
            )
    return pd.DataFrame(rows)


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


def _select_library(names: list[str], patterns: list[str]) -> tuple[str | None, list[str]]:
    exact_patterns = [re.compile(pattern) for pattern in patterns if pattern.startswith("^") and pattern.endswith("$")]
    all_patterns = [re.compile(pattern, re.I) for pattern in patterns]
    candidates = sorted({name for name in names for pattern in all_patterns if pattern.search(name)})
    exact = [name for name in candidates for pattern in exact_patterns if pattern.match(name)]
    if exact:
        return sorted(exact)[-1], candidates
    if not candidates:
        return None, []
    return sorted(candidates, key=lambda name: (-_latest_year(name), len(name), name))[0], candidates


def _latest_year(name: str) -> int:
    years = [int(year) for year in re.findall(r"(20[0-9]{2})", name)]
    return max(years) if years else 0


def _download_library(cfg: dict, enrichr_dir: Path, library: str, force: bool) -> Path:
    path = enrichr_dir / f"{sanitize_id(library)}.gmt"
    if path.exists() and not force:
        return path
    url = cfg["library_download_endpoint"].format(library=quote(library))
    LOGGER.info("Downloading Enrichr library %s.", library)
    response = requests.get(url, timeout=(15, 60))
    response.raise_for_status()
    path.write_text(response.text, encoding="utf-8")
    return path


def _split_hgnc_list(value: object) -> list[str]:
    if value is None:
        return []
    text = str(value).strip()
    if not text:
        return []
    return [item.strip() for item in text.split("|") if item.strip()]


def _parse_attrs(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for key, value in re.findall(r'([A-Za-z0-9_]+)\s+"([^"]*)"', text):
        attrs[key] = value
    return attrs


def _looks_tabular(path: Path) -> bool:
    with open_text(path) as handle:
        first = handle.readline()
    return "\t" in first or "," in first
