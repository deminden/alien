"""MSigDB: remote msigdbr releases (Zenodo), local msigdbr caches and msigdbr-like tables."""

from __future__ import annotations

import logging
import re
import shutil
import zipfile
from pathlib import Path
from typing import Final

import pandas as pd

from ..cache import ResourceCache, derived_cache_path
from ..config.schema import MsigdbCacheSource, MsigdbRemoteSource, MsigdbTsvSource
from ..errors import ConfigError, ResourceError, SourceError
from ..model import Family, sanitize_term_id
from ..schema import CANONICAL_COLUMNS
from .frames import column_or, concat_memberships, empty_memberships
from .local import source_paths
from .registry import register

LOGGER = logging.getLogger("alien")
MSIGDB_URL: Final = "https://www.gsea-msigdb.org/gsea/msigdb"
MSIGDB_LICENSE: Final = "MSigDB license and attribution terms apply."

#: source tag -> (collection, subcollection, family, aspect)
MSIGDB_SOURCES: Final = {
    "REACTOME": ("C2", "CP:REACTOME", Family.BIOLOGY_PROCESS_PATHWAY, "pathway"),
    "WIKIPATHWAYS": ("C2", "CP:WIKIPATHWAYS", Family.BIOLOGY_PROCESS_PATHWAY, "pathway"),
    "KEGG_MEDICUS": ("C2", "CP:KEGG_MEDICUS", Family.BIOLOGY_PROCESS_PATHWAY, "pathway"),
    "GOBP": ("C5", "GO:BP", Family.BIOLOGY_PROCESS_PATHWAY, "biological_process"),
    "GOMF": ("C5", "GO:MF", Family.BIOLOGY_FUNCTION_LOCATION, "molecular_function"),
    "GOCC": ("C5", "GO:CC", Family.BIOLOGY_FUNCTION_LOCATION, "cellular_component"),
    "HPO": ("C5", "HPO", Family.DISEASE_PHENOTYPE, "phenotype"),
    "C6": ("C6", "", Family.CANCER_DEPENDENCY_STATE, "oncogenic_signature"),
    "C9": ("C9", "", Family.CANCER_DEPENDENCY_STATE, "dependency_signature"),
    "C4_3CA": ("C4", "3CA", Family.CANCER_DEPENDENCY_STATE, "cancer_signature"),
    "C4_CGN": ("C4", "CGN", Family.CANCER_DEPENDENCY_STATE, "cancer_signature"),
    "C4_CM": ("C4", "CM", Family.CANCER_DEPENDENCY_STATE, "cancer_signature"),
}
#: Built-in release table; new releases can be configured with ``release:`` without a code change.
MSIGDB_RELEASES: Final = {
    "2026.1": {
        "zip_url": "https://zenodo.org/records/18968178/files/msigdb.2026.1.zip?download=1",
        "zip_md5": "512ba99c6827141a9d471972b812d4ac",
        "zip_name": "msigdb.2026.1.zip",
        "summary_rds": "msigdb.2026.1.summary.rds",
    }
}


def normalize_msigdb(
    raw: pd.DataFrame, tag: str, collection: str, subcollection: str, family: str, aspect: str
) -> pd.DataFrame:
    """msigdbr-shaped rows (gs_name, gene_symbol, ...) to canonical memberships."""
    if raw.empty:
        return empty_memberships()
    return pd.DataFrame(
        {
            "term_id": "MSIGDB_" + tag + "__" + raw["gs_name"].map(sanitize_term_id),
            "original_name": raw["gs_name"],
            "display_name": raw["gs_name"],
            "description": column_or(raw, "gs_description"),
            "source": "MSigDB",
            "source_tag": str(tag),
            "collection": column_or(raw, "gs_collection", collection).replace("", collection),
            "subcollection": column_or(raw, "gs_subcollection", subcollection).replace("", subcollection),
            "family": str(family),
            "aspect": str(aspect),
            "db_version": column_or(raw, "db_version"),
            "source_url": MSIGDB_URL,
            "source_license_note": MSIGDB_LICENSE,
            "gene_symbol": raw["gene_symbol"],
            "gene_ensembl_from_source": column_or(raw, "ensembl_gene"),
            "metadata_json": "{}",
        },
        columns=list(CANONICAL_COLUMNS),
    )


@register("msigdb_tsv")
def load_msigdb_tsv(spec: MsigdbTsvSource, cache: ResourceCache) -> pd.DataFrame:
    return concat_memberships(read_msigdb_like_table(path, spec) for path in source_paths(spec, cache.root))


def read_msigdb_like_table(path: Path, spec: MsigdbTsvSource) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    name_col, symbol_col = spec.term_name_column, spec.symbol_column
    if name_col not in raw or symbol_col not in raw:
        raise SourceError(f"{path} must contain {name_col} and {symbol_col} columns.")
    renamed = raw.rename(columns={name_col: "gs_name", symbol_col: "gene_symbol"}).copy()
    if "ensembl_gene" not in renamed:
        gene_id_col = spec.gene_id_column
        renamed["ensembl_gene"] = renamed[gene_id_col] if gene_id_col and gene_id_col in renamed else ""
    return normalize_msigdb(
        renamed,
        spec.source_tag or spec.name or "MSIGDB",
        spec.collection or "",
        spec.subcollection or "",
        spec.family or Family.BIOLOGY_PROCESS_PATHWAY,
        spec.aspect or "library",
    )


@register("msigdb_cache")
def load_msigdb_cache(spec: MsigdbCacheSource, cache: ResourceCache) -> pd.DataFrame:
    path = spec.path or cache.root / "msigdb"
    frames = []
    for tag, (collection, subcollection, family, aspect) in MSIGDB_SOURCES.items():
        if tag == "C4_CM" and not spec.include_c4_cm:
            continue
        file_path = path / f"msigdbr_{tag}.tsv.gz"
        if not file_path.exists():
            raise SourceError(f"Missing MSigDB cache file: {file_path}")
        table_spec = MsigdbTsvSource(
            type="msigdb_tsv",
            source_tag=tag,
            collection=collection,
            subcollection=subcollection,
            family=family,
            aspect=aspect,
        )
        frames.append(read_msigdb_like_table(file_path, table_spec))
    return concat_memberships(frames)


def release_info(spec: MsigdbRemoteSource) -> dict[str, str]:
    custom = {key: str(value) for key, value in spec.release.model_dump(exclude_none=True).items()}
    # A custom archive URL starts a new release record; built-in names and checksums must not leak into it.
    release = {} if "zip_url" in custom else dict(MSIGDB_RELEASES.get(spec.version, {}))
    release.update(custom)
    if not release and spec.zip_url and spec.zip_md5:
        release = {"zip_url": spec.zip_url, "zip_md5": spec.zip_md5}
    if "zip_url" not in release or not ({"zip_md5", "zip_sha256"} & set(release)):
        raise ConfigError(
            f"Unknown MSigDB release {spec.version!r}; provide release.zip_url and release.zip_sha256 (or zip_md5)."
        )
    zip_name = release.get("zip_name") or re.sub(r"\?.*$", "", Path(release["zip_url"]).name)
    release["zip_name"] = zip_name
    release.setdefault("summary_rds", zip_name.replace(".zip", ".summary.rds"))
    return release


def _safe_extract(archive: zipfile.ZipFile, destination: Path) -> None:
    root = destination.resolve()
    for member in archive.namelist():
        if not (destination / member).resolve().is_relative_to(root):
            raise ResourceError(f"Refusing to extract {member!r} outside {destination}.")
    archive.extractall(destination)


def ensure_release(spec: MsigdbRemoteSource, cache: ResourceCache) -> Path:
    """The extracted release directory (downloading and verifying the archive when needed)."""
    release = release_info(spec)
    cache_dir = spec.cache_dir or cache.root / "msigdb_remote"
    zip_path = cache_dir / release["zip_name"]
    release_dir = cache_dir / zip_path.stem
    force = cache.force or spec.force
    summary = release_dir / release["summary_rds"]
    if force and release_dir.exists() and cache.network_allowed:
        shutil.rmtree(release_dir)
    if zip_path.exists() or not summary.exists() or force:
        cache.fetch(
            f"msigdb:{release['zip_name']}",
            release["zip_url"],
            zip_path,
            sha256=release.get("zip_sha256"),
            md5=release.get("zip_md5"),
            timeout=spec.timeout_seconds,
            upstream_version=spec.version,
        )
    if not summary.exists():
        with cache.file_lock(release_dir), zipfile.ZipFile(zip_path) as archive:
            release_dir.mkdir(parents=True, exist_ok=True)
            _safe_extract(archive, release_dir)
        nested = sorted(release_dir.rglob(release["summary_rds"]))
        if not nested:
            raise ResourceError(f"Expected MSigDB summary {release['summary_rds']} was not found in {zip_path}.")
        summary = nested[0]
    return summary.parent


def read_rds_frame(path: Path) -> pd.DataFrame:
    try:
        import rdata  # noqa: PLC0415 - optional heavy import, only needed for remote MSigDB releases
    except ImportError as exc:
        raise ResourceError("Reading MSigDB remote releases requires the Python package 'rdata'.") from exc
    obj = rdata.read_rds(path)
    if not isinstance(obj, pd.DataFrame):
        raise ResourceError(f"Expected an RDS data frame in {path}, got {type(obj).__name__}.")
    return obj.fillna("").astype(str)


def _matches_subcollection(values: pd.Series, subcollection: str) -> pd.Series:
    text = values.astype(str)
    return text.eq(subcollection) | text.str.replace(r".*:", "", regex=True).eq(subcollection)


def _filter_rows(frame: pd.DataFrame, spec: MsigdbRemoteSource, *, species_column: str) -> pd.DataFrame:
    if species_column in frame:
        frame = frame[frame[species_column].astype(str).str.upper().eq(spec.db_species.upper())]
    if spec.collection and "gs_collection" in frame:
        frame = frame[frame["gs_collection"].astype(str).eq(spec.collection)]
    if spec.subcollection and "gs_subcollection" in frame:
        frame = frame[_matches_subcollection(frame["gs_subcollection"], spec.subcollection)]
    return frame


def _read_summary(release_dir: Path, release: dict[str, str], force: bool) -> pd.DataFrame:
    rds = release_dir / release["summary_rds"]
    summary_cache = release_dir / "alien_cache" / f"{rds.stem}.parquet"
    if summary_cache.exists() and not force:
        return pd.read_parquet(summary_cache).fillna("").astype(str)
    summary = read_rds_frame(rds)
    summary_cache.parent.mkdir(parents=True, exist_ok=True)
    summary.to_parquet(summary_cache, index=False)
    return summary


def remote_source_tag(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return "MSIGDB"
    tail = text.split(":")[-1]
    return "GO" + tail if tail in {"BP", "MF", "CC"} else sanitize_term_id(tail).upper()


def normalize_remote_frame(raw: pd.DataFrame, spec: MsigdbRemoteSource) -> pd.DataFrame:
    if "db_gene_symbol" not in raw or "gs_name" not in raw:
        raise SourceError("MSigDB remote RDS must contain db_gene_symbol and gs_name columns.")
    raw = _filter_rows(raw, spec, species_column="db_target_species").copy()
    if raw.empty:
        raise SourceError("No MSigDB remote rows matched the requested filters.")
    raw["gene_symbol"] = raw["db_gene_symbol"]
    raw["ensembl_gene"] = column_or(raw, "db_ensembl_gene")
    subcollection = column_or(raw, "gs_subcollection").astype(str).str.strip()
    collection = column_or(raw, "gs_collection").astype(str).str.strip()
    tags = subcollection.where(subcollection.ne(""), collection).replace("", "MSIGDB").map(remote_source_tag)
    frames = []
    for tag, tag_frame in raw.groupby(tags, sort=True):
        default = MSIGDB_SOURCES.get(str(tag), ("", "", Family.BIOLOGY_PROCESS_PATHWAY, "library"))
        frames.append(
            normalize_msigdb(
                tag_frame,
                str(tag),
                str(tag_frame["gs_collection"].iloc[0]) if "gs_collection" in tag_frame else "",
                str(tag_frame["gs_subcollection"].iloc[0]) if "gs_subcollection" in tag_frame else "",
                spec.family or default[2],
                spec.aspect or default[3],
            )
        )
    return concat_memberships(frames)


def _memberships_cache_path(spec: MsigdbRemoteSource, release: dict[str, str], release_dir: Path) -> Path:
    """Normalized memberships are cached per archive checksum and every option that shapes them."""
    options = {
        "zip": release.get("zip_sha256") or release.get("zip_md5", ""),
        "summary_rds": release["summary_rds"],
        "db_species": spec.db_species.upper(),
        "collection": spec.collection or "",
        "subcollection": spec.subcollection or "",
        "family": spec.family or "",
        "aspect": spec.aspect or "",
    }
    label = "__".join([options["db_species"], options["collection"] or "all", options["subcollection"] or "all"])
    return derived_cache_path(release_dir / "alien_cache", f"memberships__{sanitize_term_id(label)}", options)


@register("msigdb_remote")
def load_msigdb_remote(spec: MsigdbRemoteSource, cache: ResourceCache) -> pd.DataFrame:
    release = release_info(spec)
    release_dir = ensure_release(spec, cache)
    force = cache.force or spec.force
    cache_path = _memberships_cache_path(spec, release, release_dir)
    if cache_path.exists() and not force:
        LOGGER.info("Loading cached normalized MSigDB memberships from %s", cache_path)
        frame = pd.read_parquet(cache_path).fillna("").astype(str)
        return concat_memberships([frame])[list(CANONICAL_COLUMNS)]
    summary = _filter_rows(_read_summary(release_dir, release, force), spec, species_column="db_target_species")
    if summary.empty:
        raise SourceError("No MSigDB remote collection files matched the requested filters.")
    frames = [
        normalize_remote_frame(read_rds_frame(release_dir / name), spec)
        for name in sorted(set(summary["df_rds"].astype(str)))
    ]
    memberships = concat_memberships(frames)
    with cache.file_lock(cache_path):
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        memberships[list(CANONICAL_COLUMNS)].to_parquet(cache_path, index=False)
    return memberships
