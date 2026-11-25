"""``build()`` and ``fetch()``: the ordered stages of an ALIEN build, nothing else."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, TypedDict, Unpack

import pandas as pd
import yaml

from .cache import LOCK_FILE_NAME, LockFile, ResourceCache
from .config import AlienConfig, ConfigInput, is_preset_name, preset_path, read_layer, validate_layers
from .dryrun import dry_run_plan
from .errors import ResourceError
from .filtering import FilterResult, filter_namespaces
from .io.atomic import atomic_write_json, atomic_write_text
from .io.gmt import write_gmt
from .io.tables import write_table
from .mapping import MappingResult, map_namespaces, problematic_source_ids, target_context
from .mapping.project import Memberships
from .model import Target
from .outputs import StagedOutput, prune_to_mode, write_manifest
from .reports import manifest, provenance, qc
from .resources.annotation import prepare_targets, write_target_annotation
from .resources.ensembl_archive import build_archive_resolver
from .resources.hgnc import build_hgnc_maps, load_hgnc_table
from .resources.ncbi import load_ncbi_maps
from .schema import (
    AMBIGUOUS_COLUMNS,
    ARCHIVE_AUDIT_COLUMNS,
    MAPPING_STATUS_COLUMNS,
    NCBI_AUDIT_COLUMNS,
    NON_GENE_COLUMNS,
    OUTPUT_SCHEMA_VERSION,
    UNMAPPED_COLUMNS,
)
from .sources import check_term_id_collisions, load_source_memberships
from .version import __version__

LOGGER = logging.getLogger("alien")


class BuildOptions(TypedDict, total=False):
    """Keyword options of :func:`build`; each overrides the matching config value."""

    outdir: str | Path | None
    source_dir: str | Path | None
    workers: int | None
    output_mode: str | None
    overrides: Sequence[ConfigInput] | None
    dry_run: bool
    force: bool
    force_download: bool
    offline: bool
    frozen: bool
    lock: str | Path | None


@dataclass(frozen=True, eq=False)
class BuildResult:
    outdir: Path
    namespaces: tuple[str, ...]
    n_source_memberships: int
    n_source_terms: int
    warnings: tuple[str, ...]
    plan: tuple[str, ...] = ()
    gmt_paths: dict[str, Path] = field(default_factory=dict)
    manifest_path: Path | None = None
    term_counts: dict[str, int] = field(default_factory=dict)
    mapping_summary: pd.DataFrame | None = None


def resolve_build_config(
    config: ConfigInput | AlienConfig | None = None, **options: Unpack[BuildOptions]
) -> AlienConfig:
    """Validate the layered config: config or preset, overrides, then explicit options."""
    layers: list[dict[str, Any]] = []
    if isinstance(config, AlienConfig):
        layers.append(config.model_dump(mode="json", exclude_none=True))
    elif isinstance(config, str) and is_preset_name(config):
        layers.append(read_layer(preset_path(config), resolve_relative_paths=False))
    elif config is not None:
        layers.append(read_layer(config))
    layers.extend(read_layer(override) for override in options.get("overrides") or [])
    explicit: dict[str, Any] = {}
    for key, section, name in (("source_dir", "project", "source_dir"), ("outdir", "project", "outdir")):
        if (value := options.get(key)) is not None:
            explicit.setdefault(section, {})[name] = str(value)
    if (workers := options.get("workers")) is not None:
        explicit.setdefault("runtime", {})["workers"] = workers
    if (mode := options.get("output_mode")) is not None:
        explicit.setdefault("outputs", {})["mode"] = str(mode).strip().lower()
    layers.append(explicit)
    return validate_layers(layers)


def effective_config_text(config: AlienConfig, *, keep_outdir: bool = False) -> str:
    """The effective config as YAML; project.outdir is omitted by default because the file lives inside it."""
    data = config.model_dump(mode="json", exclude_none=True)
    if not keep_outdir:
        data["project"].pop("outdir", None)
    return yaml.safe_dump(data, sort_keys=False)


def default_lock_path(config: ConfigInput | AlienConfig | None) -> Path:
    if isinstance(config, str | Path) and not (isinstance(config, str) and is_preset_name(config)):
        return Path(config).resolve().parent / LOCK_FILE_NAME
    return Path.cwd() / LOCK_FILE_NAME


def make_cache(config: AlienConfig, source: ConfigInput | AlienConfig | None, options: BuildOptions) -> ResourceCache:
    frozen = bool(options.get("frozen"))
    lock_option = options.get("lock")
    lock_path = Path(lock_option) if lock_option is not None else default_lock_path(source)
    if frozen and not lock_path.exists():
        raise ResourceError(f"--frozen needs a lock file; {lock_path} does not exist. Run `alien fetch` first.")
    lock = LockFile.load(lock_path) if (frozen or lock_option is not None) else LockFile(path=lock_path)
    return ResourceCache(
        root=config.project.source_dir,
        force=bool(options.get("force_download")),
        offline=bool(options.get("offline")),
        frozen=frozen,
        lock=lock,
    )


def build(config: ConfigInput | AlienConfig | None = None, **options: Unpack[BuildOptions]) -> BuildResult:
    """Build audited GMTs for ``config`` (a preset name, YAML path, dict or AlienConfig)."""
    model = resolve_build_config(config, **options)
    if options.get("dry_run"):
        plan = dry_run_plan(model)
        for line in plan:
            LOGGER.info("%s", line)
        return BuildResult(model.project.outdir, tuple(model.namespaces()), 0, 0, (), tuple(plan))
    cache = make_cache(model, config, options)
    staged = StagedOutput(model.project.outdir, force=bool(options.get("force")))
    try:
        result = _run(model, staged.staging, staged.failed_dir / "metadata", cache)
        prune_to_mode(staged.staging, model.outputs.mode)
        write_manifest(
            staged.staging, mode=model.outputs.mode, alien_version=__version__, schema_version=OUTPUT_SCHEMA_VERSION
        )
        staged.commit()
    except Exception:
        LOGGER.error("Build failed; diagnostics kept under %s", staged.fail())
        raise
    except BaseException:
        staged.discard()
        raise
    outdir = model.project.outdir
    gmt_paths = {name: outdir / "gmt" / f"{name}.gmt" for name in result.namespaces}
    return BuildResult(
        outdir=outdir,
        namespaces=result.namespaces,
        n_source_memberships=result.n_source_memberships,
        n_source_terms=result.n_source_terms,
        warnings=result.warnings,
        gmt_paths={name: path for name, path in gmt_paths.items() if path.exists()},
        manifest_path=outdir / "MANIFEST.json",
        term_counts=result.term_counts,
        mapping_summary=result.mapping_summary,
    )


def fetch(config: ConfigInput | AlienConfig | None = None, **options: Unpack[BuildOptions]) -> Path:
    """Download everything ``config`` needs, then write or update the lock file. Returns the lock path."""
    model = resolve_build_config(config, **options)
    cache = make_cache(model, config, {**options, "frozen": False, "offline": False})
    if cache.lock.path is not None and cache.lock.path.exists():
        cache.lock = LockFile.load(cache.lock.path)
    memberships = load_source_memberships(model, cache)
    hgnc = build_hgnc_maps(load_hgnc_table(cache, model.hgnc))
    ncbi = load_ncbi_maps(model.ncbi_gene, cache)
    targets = prepare_targets(model.targets, cache, [])
    resolver = build_archive_resolver(model.ensembl_archive, cache)
    if resolver is not None:
        members = Memberships(memberships)
        for target in targets:
            known = target_context(target, model, hgnc, ncbi, resolver).known_ids
            resolver.prefetch(problematic_source_ids(members, known), target.name)
    cache.lock.save()
    assert cache.lock.path is not None
    LOGGER.info("Wrote %s with %d pinned resources", cache.lock.path, len(cache.lock.resources))
    return cache.lock.path


def _report_tables(
    memberships: pd.DataFrame, targets: Sequence[Target], mapping: MappingResult, filtered: FilterResult
) -> dict[str, tuple[pd.DataFrame, Sequence[str] | None]]:
    """Every metadata/ and qc/ table, keyed by its path inside the output directory."""
    tables = {"unmapped": mapping.unmapped, "ambiguous": mapping.ambiguous, "non_gene": mapping.non_gene}
    term_manifest = manifest.build_term_manifest(
        memberships, mapping.terms, filtered.filtered, filtered.removed_size, filtered.removed_redundancy
    )
    return {
        "metadata/source_manifest.tsv": (manifest.build_source_manifest(memberships), None),
        "metadata/term_manifest.tsv.gz": (term_manifest, None),
        "metadata/target_metadata_fallbacks.tsv": (qc.metadata_fallback_table(targets), None),
        "metadata/target_output_genes.tsv.gz": (qc.target_gene_id_table(targets, "output_genes"), None),
        "metadata/target_gene_filter.tsv.gz": (qc.target_gene_id_table(targets, "gene_filter"), None),
        "metadata/removed_terms_size_filter.tsv.gz": (filtered.removed_size, None),
        "metadata/removed_terms_redundancy.tsv.gz": (filtered.removed_redundancy, None),
        "metadata/unmapped_genes.tsv.gz": (mapping.unmapped, UNMAPPED_COLUMNS),
        "metadata/dropped_non_gene_tokens.tsv.gz": (mapping.non_gene, NON_GENE_COLUMNS),
        "metadata/ambiguous_gene_mappings.tsv.gz": (mapping.ambiguous, AMBIGUOUS_COLUMNS),
        "metadata/ensembl_archive_audit.tsv.gz": (mapping.archive_audit, ARCHIVE_AUDIT_COLUMNS),
        "metadata/ncbi_gene_rescue_audit.tsv.gz": (mapping.ncbi_audit, NCBI_AUDIT_COLUMNS),
        "qc/ensembl_archive_summary.tsv": (qc.archive_summary(mapping.archive_audit), None),
        "qc/non_gene_token_summary.tsv": (qc.non_gene_summary(mapping.non_gene), None),
        "qc/ncbi_gene_rescue_summary.tsv": (qc.ncbi_summary(mapping.ncbi_audit), None),
        "qc/target_namespace_summary.tsv": (qc.target_namespace_summary(targets, memberships), None),
        "qc/mapping_status_summary.tsv": (mapping.status_summary, MAPPING_STATUS_COLUMNS),
        "qc/mapping_summary.tsv": (
            qc.mapping_summary(memberships, mapping.terms, tables, filtered.removed_size, mapping.status_summary),
            None,
        ),
        "qc/redundancy_summary.tsv": (filtered.redundancy_summary, None),
        "qc/collection_summary.tsv": (
            qc.collection_summary(mapping.terms, filtered.size_filtered, filtered.filtered),
            None,
        ),
    }


def _run(model: AlienConfig, outdir: Path, report_dir: Path, cache: ResourceCache) -> BuildResult:
    metadata = outdir / "metadata"
    warnings: list[str] = []
    effective = effective_config_text(model)
    atomic_write_text(metadata / "effective_config.yml", effective)
    LOGGER.info("Starting ALIEN build (%d worker process(es)); sources cached in %s", model.runtime.workers, cache.root)

    memberships = load_source_memberships(model, cache)
    LOGGER.info("Loaded %d memberships across %d terms", len(memberships), memberships["term_id"].nunique())
    check_term_id_collisions(memberships, model.term_id_collisions.action, metadata, report_dir, warnings)

    LOGGER.info("Loading mapping resources")
    hgnc_table = load_hgnc_table(cache, model.hgnc)
    targets = prepare_targets(model.targets, cache, warnings)
    for target in targets:
        write_target_annotation(target, metadata / f"gene_mapping_{target.name}.tsv.gz")
    write_table(hgnc_table, metadata / "gene_symbol_master.tsv.gz")

    LOGGER.info("Projecting gene sets to configured namespaces")
    resolver = build_archive_resolver(model.ensembl_archive, cache)
    hgnc, ncbi = build_hgnc_maps(hgnc_table), load_ncbi_maps(model.ncbi_gene, cache)
    mapping = map_namespaces(memberships, model, targets, hgnc=hgnc, ncbi=ncbi, resolver=resolver)
    qc.log_archive_audit(mapping.archive_audit, model.ensembl_archive.unresolved_print_limit)

    LOGGER.info("Applying size and redundancy filters")
    filtered = filter_namespaces(mapping.terms, model, model.runtime.workers)
    for namespace, terms in filtered.filtered.items():
        write_gmt(terms, outdir / "gmt" / f"{namespace}.gmt")
        if not terms:
            warnings.append(f"Empty output for namespace {namespace}.")

    LOGGER.info("Writing metadata and QC reports")
    reports = _report_tables(memberships, targets, mapping, filtered)
    for name, (frame, columns) in reports.items():
        write_table(frame, outdir / name, columns)
    atomic_write_json(
        outdir / "qc" / "collection_summary.json", reports["qc/collection_summary.tsv"][0].to_dict("records")
    )
    atomic_write_text(outdir / "qc" / "warnings.txt", "".join(f"{warning}\n" for warning in warnings))
    record = provenance.build_provenance(model, targets, cache, warnings, effective)
    atomic_write_json(metadata / "source_provenance.json", record)
    LOGGER.info("Done. Wrote ALIEN build outputs.")
    return BuildResult(
        outdir=outdir,
        namespaces=tuple(filtered.filtered),
        n_source_memberships=len(memberships),
        n_source_terms=int(memberships["term_id"].nunique()),
        warnings=tuple(warnings),
        term_counts={namespace: len(terms) for namespace, terms in filtered.filtered.items()},
        mapping_summary=reports["qc/mapping_summary.tsv"][0],
    )
