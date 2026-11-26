# Changelog

All notable changes to ALIEN are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses [Semantic Versioning](https://semver.org/).

## [0.2.0] - Unreleased

Hardening release: correct audit trail, reproducible and pinnable builds, typed configuration, and a modular code base. GMT contents are unchanged except where a fix below says so.

### Changed (outputs)

- Output schema version 2, recorded in `MANIFEST.json` and provenance.
- `metadata/term_manifest.tsv.gz` is long format: one row per (term, namespace) with `size_after_mapping`, `kept`, `removed_reason`, `representative_term_id` and `redundancy_cluster_id`. The per-row `created_at` column is gone.
- `removed_terms_*` tables: `mapped_size_before_output_genes` and `mapped_size_after_output_genes` (always equal) are replaced by `mapped_size`.
- `qc/mapping_summary.tsv` counts terms lost to size filters per source tag instead of repeating the family total.
- `metadata/source_provenance.json` is assembled from the run: ALIEN version, git describe, all package versions, platform, config sha256, every input file and remote resource with sha256, redundancy settings, worker count, cache mode. Hard-coded release and URL fields are gone.
- `metadata/effective_config.yml` lists every validated option and omits `project.outdir`.
- New `MANIFEST.json` in every output directory.
- `symbols.gmt` applies `manual_symbol_repairs` (previously only Ensembl namespaces did).
- Previous symbols shared by several genes are reported in `ambiguous_gene_mappings.tsv.gz` (stage `hgnc_previous_symbol`) instead of `unmapped_genes.tsv.gz`.
- Output-gene rows that are not Ensembl IDs keep an empty `ensembl_gene_id` in `target_output_genes.tsv.gz`.

### Fixed

- The term manifest attributed one namespace's removal reason to every namespace (5.1).
- Fourteen documented options had no effect: `exact_duplicate_removal`, `apply_per_family` and `use_hgnc_alias_mapping` now work; the rest are rejected with a message naming the replacement (5.2).
- Typos in configuration keys were silently ignored (5.3).
- Manual symbol repairs were skipped for the symbols namespace and could not override an alias (5.4).
- Ambiguous previous symbols were missing from the ambiguity audit (5.5).
- Provenance contained constants that did not describe the build (5.6).
- Relative paths resolved against the working directory instead of the config file (5.7).
- Failed `minimal`/`gmt` builds deleted the diagnostics the error message pointed to (5.8).
- Rebuilding left stale outputs and a crash could leave a half-written output directory (5.9).
- `--dry-run` validated nothing (5.10).
- Outputs were not byte-for-byte reproducible (gzip timestamps, per-row timestamps) (5.11).
- Duplicated and misleading QC columns (5.12).
- The Ensembl archive cache defaulted to `data/gmt_sources` (5.13).
- A custom MSigDB `release.zip_url` inherited the built-in archive name and checksum.
- Extracting a release archive could write outside its directory (zip-slip).

### Changed (resources)

- HGNC downloads come from HGNC's Google Cloud Storage bucket first; the EBI FTP path, which now returns 404, is a fallback.

### Added

- Typed configuration (pydantic v2) with `extends:`, `outputs.mode`, `redundancy.linkage` (`single` or `representative`), `hgnc.url` and `hgnc.archive_date`; errors name the YAML path; JSON Schema in `docs/config.schema.json`.
- `alien fetch`, `alien build --frozen/--offline/--lock`, and `alien.lock.json` pinning every remote resource and Enrichr library resolution.
- `alien validate`, `alien verify`, `alien config show|schema`, `alien --version`, `python -m alien`, `-v/-q/--log-file`, `--force`; exit codes 2/3/4/1 by error class.
- Plugin source types through the `alien.sources` entry-point group (`type: plugin`).
- `BuildResult.gmt_paths`, `manifest_path`, `term_counts` and `mapping_summary`; `alien.fetch()`.
- `py.typed`, `docs/config.md` and `docs/outputs.md` (generated from the code), descriptions for every option in the JSON Schema, `docs/architecture.md`, `CONTRIBUTING.md`, `CITATION.cff`.

### Performance

- Redundancy candidate pairs come from a prefix-filter similarity join: 50,000 terms in about 7 s and 0.5 GB (0.1.x: 10,000 terms took 99 s and 3.4 GB). Results are identical.
- Mapping resolves each distinct (symbol, source Ensembl ID) key once per target; the mapping process pool is gone.
- Parsed MSigDB collection files (`.rds`) are cached by content hash, so sources that share a collection (for example several C2 subcollections) parse it once; a first pathways fetch is several minutes faster.
- Parsed GTF gene tables are cached by content hash; NCBI maps and normalized MSigDB memberships are keyed by input hashes and options.

### Removed / deprecated

- `default_config.yml`, `alien.defaults` and the implicit default sources (MSigDB C2 and a `local_symbolic_libraries/*.gmt` glob); a config without sources is an error.
- The fourteen inert options listed above (now rejected) and `mapping.map_ensembl`/`map_symbols` (use `map_namespaces`).
- Deprecated, removed in 0.3: `alien.presets` (use `alien.config`), `alien.sources.build_hgnc_maps` and other moved helpers, `alien.mapping.resolve_current_symbol`, the target keys `annotation_gtf` and `gencode_version`, and the gene-list key spellings `file`, `column`, `gene_column`, `field`, `gene_symbol_column` and `genes`.
- `load_config()` and `load_preset_config()` return an `AlienConfig` instead of a dict.
- The library no longer configures the root logger; the CLI does.
- Dependency floors lowered to pandas 2.2, pyarrow 15, requests 2.31 and urllib3 2.0; new dependencies pydantic (2.7 or later) and filelock.
- The R analysis scripts and reports moved to `examples/analysis/` and are excluded from the sdist.

## [0.1.6]

- Bundled presets and output modes.

## [0.1.5]

- First PyPI release.
