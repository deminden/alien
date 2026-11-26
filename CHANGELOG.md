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

### Fixed (output review)

- GMT and Enrichr term names were sanitized before they were stored, so `original_name` and `display_name` held the term ID instead of the source spelling, two names that differ only in punctuation or non-ASCII letters were merged without a collision report, and multi-word `broad_disease_regex` patterns never matched. Names are now kept as written; names that sanitize to one `term_id` go through `term_id_collisions` (kept apart by renaming by default, see below). A name listed twice in one GMT is logged.
- Symbols that differ from HGNC only in letter case (`C1ORF112` for `C1orf112`) were unmapped. They are now resolved after every exact match fails and reported with the new mapping status `hgnc_case_insensitive_symbol`; ties are written to `ambiguous_gene_mappings.tsv.gz`. Set `gene_mapping.case_insensitive_symbols: false` for the old behaviour.
- GMT members with a numeric weight (`TP53,1.0`) were unmapped; the weight is now dropped.
- `--frozen` failed after any other build added lookups to the shared Ensembl archive cache. The lock now pins, per target, a digest of the lookups that target reads (`ensembl_archive:lookups:<target>`); lock files written before this change need one `alien fetch`.
- Refusing to overwrite a directory that is not an ALIEN build exits with code 2 instead of reporting an unexpected error.
- Derived Parquet caches (GTF gene tables, MSigDB summaries and memberships) are written atomically.
- A `;` inside a GMT description value is written as `,` so the `key=value;...` format stays parseable.

### Changed (gene-set lists)

From the 2026-11-26 review of real builds (GENCODE v47, all presets).

- **Presets.** `cancer_dependency` is replaced by `cancer_programs` (MSigDB C6, C4 3CA, C9; aliases `cancer` and
  `cancer_dependency`) and the opt-in `cell_line_signatures` (DepMap, CCLE, NCI-60; alias `cell_lines`). The old preset
  was 70% per-cell-line terms (85% with C4:CGN); in the TCGA lung example 72% of its terms were significant and the
  top hits were haematopoietic cell lines (immune infiltration) and GNF2 proliferation neighbourhoods. C4:CGN is no
  longer bundled. `pathways` now starts with MSigDB Hallmark, which also leads its source priority.
- **Collisions.** `term_id_collisions.action` defaults to `rename`: when two source names sanitize to one term ID, the
  first keeps it and the other gets a stable `__<hash>` suffix; `term_id_collisions.tsv` gains `resolved_term_id`. With
  the previous `error` default the bundled `disease_phenotype` preset failed (DisGeNET `beta Thalassemia` and
  `beta^+^ Thalassemia`).
- **Term exclusions.** New `filtering.exclude_terms`: named rules remove terms by name before size filtering, with
  reason `excluded_term` and the rule in the new `exclusion_rule` column of `removed_terms_size_filter`. Bundled rules:
  KEGG MEDICUS `VARIANT`/`PATHOGEN` pathways (254 terms), pathways named after one pathogen (Reactome 55,
  WikiPathways 41), HPO onset/inheritance/course modifiers (45), DisGeNET animal-model terms (22).
- **Broad disease terms.** The default `broad_disease_regex` also removes top-level Disease Ontology classes
  ("immune system disease", "disease of cellular proliferation"), generic processes ("inflammation", "neoplasm
  metastasis") and top-level HPO organ-system categories; names are also matched with underscores read as spaces, so
  the patterns reach MSigDB HPO names. The `disease_phenotype` preset no longer overrides the list (95 terms removed
  instead of 13).
- **Short aliases.** For text-mined sources (`gene_mapping.short_alias_source_tags`, default DISGENET and Jensen
  text-mining) alias-only matches shorter than `gene_mapping.min_alias_length` (4) are rejected with reason
  `short_alias_rejected`: `HCC` in liver-cancer terms had become HYCC1, `CCT` FLVCR2, `NT3` SORT1. The disease preset
  keeps the two checked exceptions `HRX` and `ST2` as manual repairs.
- **Legacy mitochondrial names.** `ND1`..`ND6`, `ND4L`, `CYTB`, `ATP6`, `ATP8`, `COX3`, `RNR1/2` and the tRNA codes
  (`TRNK`, `TRNL1`, ...) map to their `MT-` genes (status `mitochondrial_legacy_symbol`, option
  `gene_mapping.mitochondrial_legacy_symbols`). `TRNP` was mapped to the unrelated alias target TRNP1 and `ND1` was
  ambiguous; `COX1`/`COX2` stay ambiguous (PTGS1/PTGS2).
- **Spreadsheet dates.** `1-Sep`, `MAR-6` and similar are repaired to SEPTIN/MARCHF genes (status `excel_date_symbol`,
  option `gene_mapping.excel_date_symbols`); `1-Mar`/`2-Mar` stay ambiguous with MTARC1/MTARC2.
- **`LOC<GeneID>` tokens** are resolved through the NCBI Gene ID they encode, following replaced IDs (status
  `ncbi_gene_id`, for every source because it is an ID lookup). CCLE mapping rises from 91.3% to 96.7%.
- **Audit.** New `metadata/symbol_rescues.tsv.gz`: each distinct token mapped by anything other than a source Ensembl
  ID or an exact current symbol, with its status, symbol, Ensembl ID and how many memberships and terms it reached.
- C4 3CA/CGN/CM and Hallmark sets loaded with `msigdb_remote` get their family and aspect without an explicit
  `family:`; the default `source_priority` lists the tags the remote loader produces (`3CA`, `CGN`, `CM`, `H`).
- Output schema version 3.

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
