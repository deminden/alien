# Usage

ALIEN builds combined GMT files from configured source libraries and projects them into configured target namespaces. A build can start from a bundled preset or from one YAML file.

## Run

```bash
alien build examples/pathways.yml --workers 4
alien build --config examples/pathways.yml --workers 4
```

The Python API is equivalent:

```python
from alien import build, load_preset_config

result = build("examples/pathways.yml", workers=4)
print(result.namespaces)
cfg = load_preset_config("cancer")
```

Use a bundled preset when you want ALIEN's packaged default for a common library family:

```bash
alien build pathways
alien build function
alien build disease
alien build cancer
```

The bundled preset names are `pathways`, `function_location`, `disease_phenotype`, and `cancer_dependency`. Short aliases are `function`, `disease`, and `cancer`. `alien presets` lists the installed presets and their package paths. Presets are packaged with ALIEN for offline config resolution and are GTEx/Gencode47-only by default.

Choose the GENCODE version your expression data was quantified on, not the newest release. GTEx v11 uses GENCODE v47, which is why the presets target it. A newer annotation adds gene IDs that can never appear in older data: they drop out of enrichment and shift the effective gene universe. If your counts use another release, override the target's `annotation.version` (and its `name`).

Use [examples/](../examples/) when you want the editable, study-specific YAMLs used by the included report scripts:

```bash
alien build examples/cancer_dependency.yml
alien build --config examples/cancer_dependency.yml
```

`--config` always treats its value as a filesystem path. Positional build input is resolved as a preset name/alias first and otherwise as a YAML path.

The bundled targets use full GENCODE v47 annotations, so they do not require a local GTEx expression/count matrix just to build the GMT. The example YAMLs intentionally keep hardcoded `data/gtex/...` output-gene paths for reproducing the included analyses.

TCGA-restricted targets need an output-gene list. The recount3 GTF is optional unless you configure it as a metadata fallback for measured IDs missing from the primary GENCODE annotation. Put these under `data/recount3/` or override their paths:

```text
data/recount3/tcga_gencode_v29_output_genes.tsv.gz
data/recount3/human.gene_sums.G029.gtf.gz
```

Overlay YAMLs are applied after the preset/config and before CLI overrides:

```bash
alien build pathways --override my-study.yml --source-dir scratch/alien_sources --outdir results/my_gmt
```

Mappings merge recursively. Lists replace, so overriding one target means supplying the full `targets` list you want to use:

```yaml
project:
  source_dir: scratch/alien_sources
  outdir: results/my_gmt

targets:
  - name: my_study_gencode47
    type: ensembl_gtf
    annotation:
      source: GENCODE
      version: "47"
    output_genes:
      path: data/my_study/expression.tsv.gz
      id_column: feature_id
```

For TCGA/recount3 cancer output, save this as `tcga-recount3.yml` and run `alien build cancer --override tcga-recount3.yml`:

```yaml
targets:
  - name: tcga_recount3_gencode29
    type: ensembl_gtf
    annotation:
      source: GENCODE
      version: "29"
      metadata_fallbacks:
        - source: recount3
          version: G029
          path: data/recount3/human.gene_sums.G029.gtf.gz
    output_genes:
      path: data/recount3/tcga_gencode_v29_output_genes.tsv.gz
      id_column: Ensembl_gene_ID
```

## Configuration Layers, Validation and Migration

A build's configuration is assembled from layers, lowest first:

1. Built-in defaults (defined once, on the typed configuration models).
2. The `extends:` chain of the config file. `extends:` names a bundled preset or another YAML file (relative to the declaring file); chains may be nested and loops are rejected.
3. The config file or preset itself.
4. Each `--override` file (or `overrides=[...]` entry) in order.
5. Explicit arguments: `--source-dir`, `--outdir`, `--workers`, `--output-mode`.

Mappings merge recursively and lists replace. Four keys extend a list instead of replacing it and are applied after merging: `gene_mapping.extra_non_gene_tokens`, `gene_mapping.extra_non_gene_token_regex`, `filtering.extra_broad_disease_regex`, and `source_priority.extra.<family>`.

Every option is validated. Unknown keys and out-of-range values fail with the YAML path of the offending entry, for example `targets[1].annotation.versoin: unknown option (check the spelling)`. A config must define at least one source. `alien config show [PRESET|FILE] [--override FILE ...]` prints the effective configuration; `alien config schema` prints the JSON Schema, which is also committed as [config.schema.json](config.schema.json) for editor completion (for example with the YAML language server: `# yaml-language-server: $schema=https://raw.githubusercontent.com/deminden/alien/main/docs/config.schema.json`).

`outputs.mode` (`full`, `minimal` or `gmt`) sets the output mode; `--output-mode` overrides it.

### Migration from 0.1.x

| Old key | New key | Status |
|---|---|---|
| `targets[].annotation_gtf` | `targets[].annotation.path` | deprecated, accepted with a warning; removed in 0.3 |
| `targets[].gencode_version` | `targets[].annotation.source: GENCODE` + `annotation.version` | deprecated, accepted with a warning; removed in 0.3 |
| `output_genes.file` / `gene_filter.file` | `.path` | deprecated; removed in 0.3 |
| `output_genes.column` / `.gene_column` / `.field` | `.id_column` | deprecated; removed in 0.3 |
| `output_genes.gene_symbol_column` | `.symbol_column` | deprecated; removed in 0.3 |
| `output_genes.genes` | `.ids` | deprecated; removed in 0.3 |
| `redundancy.apply_after_output_gene_projection` | none (always true) | removed, rejected |
| `redundancy.nested_overlap_filter` | none (never implemented) | removed, rejected |
| `filtering.remove_empty_after_output_gene_projection` | `filtering.min_size_default` | removed, rejected |
| `gene_mapping.strip_ensembl_version` | none (always stripped) | removed, rejected |
| `gene_mapping.allow_ambiguous_alias_mapping` | none (ambiguous aliases are always dropped and audited) | removed, rejected |
| `gene_mapping.prefer_gencode_exact_gene_name` | none (always preferred) | removed, rejected |
| `gene_mapping.prefer_current_hgnc_symbol` | none (always preferred) | removed, rejected |
| `gene_mapping.drop_unmapped_from_ensg_gmt` | none (always dropped) | removed, rejected |
| `ensembl_archive.unresolved_example_threshold` | `ensembl_archive.unresolved_print_limit` | removed, rejected |
| default MSigDB C2 source and `local_symbolic_libraries/*.gmt` glob | explicit `sources:` | removed; a config with no sources is an error |
| `--output-mode` only | `outputs.mode` | new; the flag still overrides |
| none | `extends:` | new |
| none | `redundancy.linkage` (`single` or `representative`) | new; default `single` keeps 0.1.x results |

`load_config()` and `load_preset_config()` now return a typed `AlienConfig` (attribute access, `model_dump()` for a dict) instead of a plain dict.

## Minimal Config

```yaml
project:
  # Download/cache root for source libraries and mapping resources.
  source_dir: data/alien_sources
  # Output root. ALIEN creates gmt/, metadata/, and qc/ below this folder.
  outdir: data/alien_gmt

outputs:
  include_symbols: true

sources:
  - type: msigdb_remote
    version: "2026.1"
    db_species: HS
    collection: C2

targets:
  - name: human_gencode47
    type: ensembl_gtf
    annotation:
      source: GENCODE
      version: "47"
```

This downloads/caches the configured MSigDB release archive and the GENCODE v47 annotation if they are missing, then writes:

```text
data/alien_gmt/
  gmt/
    symbols.gmt
    human_gencode47.gmt
  metadata/
  qc/
```

Those are the two normal filesystem locations in an ALIEN config. `project.source_dir` stores downloaded source GMTs, MSigDB/Enrichr caches, annotation files, HGNC/NCBI mapping resources, and other conversion inputs. `project.outdir` is the result root; ALIEN creates `gmt/`, `metadata/`, and `qc/` inside it. The CLI `--outdir` option is only an override for `project.outdir`.

Relative paths written in a YAML config file or `--override` file (`project.source_dir`, `project.outdir`, source `path`/`paths`/`cache_dir`, target `annotation.path`, `metadata_fallbacks[].path`, `output_genes` and `gene_filter` paths) resolve against the directory of the file that declares them, so a config works no matter where `alien` is started. Bundled presets are the exception: their `data/alien_sources` and `results/...` paths stay relative to the working directory. Paths passed on the command line or as Python dicts resolve against the working directory.

`--source-dir` is a convenience override for `project.source_dir`. Use it to redirect managed downloads and caches to a scratch or shared cache location without editing the YAML.

## Output Modes

`outputs.mode` in the config (or `--output-mode`) selects which files are kept:

- `full` (default): the complete `gmt/`, `metadata/` and `qc/` tree.
- `minimal`: the GMTs plus the compact reproducibility record (effective config, provenance, source and term manifests, target gene lists, mapping, redundancy and target summaries, warnings). The exact list is in [outputs.md](outputs.md).
- `gmt`: only `gmt/*.gmt`. Use it only when the reproducibility record is preserved elsewhere.

Every mode builds in a hidden staging directory next to `project.outdir`, prunes it to the mode, writes `MANIFEST.json` (every file with sha256 and size) last, and then swaps it into place with an atomic rename, so a crash never leaves a half-updated output directory and outputs of removed targets never linger. If the build fails, the previous output is left untouched and the staging tree, including any `metadata/term_id_collisions.tsv`, is kept under `<outdir>/failed/`; error messages point there. ALIEN refuses to replace a non-empty directory that is not an earlier ALIEN build unless you pass `--force` (`force=True` in Python).

`alien verify OUTDIR` re-hashes an output directory against its `MANIFEST.json` and reports missing, changed and unexpected files.

Outputs are reproducible byte for byte: tables are written with a fixed column order and gzip members carry no timestamp. The only time-dependent value is `created_at` in `metadata/source_provenance.json`, which honours `SOURCE_DATE_EPOCH`.

## Dry Runs, Lock Files and Offline Builds

`alien build CONFIG --dry-run` validates the configuration, checks that every local input exists, lists which remote resources are cached and which would be downloaded, and prints the planned outputs. It writes nothing.

HGNC, NCBI Gene, GENCODE, MSigDB and Enrichr are remote resources that change over time. To pin them:

```bash
alien fetch my-study.yml            # download everything the config needs; write my-study's alien.lock.json
alien build my-study.yml --frozen   # no network; every cached file must match the lock byte for byte
alien build my-study.yml --offline  # no network; use whatever is cached, lock optional
```

`alien.lock.json` is written next to the config file (in the working directory for presets and Python dicts); `--lock PATH` chooses another location. For every remote resource it records the URL, sha256, size, retrieval time, ETag and Last-Modified where available, and the upstream version (MSigDB release, GENCODE version, HGNC archive date, Enrichr library name). Enrichr libraries selected by `match:` are pinned by name, so a frozen build keeps using the library chosen when the lock was written even if a newer one appears. A frozen build fails with a clear message when a cached file is missing or differs from the lock. Ensembl archive lookups are pinned too: `alien fetch` performs them and records the lookup cache.

To pin HGNC to a monthly archive instead of the current release, set `hgnc.archive_date` (for example `"2025-01-01"`), or point `hgnc.url` at any complete-set file.

Caches are safe to share between concurrent builds (Snakemake, Nextflow, SLURM arrays): every download and derived cache is written through a per-file lock and an atomic rename. Derived caches (normalized MSigDB memberships, NCBI rescue maps, parsed GTF gene tables) are keyed by the sha256 of their inputs, the options that shape them and a cache schema version, so changing an option such as `ncbi_gene.use_gene_info_synonyms` never reuses a stale cache. Ensembl archive "not found" answers are stamped with the Ensembl release that gave them and are re-checked after a new release.

Output modes only stage `project.outdir`; `project.source_dir` is always the shared cache. Point it at scratch with `--source-dir` to make downloads temporary too.

## Source Types

`msigdb_remote` reads a managed `msigdbr` release archive from Zenodo and keeps normalized ALIEN Parquet caches for repeated builds:

```yaml
sources:
  - type: msigdb_remote
    version: "2026.1"
    db_species: HS
    collection: C5
    subcollection: GO:MF
```

`msigdb_cache` reads local TSV files named like `msigdbr_REACTOME.tsv.gz`:

```yaml
sources:
  - type: msigdb_cache
    path: data/alien_sources/msigdb
    include_c4_cm: false
```

`enrichr_remote` downloads and caches Enrichr libraries by public library name:

```yaml
sources:
  - type: enrichr_remote
    libraries:
      - name: KEGG_2021_Human
        source_tag: KEGG
        family: biology_process_pathway
        aspect: pathway
      - name: ClinVar_2019
        source_tag: CLINVAR
        family: disease_phenotype
        aspect: disease
```

For versioned library names, exact names are tried first. If an exact library name is not present, provide `match` to select the latest matching library deterministically:

```yaml
sources:
  - type: enrichr_remote
    libraries:
      - name: ClinVar
        match: "^ClinVar_[0-9]{4}$"
        source_tag: CLINVAR
        family: disease_phenotype
        aspect: disease
```

Downloaded Enrichr metadata is cached at `data/alien_sources/enrichr/datasetStatistics.json`; GMT files are cached as `data/alien_sources/enrichr/<library>.gmt`. The resolved library name, match method, and regex candidate names are written to `metadata/source_manifest.tsv`.

For MSigDB remote sources, ALIEN keeps two cache layers under `project.source_dir`: the verified upstream release archive/extracted files and a derived `alien_cache/` directory with normalized Parquet memberships. The first build for a source filter performs conversion from the upstream RDS tables; later builds read the Parquet cache directly.

Cached remote resources are reused by default. To refresh MSigDB archives, Enrichr metadata/GMTs, GENCODE annotations, HGNC/NCBI resources, and Ensembl archive caches, run:

```bash
alien build --config examples/pathways.yml --force-download
```

For a single remote source, set `force: true` on that source:

```yaml
sources:
  - type: enrichr_remote
    force: true
    libraries:
      - name: ClinVar
        match: "^ClinVar_[0-9]{4}$"
        source_tag: CLINVAR
        family: disease_phenotype
        aspect: disease
```

This matters for regex-matched Enrichr libraries: ALIEN chooses the latest matching library from the cached `datasetStatistics` file unless the cache is refreshed. For publication or protocol configs, exact Enrichr library names are preferable; when regex matching is used, keep `metadata/source_manifest.tsv` with the released outputs.

Every enabled source is required. If a configured remote library cannot be resolved/downloaded, or a configured cache file is missing, the build fails rather than silently weakening the output. Use `enabled: false` to keep a source in a config file without including it in a particular build.

```yaml
sources:
  - type: enrichr_remote
    enabled: false
    libraries:
      - name: ClinVar_2019
```

NCBI rescue maps are prepared once into a JSON cache after the raw NCBI files are read. Repeated builds also reuse in-memory mapping lookups during each run, and multiple configured target namespaces can be mapped in parallel with `--workers`.

`msigdb_tsv` reads one or more MSigDB-like tables with term names and symbols:

```yaml
sources:
  - type: msigdb_tsv
    path: data/custom/reactome.tsv
    source_tag: REACTOME
    collection: C2
    family: biology_process_pathway
    aspect: pathway
```

`symbol_gmt` and `enrichr_gmt` read GMT files whose members are gene symbols:

```yaml
sources:
  - type: symbol_gmt
    paths:
      - data/custom/library.gmt
    source: Local library
    source_tag: LOCAL
    family: biology_process_pathway
    aspect: pathway
```

`canonical_tsv` reads ALIEN's normalized table format. Required columns are `term_id` and `gene_symbol`; additional metadata columns include `source`, `source_tag`, `collection`, `family`, `aspect`, `gene_id`, and `gene_id_namespace`.

## Term ID Collisions

`term_id` values are the stable identifiers used for merging, filtering, auditing, and GMT writing. Multiple gene rows for the same term are expected, but the same `term_id` must not describe two different source terms.

By default ALIEN stops if a `term_id` has conflicting source or term metadata:

```yaml
term_id_collisions:
  action: error
```

When collisions are found with `action: error`, the build stops and the audit table is kept at `<outdir>/failed/metadata/term_id_collisions.tsv` (the error message prints the path); with `warn` or `merge` it is written to `metadata/term_id_collisions.tsv` in the output. The table lists each conflicting identity, source metadata, row counts, and sample genes. The usual fix is to rename or prefix one source's term IDs.

For legacy inputs where intentional merging is acceptable, the guard can be relaxed:

```yaml
term_id_collisions:
  action: merge
```

`merge` (and `warn`) keep the previous behavior: memberships with the same `term_id` are combined and, in each namespace, the metadata of the term's first mapped membership is used downstream. ALIEN still writes the collision audit table and adds a warning.

## Targets

The implemented target adapter is `ensembl_gtf`. It projects source terms into Ensembl stable gene IDs using a target annotation GTF. The GTF must contain `gene` records with `gene_id` and `gene_name` attributes.

Each target has two separate roles:

- output gene set: the Ensembl IDs allowed in the final GMT;
- annotation helper: the GTF used to map source gene symbols to Ensembl IDs and supply audit metadata.

By default, the annotation GTF supplies both. When `output_genes` is configured, the output gene file/list supplies the final target ID set and the GTF is only a mapping/metadata helper. The `id_column` must contain Ensembl gene IDs; `symbol_column` is optional metadata, not a symbol-only target namespace. When `gene_filter` is configured, the final target ID set is the intersection of annotation IDs and the configured filter IDs.

The primary annotation is always the first mapping authority. If a dataset-specific output namespace contains Ensembl IDs absent from that primary annotation, `annotation.metadata_fallbacks` can add metadata from another GTF only for those missing output genes. Fallback IDs participate in mapping/auditing, but duplicate symbols already present in the primary annotation keep the primary mapping.

For human GENCODE releases, `source` and any numeric `version` are enough. ALIEN builds the official FTP URL as `release_<version>/gencode.v<version>.annotation.gtf.gz`:

```yaml
targets:
  - name: human_gencode47
    type: ensembl_gtf
    annotation:
      source: GENCODE
      version: "47"
```

For another Ensembl-style GTF origin, provide a local path:

```yaml
targets:
  - name: my_annotation
    type: ensembl_gtf
    annotation:
      source: MyProvider
      version: "2026-05"
      path: data/annotations/my_provider.gtf.gz
```

or a URL:

```yaml
targets:
  - name: my_remote_annotation
    type: ensembl_gtf
    annotation:
      source: MyProvider
      version: "2026-05"
      url: https://example.org/annotations/my_provider.gtf.gz
```

ALIEN caches URL-based target annotations under `data/alien_sources/targets/<target_name>/`.

If a GTF uses different attribute names, map them with `annotation.attributes`:

```yaml
targets:
  - name: custom_annotation
    type: ensembl_gtf
    annotation:
      source: MyProvider
      path: data/annotations/custom.gtf.gz
      attributes:
        gene_id: ID
        gene_name: Name
        gene_biotype: biotype
```

These keys refer to attributes in the ninth GTF field, not the fixed positional GTF columns.

Non-Ensembl target ID systems, such as Entrez, UniProt, or RefSeq output GMTs, are not supported by this adapter yet.

`output_genes` defines the Ensembl IDs that belong to a target output namespace when the downstream dataset has its own measured-gene output gene set. Normal whole-annotation builds can omit it.

Target ID set rules:

- no `output_genes` and no `gene_filter`: use all IDs from the annotation GTF.
- `output_genes`: use IDs from the configured output gene file/list.
- `gene_filter`: use IDs present in both the annotation GTF and the configured filter file/list.

```yaml
targets:
  - name: human_gencode47_study_only
    type: ensembl_gtf
    annotation:
      source: GENCODE
      version: "47"
    output_genes: data/study/count_matrix.tsv
```

If the matrix is tabular and the output gene ID column is known, use the explicit form:

```yaml
targets:
  - name: human_gencode47_study_only
    type: ensembl_gtf
    annotation:
      source: GENCODE
      version: "47"
    output_genes:
      path: data/study/expression.tsv.gz
      id_column: feature_id
```

Without `id_column`, ALIEN tries common gene-ID column names and then falls back to the first tabular column. For Parquet files, `id_column: index` reads the row index.

If the output-gene file also contains gene symbols, provide `symbol_column`. ALIEN uses those symbols as an additional mapping helper for output genes absent from the annotation GTF, and records the metadata source in `metadata/gene_mapping_<target>.tsv.gz`:

```yaml
targets:
  - name: study_expression_ids
    type: ensembl_gtf
    annotation:
      source: GENCODE
      version: "47"
    output_genes:
      path: data/study/expression.tsv.gz
      id_column: Ensembl_gene_ID
      symbol_column: gene_symbol
```

For dataset-specific annotations, add a metadata fallback instead of replacing the primary annotation. This is useful when the final output IDs come from a measured matrix, but the preferred biological mapping helper is still an official release such as GENCODE:

```yaml
targets:
  - name: tcga_recount3_gencode29
    type: ensembl_gtf
    annotation:
      source: GENCODE
      version: "29"
      metadata_fallbacks:
        - source: recount3
          version: G029
          path: data/recount3/human.gene_sums.G029.gtf.gz
    output_genes:
      path: data/recount3/tcga_gencode_v29_output_genes.tsv.gz
      id_column: Ensembl_gene_ID
```

The default fallback mode is deliberately conservative: ALIEN loads the fallback GTF, selects only output IDs missing from the primary annotation, and records the selected rows as `metadata_fallback` in `metadata/gene_mapping_<target>.tsv.gz`. The audit table `metadata/target_metadata_fallbacks.tsv` records how many output IDs were missing before and after each fallback.

`gene_filter` is the stricter alternative for true restriction behavior. It keeps only IDs present in both the annotation and the configured file/list:

```yaml
targets:
  - name: human_gencode47_filtered
    type: ensembl_gtf
    annotation:
      source: GENCODE
      version: "47"
    gene_filter:
      path: data/study/allowed_genes.tsv
      id_column: ensembl_gene_id
```

Use only one of `output_genes` or `gene_filter` for a target. ALIEN rejects configs that define both because they mean different set operations.

The same shape works from Python because `alien.build()` accepts a config dictionary:

```python
from alien import build

cfg = {
    "sources": [{"type": "symbol_gmt", "path": "data/source.gmt"}],
    "targets": [
        {
            "name": "human_gencode47_study_only",
            "type": "ensembl_gtf",
            "annotation": {"source": "GENCODE", "version": "47"},
            "output_genes": {
                "path": "data/study/expression.tsv.gz",
                "id_column": "feature_id",
                "symbol_column": "gene_symbol",
            },
        }
    ],
}

result = build(cfg, outdir="data/alien_gmt", workers=4)
```

For code-only workflows, `output_genes` can also carry inline Ensembl IDs:

```python
"output_genes": {"genes": ["ENSG00000141510.18", "ENSG000002"]}
```

## Filtering

Size filters remove terms outside configured minimum and maximum mapped sizes:

```yaml
filtering:
  min_size_default: 10
  max_size_default: 500
  max_size_disease: 500
  max_size_cancer: 500
  drop_broad_disease_terms: true
```

`broad_disease_regex` and `extra_broad_disease_regex` control broad disease-term removal.

## Redundancy

Redundancy filtering is configurable. ALIEN first collapses exact duplicate memberships, then clusters highly overlapping terms by Jaccard similarity within each namespace and family:

```yaml
redundancy:
  exact_duplicate_removal: true
  jaccard_cutoff: 0.85
  apply_per_family: true
```

`exact_duplicate_removal: false` skips the exact-membership pass, and `apply_per_family: false` compares terms across all families in a namespace (the per-family `source_priority` lists are concatenated in family order). Terms with Jaccard similarity strictly greater than `jaccard_cutoff` are treated as redundant. Representatives are chosen by `source_priority`, then larger mapped size, then names that do not start with a generic prefix ("regulation of", "process", ...), then shorter name, then term ID.

`redundancy.linkage` controls how similar pairs become clusters:

- `single` (default, the 0.1.x behaviour): similar terms are joined transitively. If A~B and B~C, all three form one cluster and one representative is kept, even when A and C are not similar to each other.
- `representative`: terms are visited in representative order; a term is removed only when it is above the cutoff against a term that was already kept, which is recorded as its `representative_term_id`. No chaining, so more terms survive.

The choice is recorded in provenance. Candidate pairs are found with a prefix-filter similarity join, so large families (tens of thousands of terms) filter in seconds.

Source priority is configured by family:

```yaml
source_priority:
  biology_process_pathway: [REACTOME, WIKIPATHWAYS, KEGG_MEDICUS, GOBP]
  biology_function_location: [GOMF, GOCC]
```

## Mapping

ALIEN maps source symbols through audited human mapping resources:

- HGNC current symbols.
- HGNC previous and alias symbols.
- NCBI Gene history/info rescue, enabled by default for configured legacy source tags.
- Ensembl archive lookup, enabled by default for source Ensembl IDs absent from the target annotation, any configured target output genes/filter, and HGNC's current Ensembl IDs.

For each membership and target the decision order is: configured non-gene token; source Ensembl ID already in the target; Ensembl archive replacement (for retired source IDs only); HGNC symbol resolution; annotation candidates for the resolved symbol (preferring output genes, then the primary annotation over fallbacks and output-gene symbol rows, then a unique protein-coding gene); the official HGNC Ensembl ID as a tie-breaker; NCBI Gene rescue (only for `ncbi_gene.allowed_source_tags`); otherwise the membership is classified as ambiguous, unmapped or not in the output genes.

HGNC symbol resolution tries, in order: an exact current symbol, `gene_mapping.manual_symbol_repairs` (a repair never overrides an exact current symbol, but does override previous-symbol and alias matches), an unambiguous previous symbol, an ambiguous previous symbol, an ambiguous alias, an unambiguous alias (`use_hgnc_alias_mapping: false` skips both alias steps). Repairs apply to the `symbols` namespace too.

Ambiguous mappings (previous symbols and aliases shared by several genes included) are dropped rather than guessed and are written to `metadata/ambiguous_gene_mappings.tsv.gz`.

Disable these fallback layers explicitly when a build must avoid managed remote mapping resources:

```yaml
ncbi_gene:
  enabled: false

ensembl_archive:
  enabled: false
```

## Outputs

Every file and column is listed in [outputs.md](outputs.md), which is generated from the code. The most useful ones:

- `gmt/<namespace>.gmt`: one GMT per namespace (`symbols` plus each target).
- `metadata/source_manifest.tsv`: one row per resolved source collection with source tag, collection, URL and license note, term and member counts, and Enrichr resolution fields.
- `metadata/term_manifest.tsv.gz`: one row per term and namespace with `size_after_mapping`, `kept`, `removed_reason` (`too_small`, `too_large`, `broad_disease_term`, `exact_duplicate`, `jaccard_duplicate` or `no_mapped_genes`), `representative_term_id` and `redundancy_cluster_id`. A term can be kept in one namespace and removed in another.
- `metadata/source_provenance.json`: what the build read and how (versions, hashes, resources, settings).
- `metadata/unmapped_genes.tsv.gz`, `ambiguous_gene_mappings.tsv.gz`, `dropped_non_gene_tokens.tsv.gz`, `ensembl_archive_audit.tsv.gz`, `ncbi_gene_rescue_audit.tsv.gz`: the mapping audit.
- `qc/mapping_summary.tsv`, `qc/mapping_status_summary.tsv`, `qc/redundancy_summary.tsv`, `qc/target_namespace_summary.tsv`, `qc/warnings.txt`.

`metadata/target_output_genes.tsv.gz` audits configured `output_genes`, including the input ID, stripped Ensembl ID, optional dataset-provided symbol, ID type, and whether the ID had annotation-helper metadata after primary annotation and metadata fallback rows were applied. Rows whose ID is not an Ensembl gene ID are kept with an empty `ensembl_gene_id` and counted in a warning. `metadata/target_metadata_fallbacks.tsv` audits fallback use per target. `qc/target_namespace_summary.tsv` summarizes `output_gene_source`, primary/fallback annotation helper sizes, annotation metadata coverage, and the number of target IDs missing annotation metadata.

## Command-Line Reference

```text
alien build CONFIG_OR_PRESET [--override FILE]... [--source-dir DIR] [--outdir DIR] [--workers N]
            [--output-mode full|minimal|gmt] [--dry-run] [--force] [--force-download]
            [--offline | --frozen] [--lock FILE]
alien fetch CONFIG_OR_PRESET [--override FILE]... [--source-dir DIR] [--lock FILE] [--force-download]
alien validate CONFIG_OR_PRESET [--override FILE]...
alien verify OUTDIR
alien config show [CONFIG_OR_PRESET] [--override FILE]...
alien config schema
alien presets
alien --version
```

Global options: `-v/--verbose` (debug logging and full tracebacks), `-q/--quiet` (warnings and errors only), `--log-file FILE`. `python -m alien` is equivalent to `alien`.

Exit codes: `0` success, `2` invalid configuration, `3` input data problem (missing or malformed source files, term-ID collisions, `alien verify` mismatches), `4` network or cache problem (download failed, offline resource missing, lock mismatch), `1` anything unexpected. Errors print one line; add `--verbose` for the traceback.

## Python API

```python
from alien import build, fetch, load_config

config = load_config("my-study.yml")  # a validated AlienConfig
fetch("my-study.yml")  # writes alien.lock.json
result = build("my-study.yml", frozen=True, workers=8, output_mode="minimal")
result.gmt_paths["gtex_v11_gencode47"]  # Path to the GMT
result.term_counts  # terms per namespace
result.mapping_summary  # qc/mapping_summary.tsv as a DataFrame
```

`build()` accepts a preset name, a YAML path, a dict or an `AlienConfig`, plus keyword options `outdir`, `source_dir`, `workers`, `output_mode`, `overrides`, `dry_run`, `force`, `force_download`, `offline`, `frozen` and `lock`. The library logs to the `alien` logger and never configures logging itself.

## Custom Source Types

A project can add its own source type without forking ALIEN. Publish a function `load(options: dict, cache) -> pandas.DataFrame` that returns the canonical membership columns (at least `term_id` and `gene_symbol`) under the `alien.sources` entry-point group:

```toml
[project.entry-points."alien.sources"]
internal_db = "my_package.alien_source:load"
```

and configure it:

```yaml
sources:
  - type: plugin
    plugin: internal_db
    options: {release: "2026-03"}
```
