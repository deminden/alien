<!-- Generated from README.md by scripts/generate_pypi_readme.py. Do not edit directly. -->

# ALIEN

**ALIEN: Audited Library Integration for External Namespaces** is a fully Python-based tool for building namespace-specific GMT libraries for human gene-set workflows.

ALIEN is centered on one job: take configured source libraries, normalize their gene memberships into a canonical table, project them into configured target namespaces, and write combined GMT files with audit metadata.

## Install

Install ALIEN from PyPI:

```bash
pip install bioalien
```

The PyPI distribution name is `bioalien`; the Python import and command-line tool are `alien`.

## Quick Start

Run ALIEN with a YAML config:

```bash
alien build config.yml --workers 16
```

The config controls which source libraries are used, which target namespace is written, where downloaded resources are cached (`project.source_dir`), and where outputs are written (`project.outdir`). The main outputs are:

- `gmt/<target_namespace>.gmt`: combined GMT for each configured target namespace.
- `metadata/`: source manifests, term manifests, gene mapping audits, filtering logs, and provenance.
- `qc/`: collection, mapping, redundancy, target coverage, and warning summaries.

## Presets

For a quick GENCODE v47 build without writing a YAML first, use a bundled preset:

```bash
alien build pathways --workers 16
alien build cancer --output-mode minimal  # keep only GMTs plus minimal reproducibility metadata
```

Available preset aliases:

- `pathways`: Reactome, WikiPathways, KEGG MEDICUS, and GO biological process terms.
- `function`: GO molecular function and cellular component terms.
- `disease`: HPO, DisGeNET, ClinVar, GWAS Catalog, and Jensen disease libraries.
- `cancer`: cancer and dependency signatures.

By default they write against the full GENCODE v47 (matching GTEx).

Override any part of a preset with an overlay YAML:

```bash
alien build cancer --override tcga-recount3.yml
```

Overlay mappings deep-merge; lists such as `targets` and `sources` replace the previous list. For TCGA-restricted output, the key input is an output-gene list such as `data/recount3/tcga_gencode_v29_output_genes.tsv.gz` and `data/recount3/human.gene_sums.G029.gtf.gz` as a metadata fallback.

Examples with hardcoded study paths are in [examples/](https://github.com/deminden/alien/tree/v0.1.6/examples/). Report walkthroughs are in [docs/sex_contrast/gtex_thyroid/analysis.md](https://github.com/deminden/alien/blob/v0.1.6/docs/sex_contrast/gtex_thyroid/analysis.md) and [docs/sex_contrast/tcga_lung/analysis.md](https://github.com/deminden/alien/blob/v0.1.6/docs/sex_contrast/tcga_lung/analysis.md).

See [docs/usage.md](https://github.com/deminden/alien/blob/v0.1.6/docs/usage.md) for overlay examples, output modes, and the full configuration reference.

## How It Works

ALIEN builds one canonical membership table from configured sources, audits source gene symbols against human mapping resources, then projects each term into the requested output namespaces.

The target namespace is defined in the config. For an Ensembl-style namespace, provide a target name and a GTF annotation:

```yaml
targets:
  - name: human_gencode49
    type: ensembl_gtf
    annotation:
      source: GENCODE
      version: "49"
```

This writes `gmt/human_gencode49.gmt`.

For human GENCODE releases, any numeric `version` is enough; ALIEN builds the official FTP URL and caches the GTF under `data/alien_sources/gencode/` if it is missing. You can still provide `annotation.path` to pin a local file explicitly.

`source: GENCODE` is not special to the config shape; other Ensembl-style GTF origins can use the same adapter by providing a local path or URL. Fully different target ID systems, such as Entrez or UniProt GMT output, are planned as future target adapters.

Targets separate the output gene set from the annotation helper. By default the GTF supplies both. Use `output_genes` when a dataset file supplies the final Ensembl ID namespace, or `gene_filter` when the annotation namespace should be intersected with a file/list of allowed IDs. In `output_genes` builds, `id_column` contains Ensembl IDs and the optional `symbol_column` adds dataset-provided symbol metadata for IDs absent from the GTF. `annotation.metadata_fallbacks` can add secondary GTF metadata only for output genes missing from the primary annotation; primary annotation mappings keep priority.

The build keeps audit outputs beside the GMTs so each namespace projection can be traced back to source terms, symbol repairs, unmapped genes, filters, redundancy decisions, and provenance. The compact `metadata/source_manifest.tsv` table is the main record of the exact source collections used in a build; for regex-matched Enrichr libraries it records the resolved library, match method, and candidate names.

`term_id` values must identify one source term unambiguously. If two source libraries reuse the same `term_id` for different term metadata, ALIEN fails by default and writes `metadata/term_id_collisions.tsv` so the IDs can be renamed or prefixed before rebuilding.

Redundancy filtering removes exact duplicate terms and then clusters highly overlapping terms by Jaccard similarity within each namespace and family. The default cutoff is `0.85`:

```yaml
redundancy:
  jaccard_cutoff: 0.85
```

Within each redundant cluster, ALIEN keeps one representative using the configured source priority, then term size and name-based tie-breaks.

Source priority is configured per term family and controls which library wins when redundant terms overlap. For example, pathway terms can prefer `REACTOME` over broader ontology-derived terms:

```yaml
source_priority:
  biology_process_pathway: [REACTOME, WIKIPATHWAYS, KEGG_MEDICUS, GOBP]
```

See [docs/usage.md](https://github.com/deminden/alien/blob/v0.1.6/docs/usage.md) for the full configuration reference.

## Source Inputs

ALIEN 0.1.6 supports managed MSigDB and Enrichr download/cache sources and local file sources:

- `msigdb_remote`: a Python downloader/reader for configured MSigDB release archives, with normalized Parquet caches for repeated builds.
- `enrichr_remote`: a Python downloader/reader for Enrichr libraries by public library name.
- `symbol_gmt`: a GMT file whose members are gene symbols.

Additional source metadata fields are documented in [docs/usage.md](https://github.com/deminden/alien/blob/v0.1.6/docs/usage.md).

## Python usage

```python
from alien import build

result = build("config.yml", workers=16)
print(result.namespaces)
```

Preset names and config dictionaries work too:

```python
result = build("pathways", outdir="results/pathways", output_mode="minimal")
```

## MSigDB

```yaml
sources:
  - type: msigdb_remote
    version: "2026.1"
    db_species: HS
    collection: C2
```

This stores the configured MSigDB release archive under `project.source_dir/msigdb_remote/` by default, or a source-level `cache_dir` override when provided. ALIEN verifies the archive by MD5, extracts the matching RDS files, and writes normalized Parquet memberships for later builds.

Remote caches are reused by default. ALIEN also caches normalized MSigDB memberships and prepared NCBI rescue maps, so repeated large builds avoid expensive source-format conversion. Use `alien build --force-download` or per-source `force: true` to refresh managed downloads such as MSigDB, Enrichr, GENCODE, HGNC, and NCBI resources. For publication configs, prefer exact Enrichr library names and keep `metadata/source_manifest.tsv` with the released GMTs. Enabled sources are always required; use `enabled: false` to exclude a source deliberately.

## Scope

The 0.1.6 release officially supports human gene sets using HGNC symbols, Python MSigDB and Enrichr cache integration, optimized repeated-build caches, and Ensembl-style target namespaces. The code is organized so broader namespace integrations can be added later without tying the package to any single downstream analysis project.

## Contributing

Contributions are welcome. Useful areas include additional tests, documentation, source adapters, target namespace adapters, mapping-audit improvements, curated filtering/source-priority defaults, and validation against established gene-set resources. See [docs/development.md](https://github.com/deminden/alien/blob/v0.1.6/docs/development.md) for development setup and current future plans.
