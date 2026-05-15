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

Prepare a config with sources and a target annotation, then run:

```bash
alien build --config examples/cancer_dependency.yml --workers 16
```

Ready-made configs are available for GTEx v11 / GENCODE v47 and TCGA recount3 / GENCODE v29 targets. The first three configs write both GTEx and TCGA GMTs; the cancer dependency config is TCGA-only. TCGA configs use GENCODE v29 as the primary annotation helper and recount3 G029 as a metadata fallback for measured IDs missing from GENCODE.

- [examples/pathways.yml](https://github.com/deminden/alien/blob/v0.1.5/examples/pathways.yml): Reactome, WikiPathways, KEGG MEDICUS, and GO biological process terms.
- [examples/function_location.yml](https://github.com/deminden/alien/blob/v0.1.5/examples/function_location.yml): GO molecular function and cellular component terms.
- [examples/disease_phenotype.yml](https://github.com/deminden/alien/blob/v0.1.5/examples/disease_phenotype.yml): HPO, DisGeNET, ClinVar, GWAS Catalog, and Jensen disease libraries.
- [examples/cancer_dependency.yml](https://github.com/deminden/alien/blob/v0.1.5/examples/cancer_dependency.yml): cancer and dependency signatures for TCGA recount3.

The primary outputs are:

- `gmt/<target_namespace>.gmt`: combined GMT for each configured target namespace.
- `metadata/`: source manifest, term manifest, gene mapping tables, removed terms, unmapped genes, ambiguity logs, and provenance.
- `qc/`: collection, mapping, redundancy, target coverage, and warning summaries.

The two filesystem roots are configured in YAML: `project.source_dir` is for downloaded/cached source and mapping resources, while `project.outdir` is the output root containing `gmt/`, `metadata/`, and `qc/`.

Downstream enrichment reports using these GMTs are included in [docs/sex_contrast/gtex_thyroid/analysis.md](https://github.com/deminden/alien/blob/v0.1.5/docs/sex_contrast/gtex_thyroid/analysis.md) and [docs/sex_contrast/tcga_lung/analysis.md](https://github.com/deminden/alien/blob/v0.1.5/docs/sex_contrast/tcga_lung/analysis.md). Scripts to regenerate the report data, text, and figures are available in [scripts/](https://github.com/deminden/alien/tree/v0.1.5/scripts/).

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

See [docs/usage.md](https://github.com/deminden/alien/blob/v0.1.5/docs/usage.md) for the full configuration reference.

## Source Inputs

ALIEN 0.1.5 supports managed MSigDB and Enrichr download/cache sources and local file sources:

- `msigdb_remote`: a Python downloader/reader for configured MSigDB release archives, with normalized Parquet caches for repeated builds.
- `enrichr_remote`: a Python downloader/reader for Enrichr libraries by public library name.
- `symbol_gmt`: a GMT file whose members are gene symbols.

Additional source metadata fields are documented in [docs/usage.md](https://github.com/deminden/alien/blob/v0.1.5/docs/usage.md).

## Python API

```python
from alien import build

result = build("examples/pathways.yml", workers=16)
print(result.namespaces)
```

You can also pass the same configuration as a Python dictionary, including target output genes such as `output_genes: {"path": "expression.tsv.gz", "id_column": "feature_id", "symbol_column": "gene_symbol"}` for a known Ensembl ID field.

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

The 0.1.5 release officially supports human gene sets using HGNC symbols, Python MSigDB and Enrichr cache integration, optimized repeated-build caches, and Ensembl-style target namespaces. The code is organized so broader namespace integrations can be added later without tying the package to any single downstream analysis project.

## Contributing

Contributions are welcome. Useful areas include additional tests, documentation, source adapters, target namespace adapters, mapping-audit improvements, curated filtering/source-priority defaults, and validation against established gene-set resources. See [docs/development.md](https://github.com/deminden/alien/blob/v0.1.5/docs/development.md) for development setup and current future plans.
