# ALIEN

**ALIEN: Audited Library Integration for External Namespaces** builds namespace-specific GMT libraries for human gene-set workflows.

ALIEN 0.1.1 is centered on one job: take configured source libraries, normalize their gene memberships into a canonical table, project them into configured target namespaces, and write combined GMT files with audit metadata.

## Install

ALIEN is not on PyPI yet. Install it from the repository:

```bash
git clone https://github.com/deminden/alien.git
cd alien
python -m pip install -e .
```

## Quick Start

Prepare a config with sources and a target annotation, then run:

```bash
alien build --config configs/production.yml --outdir data/alien_gmt --workers 4
```

The primary outputs are:

- `gmt/<target_namespace>.gmt`: combined GMT for each namespace, including `symbols.gmt` when enabled.
- `metadata/`: term manifest, gene mapping tables, removed terms, unmapped genes, ambiguity logs, and provenance.
- `qc/`: collection, mapping, redundancy, target coverage, and warning summaries.

## How It Works

ALIEN builds one canonical membership table from configured sources, audits source gene symbols against human mapping resources, then projects each term into the requested output namespaces.

The target namespace is defined in the config. For an Ensembl-style namespace, provide a target name and a GTF annotation:

```yaml
outputs:
  include_symbols: true

targets:
  - name: human_gencode47
    type: ensembl_gtf
    annotation:
      source: GENCODE
      version: "47"
```

This writes `gmt/human_gencode47.gmt`. When `include_symbols` is enabled, ALIEN also writes `gmt/symbols.gmt`.

For built-in GENCODE releases, `version` is enough; ALIEN downloads and caches the GTF under `data/alien_sources/gencode/` if it is missing. You can still provide `annotation.path` to pin a local file explicitly.

`source: GENCODE` is not special to the config shape; it is only the annotation origin used in this example. Another Ensembl-style GTF origin can use the same adapter:

```yaml
targets:
  - name: my_ensembl_namespace
    type: ensembl_gtf
    annotation:
      source: MyAnnotationProvider
      version: "2026-05"
      path: data/alien_sources/my_provider/genes.gtf.gz
```

Fully different target ID systems, such as Entrez or UniProt GMT output, are planned as future target adapters.

Advanced use: add `restrict_to: path/to/matrix_or_gene_list.tsv` only when you deliberately want to narrow a target GMT to genes present in one dataset. Normal builds should omit it and use all genes from the target annotation.

The build keeps audit outputs beside the GMTs so each namespace projection can be traced back to source terms, symbol repairs, unmapped genes, filters, redundancy decisions, and provenance.

## Source Inputs

ALIEN 0.1.1 supports managed MSigDB download/cache sources and local file sources:

- `msigdb_cache`: a directory of `msigdbr_<SOURCE_TAG>.tsv.gz` files.
- `msigdb_remote`: a Python downloader/reader for the current `msigdbr` Zenodo release cache.
- `msigdb_tsv`: a MSigDB-like TSV with term names, gene symbols, and optional source Ensembl IDs.
- `symbol_gmt`: a GMT file whose members are gene symbols.
- `canonical_tsv`: ALIEN’s normalized membership schema.

Canonical tables require:

```text
term_id    gene_symbol
```

Recommended optional fields are `original_name`, `display_name`, `description`, `source`, `source_tag`, `collection`, `family`, `aspect`, `gene_id`, and `gene_id_namespace`. If `gene_id_namespace` identifies Ensembl IDs, ALIEN audits those IDs before falling back through symbol mapping.

## Python API

```python
from alien import build

result = build("configs/production.yml", outdir="data/alien_gmt", workers=4)
print(result.namespaces)
```

Only the small facade API is stable for 0.1.1. Lower-level modules are importable for experimentation but may change while the package grows.

## MSigDB Cache Preparation

ALIEN can read local MSigDB-style caches or download the current `msigdbr` release archive directly from Python with a source like:

```yaml
sources:
  - type: msigdb_remote
    version: "2026.1"
    db_species: HS
    collection: C2
```

The remote archive is cached under `data/alien_sources/msigdb_remote/`, verified by MD5, and read from the same collection-specific RDS files used by the R package.

The older R helper is kept for compatibility and comparison tests:

```bash
Rscript scripts/fetch_msigdb.R data/alien_sources/msigdb FALSE 26.1.0
```

## Scope

The 0.1.1 release officially supports human gene sets using HGNC symbols, Python MSigDB cache integration, and Ensembl-style target namespaces. The code is organized so broader namespace integrations can be added later without tying the package to any single downstream analysis project.
