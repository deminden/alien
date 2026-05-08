# Usage

ALIEN builds combined GMT files from configured source libraries and projects them into configured target namespaces. A build is controlled by one YAML file.

## Run

```bash
alien build --config configs/production.yml --workers 4
```

The Python API is equivalent:

```python
from alien import build

result = build("configs/production.yml", workers=4)
print(result.namespaces)
```

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

## Source Types

`msigdb_remote` reads a managed `msigdbr` release archive from Zenodo:

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

Downloaded Enrichr metadata is cached at `data/alien_sources/enrichr/datasetStatistics.json`; GMT files are cached as `data/alien_sources/enrichr/<library>.gmt`.

Cached remote resources are reused by default. To refresh MSigDB archives, Enrichr metadata/GMTs, GENCODE annotations, HGNC/NCBI resources, and Ensembl archive caches, run:

```bash
alien build --config configs/production.yml --force-download
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

This matters for regex-matched Enrichr libraries: ALIEN chooses the latest matching library from the cached `datasetStatistics` file unless the cache is refreshed.

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

`canonical_tsv` reads ALIEN's normalized table format. Required columns are `term_id` and `gene_symbol`; optional metadata columns include `source`, `source_tag`, `collection`, `family`, `aspect`, `gene_id`, and `gene_id_namespace`.

## Targets

The implemented target adapter is `ensembl_gtf`. It projects source terms into Ensembl stable gene IDs using a target annotation GTF. The GTF must contain `gene` records with `gene_id` and `gene_name` attributes.

For built-in GENCODE releases, `source` and `version` are enough:

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

Advanced narrowing: `restrict_to` can limit a target GMT to Ensembl IDs present in a dataset matrix or gene list. Normal builds should omit it.

```yaml
targets:
  - name: human_gencode47_study_only
    type: ensembl_gtf
    annotation:
      source: GENCODE
      version: "47"
    restrict_to: data/study/count_matrix.tsv
```

If the matrix is tabular and the gene ID column is known, use the explicit form:

```yaml
targets:
  - name: human_gencode47_study_only
    type: ensembl_gtf
    annotation:
      source: GENCODE
      version: "47"
    restrict_to:
      path: data/study/expression.tsv.gz
      column: feature_id
```

Without `column`, ALIEN tries common gene-ID column names and then falls back to the first tabular column. For Parquet files, `column: index` reads the row index.

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
            "restrict_to": {"path": "data/study/expression.tsv.gz", "column": "feature_id"},
        }
    ],
}

result = build(cfg, outdir="data/alien_gmt", workers=4)
```

For code-only workflows, `restrict_to` can also carry inline Ensembl IDs:

```python
"restrict_to": {"ids": ["ENSG00000141510.18", "ENSG000002"]}
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

Terms with Jaccard similarity greater than `jaccard_cutoff` are treated as redundant. One representative is kept using `source_priority`, mapped term size, and deterministic name tie-breaks.

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
- Optional NCBI Gene history/info rescue.
- Optional Ensembl archive lookup for source Ensembl IDs absent from the target IDs.

Ambiguous mappings are dropped rather than guessed and are written to audit tables.

## Outputs

Main GMT outputs:

```text
gmt/<target_namespace>.gmt
gmt/symbols.gmt
```

Main audit outputs:

```text
metadata/term_manifest.tsv.gz
metadata/gene_mapping_<target>.tsv.gz
metadata/removed_terms_size_filter.tsv.gz
metadata/removed_terms_redundancy.tsv.gz
metadata/unmapped_genes.tsv.gz
metadata/ambiguous_gene_mappings.tsv.gz
metadata/source_provenance.json
qc/redundancy_summary.tsv
qc/mapping_summary.tsv
qc/warnings.txt
```
