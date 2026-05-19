# ALIEN 0.1.6 Method Note

ALIEN builds gene-set GMT files when source libraries and analysis namespaces do not line up cleanly. Source collections often mix current symbols, retired symbols, aliases, source Ensembl IDs, and occasional non-gene tokens. Downstream workflows, meanwhile, may require a specific target annotation release.

ALIEN treats integration as an audited projection problem. Every input source is normalized into one membership table, mapped through deterministic human gene resources, filtered, deduplicated, and written with QC records that explain what happened.

## Design

ALIEN 0.1.6 accepts these source styles:

- Managed MSigDB release archives through the Python `msigdb_remote` reader.
- Managed Enrichr libraries through the Python `enrichr_remote` reader.
- Local MSigDB-like tables with symbols and additional source identifiers.
- Symbol-only GMT libraries, including local Enrichr-style exports.

The builder creates a current-symbol namespace for display/review and any number of configured target namespaces. In 0.1.6 the implemented target adapter is `ensembl_gtf`: each target separates the output gene set from the annotation helper. By default, one Ensembl-style annotation GTF provides both; with `output_genes`, a dataset file/list defines the final output gene namespace and the GTF is primarily a symbol/metadata helper; with `gene_filter`, the final namespace is the intersection of annotation IDs and configured filter IDs. The GTF can come from GENCODE by version or from another annotation origin via a local path or URL, as long as gene ID and symbol attributes are present or mapped in config. In `output_genes`, `id_column` must contain Ensembl IDs; `symbol_column` is optional helper metadata for output genes absent from the GTF.

For dataset-specific output namespaces, the primary annotation can receive metadata fallbacks without being replaced. `annotation.metadata_fallbacks` loads a secondary GTF, adds only output IDs missing from the primary annotation, and records those rows as fallback metadata. When the same symbol is present in both primary and fallback annotations, the primary annotation remains the preferred mapping candidate. This keeps official releases such as GENCODE as the mapping authority while allowing resources such as recount3 to explain measured IDs absent from the official release file.

## Auditing

Mapping uses HGNC current, previous, and alias symbols first. By default, source Ensembl IDs that are absent from the target annotation, any configured target output genes/filter, and HGNC's current Ensembl IDs are checked against the Ensembl archive cache. NCBI Gene history/info rescue is also enabled by default, but restricted to configured legacy source tags. Both fallback layers can be disabled in config. Ambiguous mappings are not guessed; they are written to audit tables.

Before mapping, ALIEN checks that each `term_id` refers to one source-term identity. Conflicting reuse of the same identifier across libraries is treated as an input error by default, because otherwise memberships would be merged while one term's metadata silently wins. When this happens, `metadata/term_id_collisions.tsv` records the conflicting identities for repair.

Filtering records dropped non-gene tokens, broad configured terms, size-filtered terms, and redundancy removals. Redundancy filtering first collapses exact duplicate memberships, then clusters terms whose Jaccard similarity is above the configured cutoff within each namespace and family. A representative term is chosen by source priority first, then mapped size and deterministic name tie-breaks.

Source priority is configured per term family because different biological domains have different preferred authorities. For example, pathway outputs can prefer Reactome over broader ontology-derived terms, while disease outputs can prefer curated disease sources over noisier text-mined libraries. Priority only decides between redundant alternatives; it does not remove non-overlapping terms from lower-priority sources.

The primary GMT output is intentionally simple: one combined GMT per namespace. The surrounding metadata explains how each term reached that output.

## Reproducibility

ALIEN writes a compact source manifest, source provenance, package versions, target coverage summaries, mapping summaries, and warning logs. `metadata/source_manifest.tsv` is the main article/protocol record for exact source collections: it lists each resolved source collection once, with term/member counts and Enrichr exact/regex resolution details when relevant. `metadata/target_metadata_fallbacks.tsv` records any secondary annotation files used to add metadata for missing output IDs. Enabled sources are required; missing source libraries stop the build rather than silently weakening the output. The source cache remains outside the package so license-sensitive resources can be managed by each user.

## Cache And Performance

Managed downloads are treated as reproducible source artifacts, while ALIEN keeps its own derived caches for repeated builds. MSigDB remote releases are verified and unpacked from the official archive, then matching RDS tables are normalized once into ALIEN Parquet membership caches. NCBI Gene rescue maps are also stored in a prepared JSON form after first construction.

Large builds reuse those derived caches unless `--force-download` or source-level `force: true` is requested. Namespace mapping uses cached symbol/source-ID resolutions and can map multiple target namespaces in parallel when workers are available. Redundancy filtering uses deterministic set logic and union-find clustering rather than a graph dependency. These optimizations are intended to change runtime, not output semantics; audit tables and GMT outputs remain deterministic for the same inputs and configuration.

## Limitations

ALIEN 0.1.6 is scoped to human gene-set integration with HGNC symbols, MSigDB and Enrichr local/remote caches, optimized repeated-build caches, and Ensembl-style target namespaces. Non-Ensembl target ID systems, such as Entrez or UniProt GMT output, are future target adapters.
