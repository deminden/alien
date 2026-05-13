# ALIEN 0.1.3 Method Note

ALIEN builds gene-set GMT files when source libraries and analysis namespaces do not line up cleanly. Source collections often mix current symbols, retired symbols, aliases, source Ensembl IDs, and occasional non-gene tokens. Downstream workflows, meanwhile, may require a specific target annotation release.

ALIEN treats integration as an audited projection problem. Every input source is normalized into one membership table, mapped through deterministic human gene resources, filtered, deduplicated, and written with QC records that explain what happened.

## Design

ALIEN 0.1.3 accepts these source styles:

- Managed MSigDB release archives through the Python `msigdb_remote` reader.
- Managed Enrichr libraries through the Python `enrichr_remote` reader.
- Local MSigDB-like tables with symbols and optional source identifiers.
- Symbol-only GMT libraries, including local Enrichr-style exports.

The builder creates a current-symbol namespace for display/review and any number of configured target namespaces. In 0.1.3 the implemented target adapter is `ensembl_gtf`: each target is defined by one Ensembl-style annotation GTF plus, when needed, either a `gene_universe` file/list that defines the output ID namespace or a `gene_filter` file/list that restricts the annotation namespace by intersection. The GTF can come from GENCODE by version or from another annotation origin via a local path or URL, as long as gene ID and symbol attributes are present or mapped in config. When `gene_universe` is present, the GTF is primarily a symbol and metadata helper; the configured universe owns the output ID set.

## Auditing

Mapping uses HGNC current, previous, and alias symbols first. By default, source Ensembl IDs that are absent from the target annotation, any configured target gene universe/filter, and HGNC's current Ensembl IDs are checked against the Ensembl archive cache. NCBI Gene history/info rescue is also enabled by default, but restricted to configured legacy source tags. Both fallback layers can be disabled in config. Ambiguous mappings are not guessed; they are written to audit tables.

Before mapping, ALIEN checks that each `term_id` refers to one source-term identity. Conflicting reuse of the same identifier across libraries is treated as an input error by default, because otherwise memberships would be merged while one term's metadata silently wins. When this happens, `metadata/term_id_collisions.tsv` records the conflicting identities for repair.

Filtering records dropped non-gene tokens, broad configured terms, size-filtered terms, and redundancy removals. Redundancy filtering first collapses exact duplicate memberships, then clusters terms whose Jaccard similarity is above the configured cutoff within each namespace and family. A representative term is chosen by source priority first, then mapped size and deterministic name tie-breaks.

Source priority is configured per term family because different biological domains have different preferred authorities. For example, pathway outputs can prefer Reactome over broader ontology-derived terms, while disease outputs can prefer curated disease sources over noisier text-mined libraries. Priority only decides between redundant alternatives; it does not remove non-overlapping terms from lower-priority sources.

The primary GMT output is intentionally simple: one combined GMT per namespace. The surrounding metadata explains how each term reached that output.

## Reproducibility

ALIEN writes source provenance, package versions, target coverage summaries, mapping summaries, and warning logs. The source cache remains outside the package so license-sensitive resources can be managed by each user.

## Cache And Performance

Managed downloads are treated as reproducible source artifacts, while ALIEN keeps its own derived caches for repeated builds. MSigDB remote releases are verified and unpacked from the official archive, then matching RDS tables are normalized once into ALIEN Parquet membership caches. NCBI Gene rescue maps are also stored in a prepared JSON form after first construction.

Large builds reuse those derived caches unless `--force-download` or source-level `force: true` is requested. Namespace mapping uses cached symbol/source-ID resolutions and can map multiple target namespaces in parallel when workers are available. Redundancy filtering uses deterministic set logic and union-find clustering rather than a graph dependency. These optimizations are intended to change runtime, not output semantics; audit tables and GMT outputs remain deterministic for the same inputs and configuration.

## Limitations

ALIEN 0.1.3 is scoped to human gene-set integration with HGNC symbols, MSigDB and Enrichr local/remote caches, optimized repeated-build caches, and Ensembl-style target namespaces. Non-Ensembl target ID systems, such as Entrez or UniProt GMT output, are future target adapters.
