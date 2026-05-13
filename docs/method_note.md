# ALIEN 0.1.2 Method Note

ALIEN builds gene-set GMT files when source libraries and analysis namespaces do not line up cleanly. Source collections often mix current symbols, retired symbols, aliases, source Ensembl IDs, and occasional non-gene tokens. Downstream workflows, meanwhile, may require a specific target annotation release.

ALIEN treats integration as an audited projection problem. Every input source is normalized into one membership table, mapped through deterministic human gene resources, filtered, deduplicated, and written with QC records that explain what happened.

## Design

ALIEN 0.1.2 accepts three source styles:

- Managed MSigDB release archives through the Python `msigdb_remote` reader.
- Managed Enrichr libraries through the Python `enrichr_remote` reader.
- Local MSigDB-like tables with symbols and optional source identifiers.
- Symbol-only GMT libraries, including local Enrichr-style exports.

The builder creates a current-symbol namespace for display/review and any number of configured target namespaces. In 0.1.2 the implemented target adapter is `ensembl_gtf`: each target is defined by one Ensembl-style annotation GTF. The GTF can come from GENCODE by version or from another annotation origin via a local path or URL, as long as gene ID and symbol attributes are present or mapped in config. The optional `restrict_to` field is an advanced dataset-specific narrowing filter, not the main target path.

## Auditing

Mapping uses HGNC current, previous, and alias symbols first. Source Ensembl IDs that are absent from the target IDs can be checked against the Ensembl archive cache. NCBI Gene history/info can be enabled as a lower-confidence rescue layer for configured legacy sources. Ambiguous mappings are not guessed; they are written to audit tables.

Before mapping, ALIEN checks that each `term_id` refers to one source-term identity. Conflicting reuse of the same identifier across libraries is treated as an input error by default, because otherwise memberships would be merged while one term's metadata silently wins. When this happens, `metadata/term_id_collisions.tsv` records the conflicting identities for repair.

Filtering records dropped non-gene tokens, broad configured terms, size-filtered terms, and redundancy removals. Redundancy filtering first collapses exact duplicate memberships, then clusters terms whose Jaccard similarity is above the configured cutoff within each namespace and family. A representative term is chosen by source priority first, then mapped size and deterministic name tie-breaks.

Source priority is configured per term family because different biological domains have different preferred authorities. For example, pathway outputs can prefer Reactome over broader ontology-derived terms, while disease outputs can prefer curated disease sources over noisier text-mined libraries. Priority only decides between redundant alternatives; it does not remove non-overlapping terms from lower-priority sources.

The primary GMT output is intentionally simple: one combined GMT per namespace. The surrounding metadata explains how each term reached that output.

## Reproducibility

ALIEN writes source provenance, package versions, target coverage summaries, mapping summaries, and warning logs. The source cache remains outside the package so license-sensitive resources can be managed by each user.

## Limitations

ALIEN 0.1.2 is scoped to human gene-set integration with HGNC symbols, MSigDB and Enrichr local/remote caches, and Ensembl-style target namespaces. Non-Ensembl target ID systems, such as Entrez or UniProt GMT output, are future target adapters.
