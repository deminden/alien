# ALIEN 0.1.0 Method Note

ALIEN builds gene-set GMT files when source libraries and analysis namespaces do not line up cleanly. Source collections often mix current symbols, retired symbols, aliases, source Ensembl IDs, and occasional non-gene tokens. Downstream expression matrices, meanwhile, may use a specific Ensembl release or a dataset-specific gene universe.

ALIEN treats integration as an audited projection problem. Every input source is normalized into one membership table, mapped through deterministic human gene resources, filtered, deduplicated, and written with QC records that explain what happened.

## Design

ALIEN 0.1.0 accepts two source styles:

- MSigDB-like tables with symbols and optional source identifiers.
- Symbol-only GMT libraries, including local Enrichr-style exports.

The builder creates a current-symbol namespace for display/review and any number of configured Ensembl target namespaces. Each target defines its own annotation GTF and optional gene-universe file. When a target universe is supplied, matrix row IDs are the final membership authority; annotation files are helpers for symbol and metadata lookup.

## Auditing

Mapping uses HGNC current, previous, and alias symbols first. Source Ensembl IDs that are absent from a target universe can be checked against the Ensembl archive cache. NCBI Gene history/info can be enabled as a lower-confidence rescue layer for configured legacy sources. Ambiguous mappings are not guessed; they are written to audit tables.

Filtering records dropped non-gene tokens, broad configured terms, size-filtered terms, and redundancy removals. The primary GMT output is intentionally simple: one combined GMT per namespace. The surrounding metadata explains how each term reached that output.

## Reproducibility

ALIEN writes source provenance, package versions, target-universe summaries, mapping summaries, and warning logs. The source cache remains outside the package so license-sensitive resources can be managed by each user.

## Limitations

ALIEN 0.1.0 is scoped to human gene-set integration with HGNC symbols and Ensembl target namespaces. MSigDB downloading is still external; a Python-based downloader is planned for 0.1.1.
