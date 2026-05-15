#!/usr/bin/env Rscript

# GTEx v11 Thyroid sex contrast with DESeq2 and multilevel rsfgsea
# The script assumes the R environment and rsfgsea binary are already available
# Install rsfgsea with: cargo install rsfgsea

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(BiocParallel)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file, mustWork = FALSE)) else "scripts"
source(file.path(script_dir, "sex_contrast_common.R"))

# Main workflow
main <- function() {
# Parse command-line options
default_gmt <- paste(
  c(
    "results/alien_outputs/pathways/gmt/gtex_v11_gencode47.gmt",
    "results/alien_outputs/disease_phenotype/gmt/gtex_v11_gencode47.gmt",
    "results/alien_outputs/function_location/gmt/gtex_v11_gencode47.gmt"
  ),
  collapse = ","
)
default_build_gmt_configs <- paste(
  c(
    "examples/pathways.yml",
    "examples/disease_phenotype.yml",
    "examples/function_location.yml"
  ),
  collapse = ","
)
args <- parse_args(list(
  "data-dir" = "data/gtex",
  "outdir" = "docs/sex_contrast/gtex_thyroid",
  "results-dir" = "results/sex_contrast/gtex_thyroid",
  "gmt" = default_gmt,
  "build-gmt-configs" = default_build_gmt_configs,
  "alien" = Sys.getenv("ALIEN_BIN", unset = "alien"),
  "tcga-recount3-gtf" = "data/recount3/human.gene_sums.G029.gtf.gz",
  "tcga-output-genes" = "data/recount3/tcga_gencode_v29_output_genes.tsv.gz",
  "tcga-recount3-gtf-url" = "https://recount-opendata.s3.amazonaws.com/recount3/release/human/annotations/gene_sums/human.gene_sums.G029.gtf.gz",
  "tissue" = "Thyroid",
  "workers" = "1",
  "min-count" = "10",
  "min-size" = "10",
  "max-size" = "500",
  "nperm-simple" = "1000",
  "eps" = "1e-50",
  "seed" = "42",
  "rsfgsea" = Sys.getenv("RSFGSEA_BIN", unset = "rsfgsea"),
  "gene-reads-url" = "https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/counts-by-tissue/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_reads.parquet",
  "sample-attributes-url" = "https://storage.googleapis.com/adult-gtex/annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt",
  "subject-phenotypes-url" = "https://storage.googleapis.com/adult-gtex/annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt"
))

if (isTRUE(args[["help"]])) {
  cat(
    "Usage: Rscript scripts/gtex_thyroid_sex_contrast.R [options]\n\n",
    "Main options:\n",
    "  --data-dir PATH        GTEx cache directory [data/gtex]\n",
    "  --outdir PATH          Report and plot directory [docs/sex_contrast/gtex_thyroid]\n",
    "  --results-dir PATH     TSV/rank output directory [results/sex_contrast/gtex_thyroid]\n",
    "  --gmt PATHS            Comma-separated ALIEN GMT files for GTEx target namespace\n",
    "  --build-gmt-configs PATHS  Comma-separated ALIEN configs used when default GMTs are missing\n",
    "  --alien PATH|NAME      alien CLI used to build missing default GMTs [alien]\n",
    "  --tcga-output-genes PATH  recount3 G029 output-gene list for shared example configs\n",
    "  --rsfgsea PATH|NAME    rsfgsea binary; PATH is optional when cargo install put it on PATH\n",
    "  --workers N            DESeq2/rsfgsea workers [1]\n",
    "  --nperm-simple N       Simple-stage permutations used by multilevel rsfgsea [1000]\n",
    sep = ""
  )
  quit(status = 0)
}

# Validate analysis parameters
workers <- as_integer(args[["workers"]], "--workers")
min_count <- as_integer(args[["min-count"]], "--min-count")
min_size <- as_integer(args[["min-size"]], "--min-size")
max_size <- as_integer(args[["max-size"]], "--max-size")
nperm_simple <- as_integer(args[["nperm-simple"]], "--nperm-simple")
seed <- as_integer(args[["seed"]], "--seed")
eps <- as.numeric(args[["eps"]])
if (!is.finite(eps) || eps <= 0) {
  stop("--eps must be a positive number.", call. = FALSE)
}
gmt_paths <- parse_paths(args[["gmt"]])
if (length(gmt_paths) == 0) {
  stop("--gmt must contain at least one GMT path.", call. = FALSE)
}
build_gmt_configs <- parse_paths(args[["build-gmt-configs"]])
if (!identical(args[["gmt"]], default_gmt) && identical(args[["build-gmt-configs"]], default_build_gmt_configs)) {
  build_gmt_configs <- character()
}
design_formula <- "~ AGE + DTHHRDY + SMRIN + SMTSISCH + SMCENTER + sex"
design <- ~ AGE + DTHHRDY + SMRIN + SMTSISCH + SMCENTER + sex

# Set folder paths
data_dir <- args[["data-dir"]]
outdir <- args[["outdir"]]
results_dir <- args[["results-dir"]]
tables_dir <- file.path(results_dir, "tables")
plots_dir <- outdir
ranks_dir <- file.path(results_dir, "ranks")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ranks_dir, recursive = TRUE, showWarnings = FALSE)

gene_reads_file <- file.path(data_dir, basename(args[["gene-reads-url"]]))
sample_attributes_file <- file.path(data_dir, basename(args[["sample-attributes-url"]]))
subject_phenotypes_file <- file.path(data_dir, basename(args[["subject-phenotypes-url"]]))

# Load GTEx metadata and source files
message("1/7 Preparing GTEx inputs")
download_if_missing(args[["gene-reads-url"]], gene_reads_file)
download_if_missing(args[["sample-attributes-url"]], sample_attributes_file)
download_if_missing(args[["subject-phenotypes-url"]], subject_phenotypes_file)
if (any(!file.exists(gmt_paths)) && length(build_gmt_configs) > 0) {
  ensure_recount3_g029_output_genes(args[["tcga-recount3-gtf"]], args[["tcga-output-genes"]], args[["tcga-recount3-gtf-url"]])
}
ensure_gmt_files(gmt_paths, build_gmt_configs, workers, args[["alien"]])

sample_attributes <- fread(sample_attributes_file)
subject_phenotypes <- fread(subject_phenotypes_file)
stop_if_missing_columns(sample_attributes, c("SAMPID", "SMTSD"), "GTEx sample attributes")

sample_attributes[, SUBJID := sub("^(GTEX-[^-]+)-.*$", "\\1", SAMPID)]
sex_column <- find_sex_column(subject_phenotypes)
subject_phenotypes[, sex := normalize_sex(get(sex_column))]
subject_phenotypes <- subject_phenotypes[!is.na(sex), .(SUBJID, sex, AGE, DTHHRDY)]

# Select Thyroid samples
metadata <- merge(
  sample_attributes[SMTSD == args[["tissue"]]],
  subject_phenotypes,
  by = "SUBJID"
)

parquet_columns <- ParquetFileReader$create(gene_reads_file)$GetSchema()$names
sample_ids <- intersect(metadata$SAMPID, parquet_columns)
metadata <- metadata[SAMPID %in% sample_ids]
metadata <- metadata[!duplicated(SAMPID)]
design_columns <- c("AGE", "DTHHRDY", "SMRIN", "SMTSISCH", "SMCENTER", "sex")
stop_if_missing_columns(metadata, design_columns, "GTEx metadata")
metadata <- metadata[complete_design_rows(metadata, design_columns)]

if (nrow(metadata) < 6 || length(unique(metadata$sex)) < 2) {
  stop("Need at least two sex groups and enough GTEx samples for the contrast.", call. = FALSE)
}

message(sprintf(
  "Using %s %s samples: %s Female, %s Male",
  nrow(metadata),
  args[["tissue"]],
  sum(metadata$sex == "Female"),
  sum(metadata$sex == "Male")
))

sample_columns <- metadata$SAMPID
deseq_path <- file.path(tables_dir, "gtex_thyroid_female_vs_male_deseq2.tsv")
ranks_path <- file.path(ranks_dir, "gtex_thyroid_female_vs_male.ranks.tsv")
summary_path <- file.path(tables_dir, "gtex_thyroid_run_summary.tsv")
cache_expected <- list(
  design_formula = design_formula,
  tissue = args[["tissue"]],
  min_count = min_count
)
reuse_deseq <- all(file.exists(c(deseq_path, ranks_path))) && summary_matches(summary_path, cache_expected)
force_enrichment <- !reuse_deseq

if (reuse_deseq) {
# Reuse heavy DESeq2 outputs when available
  message("2/7 Reusing existing DESeq2 table and rank file")
  deseq_table <- fread(deseq_path)
  ranks <- fread(ranks_path, header = FALSE, col.names = c("gene", "score"))
} else {
# Read counts only for selected samples
  message("2/7 Reading Thyroid counts from parquet")
  counts_dt <- as.data.table(read_parquet(gene_reads_file, col_select = tidyselect::all_of(c("Name", sample_columns))))
  stop_if_missing_columns(counts_dt, "Name", "GTEx gene reads parquet")

  gene_ids_versioned <- counts_dt$Name
  count_matrix <- round(as.matrix(counts_dt[, sample_columns, with = FALSE]))
  storage.mode(count_matrix) <- "integer"
  rownames(count_matrix) <- gene_ids_versioned
  colnames(count_matrix) <- sample_columns

  coldata <- data.frame(
    row.names = sample_columns,
    AGE = factor(metadata$AGE),
    DTHHRDY = factor(metadata$DTHHRDY),
    SMRIN = as.numeric(metadata$SMRIN),
    SMTSISCH = as.numeric(metadata$SMTSISCH),
    SMCENTER = factor(metadata$SMCENTER),
    sex = factor(metadata$sex, levels = c("Male", "Female"))
  )

# Run DESeq2 sex contrast on raw counts
  message(sprintf("3/7 Running DESeq2 Female vs Male contrast with design %s", design_formula))
  assert_full_rank_design(design, coldata)
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = design)
  dds <- dds[rowSums(counts(dds)) >= min_count, ]
  bpparam <- if (workers > 1) MulticoreParam(workers) else SerialParam()
  dds <- DESeq(dds, quiet = TRUE, parallel = workers > 1, BPPARAM = bpparam)
  res <- results(dds, contrast = c("sex", "Female", "Male"))

  deseq_table <- as.data.table(as.data.frame(res), keep.rownames = "gene_id_versioned")
  deseq_table[, gene_id := strip_ensembl_version(gene_id_versioned)]
  setcolorder(deseq_table, c("gene_id", "gene_id_versioned", setdiff(names(deseq_table), c("gene_id", "gene_id_versioned"))))
  fwrite(deseq_table, deseq_path, sep = "\t")

# Write ranked statistic for rsfgsea
  message("4/7 Writing ranked genes for rsfgsea")
  ranks <- deseq_table[
    is.finite(stat) & !is.na(gene_id) & nzchar(gene_id),
    .(gene = gene_id, score = stat)
  ][order(-score)]
  ranks <- ranks[!duplicated(gene)]
  fwrite(ranks, ranks_path, sep = "\t", col.names = FALSE)
}

fwrite(data.table(
  metric = c("samples", "female_samples", "male_samples", "ranked_genes", "design_formula", "tissue", "min_count"),
  value = c(nrow(metadata), sum(metadata$sex == "Female"), sum(metadata$sex == "Male"), nrow(ranks), design_formula, args[["tissue"]], min_count)
), summary_path, sep = "\t")

# Run enrichment for each ALIEN GMT collection
message("5/7 Running multilevel rsfgsea")
rsfgsea_bin <- resolve_rsfgsea(args[["rsfgsea"]])
enrichment_runs <- list()
plot_theme_specs <- list(
  pathways = list(
    theme_spec("immune and antigen presentation", "IMMUN|ANTIGEN|MHC|B_CELL|B CELL|T_CELL|T CELL|LYMPHOID|LEUKOCYTE|CYTOKINE")
  ),
  disease_phenotype = list(
    theme_spec("thyroid disease and autoimmunity", "THYROID|GRAVES|HASHIMOTO|AUTOIMMUN")
  )
)
for (gmt_path in gmt_paths) {
  label <- gmt_label(gmt_path)
  enrichment_path <- file.path(tables_dir, sprintf("gtex_thyroid_%s_female_vs_male_rsfgsea.tsv", label))
  if (file.exists(enrichment_path) && !force_enrichment) {
    message(sprintf("Reusing existing rsfgsea table for %s", label))
  } else {
    message(sprintf("Running rsfgsea for %s", label))
    run_rsfgsea(
      rsfgsea_bin = rsfgsea_bin,
      ranks_path = ranks_path,
      gmt_path = gmt_path,
      output_path = enrichment_path,
      min_size = min_size,
      max_size = max_size,
      nperm_simple = nperm_simple,
      eps = eps,
      seed = seed,
      workers = workers
    )
  }

# Build plots and summaries
  message(sprintf("6/7 Plotting enrichment summaries for %s", label))
  enrichment <- read_enrichment_table(enrichment_path)
  theme_specs <- plot_theme_specs[[label]]
  if (is.null(theme_specs)) theme_specs <- list()
  plot_prefix <- plot_slug(label)
  cleanup_running_score_plots(plots_dir, plot_prefix)
  running_plots <- character()
  if (length(theme_specs) > 0) {
    gmt_terms <- read_gmt(gmt_path)
    running_plots <- plot_representative_running_scores(
      enrichment = enrichment,
      ranks = ranks,
      gmt_terms = gmt_terms,
      plots_dir = plots_dir,
      prefix = plot_prefix,
      theme_specs = theme_specs
    )
  }

  enrichment_runs[[label]] <- list(
    label = label,
    gmt = gmt_path,
    enrichment = enrichment,
    enrichment_path = enrichment_path,
    figure_paths = running_plots
  )
}

# Write Markdown report with links to result tables
message("7/7 Plotting global enrichment summaries and writing Markdown report")
collection_plot <- file.path(plots_dir, "collection_significance.png")
plot_collection_significance(enrichment_runs, collection_plot)

report_path <- file.path(outdir, "analysis.md")
write_report(
  path = report_path,
  title = "GTEx Thyroid Sex Contrast",
  summary_lines = c(
    "GTEx v11 thyroid samples were contrasted by sex using DESeq2 on raw counts, adjusting for age group, Hardy death class, RNA integrity, tissue ischemic time, and sequencing center. Genes were ranked by the Wald statistic for Female vs Male and tested against three ALIEN GTEx/Gencode GMT collections with multilevel rsfgsea: pathways, disease/phenotype, and function/location. Cancer dependency is intentionally not used for GTEx.",
    "",
    sprintf("- Samples: %s total (%s Female, %s Male).", nrow(metadata), sum(metadata$sex == "Female"), sum(metadata$sex == "Male")),
    sprintf("- Ranked genes: %s.", nrow(ranks)),
    sprintf("- DESeq2 design: `%s`.", design_formula),
    format_gmt_summary(gmt_paths),
    sprintf("- rsfgsea: `%s` in multilevel mode.", basename(rsfgsea_bin)),
    "- Direction: positive NES means Female > Male; negative NES means Male > Female."
  ),
  enrichment_runs = enrichment_runs,
  deseq_path = deseq_path,
  ranks_path = ranks_path,
  collection_plot_path = collection_plot,
  outdir = outdir,
  collection_note = "The figure below gives a collection-level view of significant enrichment. Detailed term-level results are shown in the tables; figures are limited to the clearest thyroid immune and autoimmune examples."
)

message(sprintf("Done: %s", report_path))
}


main()
