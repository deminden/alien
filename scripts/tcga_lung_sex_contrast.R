#!/usr/bin/env Rscript

# TCGA lung cancer sex contrast with DESeq2 and multilevel rsfgsea
# The script assumes the R environment and rsfgsea binary are already available
# Install rsfgsea with: cargo install rsfgsea

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(BiocParallel)
  library(recount3)
  library(SummarizedExperiment)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file, mustWork = FALSE)) else "scripts"
source(file.path(script_dir, "sex_contrast_common.R"))

options(timeout = 3600)

# Main workflow
main <- function() {
# Parse command-line options
default_gmt <- paste(
  c(
    "results/alien_outputs/pathways/gmt/tcga_recount3_gencode29.gmt",
    "results/alien_outputs/disease_phenotype/gmt/tcga_recount3_gencode29.gmt",
    "results/alien_outputs/function_location/gmt/tcga_recount3_gencode29.gmt",
    "results/alien_outputs/cancer_dependency/gmt/tcga_recount3_gencode29.gmt"
  ),
  collapse = ","
)
default_build_gmt_configs <- paste(
  c(
    "examples/pathways.yml",
    "examples/disease_phenotype.yml",
    "examples/function_location.yml",
    "examples/cancer_dependency.yml"
  ),
  collapse = ","
)
args <- parse_args(list(
  "cache-dir" = "data/tcga/recount3_cache_lung_sex",
  "outdir" = "docs/sex_contrast/tcga_lung",
  "results-dir" = "results/sex_contrast/tcga_lung",
  "gmt" = default_gmt,
  "build-gmt-configs" = default_build_gmt_configs,
  "alien" = Sys.getenv("ALIEN_BIN", unset = "alien"),
  "tcga-recount3-gtf" = "data/recount3/human.gene_sums.G029.gtf.gz",
  "tcga-output-genes" = "data/recount3/tcga_gencode_v29_output_genes.tsv.gz",
  "tcga-recount3-gtf-url" = "https://recount-opendata.s3.amazonaws.com/recount3/release/human/annotations/gene_sums/human.gene_sums.G029.gtf.gz",
  "projects" = "LUAD,LUSC",
  "sample-type" = "Primary Tumor",
  "annotation" = "gencode_v29",
  "recount3-url" = "https://recount-opendata.s3.amazonaws.com/recount3/release",
  "workers" = "1",
  "min-count" = "10",
  "min-size" = "10",
  "max-size" = "500",
  "nperm-simple" = "1000",
  "eps" = "1e-50",
  "seed" = "42",
  "rsfgsea" = Sys.getenv("RSFGSEA_BIN", unset = "rsfgsea")
))

if (isTRUE(args[["help"]])) {
  cat(
    "Usage: Rscript scripts/tcga_lung_sex_contrast.R [options]\n\n",
    "Main options:\n",
    "  --cache-dir PATH       recount3 cache directory [data/tcga/recount3_cache_lung_sex]\n",
    "  --outdir PATH          Report and plot directory [docs/sex_contrast/tcga_lung]\n",
    "  --results-dir PATH     TSV/rank output directory [results/sex_contrast/tcga_lung]\n",
    "  --gmt PATHS            Comma-separated ALIEN GMT files for TCGA target namespace\n",
    "  --build-gmt-configs PATHS  Comma-separated ALIEN configs used when default GMTs are missing\n",
    "  --alien PATH|NAME      alien CLI used to build missing default GMTs [alien]\n",
    "  --tcga-output-genes PATH  recount3 G029 output-gene list for shared example configs\n",
    "  --projects IDS         Comma-separated TCGA projects [LUAD,LUSC]\n",
    "  --sample-type NAME     TCGA sample type filter; use 'all' to disable [Primary Tumor]\n",
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
projects <- toupper(trimws(unlist(strsplit(args[["projects"]], ","))))
projects <- projects[nzchar(projects)]

if (length(projects) == 0) {
  stop("--projects must contain at least one TCGA project ID.", call. = FALSE)
}
design_formula <- if (length(projects) > 1) "~ project + age_at_diagnosis + sex" else "~ age_at_diagnosis + sex"

# Set folder paths
outdir <- args[["outdir"]]
results_dir <- args[["results-dir"]]
tables_dir <- file.path(results_dir, "tables")
plots_dir <- outdir
ranks_dir <- file.path(results_dir, "ranks")
dir.create(args[["cache-dir"]], recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ranks_dir, recursive = TRUE, showWarnings = FALSE)
if (any(!file.exists(gmt_paths)) && length(build_gmt_configs) > 0) {
  ensure_recount3_g029_output_genes(args[["tcga-recount3-gtf"]], args[["tcga-output-genes"]], args[["tcga-recount3-gtf-url"]])
}
ensure_gmt_files(gmt_paths, build_gmt_configs, workers, args[["alien"]])

read_cached_tcga_project <- function(cache_dir, project_id) {
  files <- list.files(file.path(cache_dir, project_id), full.names = TRUE)
  metadata_files <- grep("tcga\\.(tcga|recount_project|recount_qc|recount_seq_qc)\\.", files, value = TRUE)
  counts_file <- grep(sprintf("tcga\\.gene_sums\\.%s\\.G029\\.gz$", project_id), files, value = TRUE)[1]
  if (length(metadata_files) == 0 || is.na(counts_file)) {
    stop(sprintf("Cached recount3 files are incomplete for %s.", project_id), call. = FALSE)
  }
  names(metadata_files) <- basename(metadata_files)
  metadata <- as.data.table(recount3:::read_metadata(metadata_files))
  counts <- recount3:::read_counts(counts_file, samples = metadata$external_id)
  scale_factor <- recount3:::compute_scale_factors(metadata)
  counts <- round(counts * matrix(rep(scale_factor[colnames(counts)], each = nrow(counts)), ncol = ncol(counts)))
  list(counts = counts, metadata = metadata)
}

deseq_path <- file.path(tables_dir, "tcga_lung_female_vs_male_deseq2.tsv")
ranks_path <- file.path(ranks_dir, "tcga_lung_female_vs_male.ranks.tsv")
summary_path <- file.path(tables_dir, "tcga_lung_run_summary.tsv")
cache_expected <- list(
  design_formula = design_formula,
  projects = paste(projects, collapse = ","),
  sample_type = args[["sample-type"]],
  annotation = args[["annotation"]],
  min_count = min_count
)
reuse_deseq <- all(file.exists(c(deseq_path, ranks_path))) && summary_matches(summary_path, cache_expected)
force_enrichment <- !reuse_deseq

if (reuse_deseq) {
# Reuse heavy DESeq2 outputs when available
  message("1/7 Reusing existing TCGA DESeq2 table, rank file, and run summary")
  deseq_table <- fread(deseq_path)
  ranks <- fread(ranks_path, header = FALSE, col.names = c("gene", "score"))
  run_summary <- fread(summary_path)
} else {
# Resolve requested recount3 projects
message("1/7 Resolving TCGA projects in recount3")
annotation_choices <- annotation_options("human")
if (!(args[["annotation"]] %in% annotation_choices)) {
  stop(sprintf(
    "Annotation '%s' is not available in recount3. Available human annotations: %s",
    args[["annotation"]],
    paste(annotation_choices, collapse = ", ")
  ), call. = FALSE)
}

project_cache <- recount3_cache(file.path(args[["cache-dir"]], "projects"))
human_projects <- tryCatch(
  available_projects(recount3_url = args[["recount3-url"]], bfc = project_cache),
  error = function(error) {
    warning(
      sprintf("Could not refresh recount3 project index; using cached TCGA project homes for %s.", paste(projects, collapse = ", ")),
      call. = FALSE
    )
    data.frame(
      project = projects,
      organism = "human",
      file_source = "tcga",
      project_type = "data_sources",
      project_home = "data_sources/tcga",
      stringsAsFactors = FALSE
    )
  }
)
tcga_projects <- subset(
  human_projects,
  organism == "human" &
    file_source == "tcga" &
    project_type == "data_sources" &
    project %in% projects
)
missing_projects <- setdiff(projects, tcga_projects$project)
if (length(missing_projects) > 0) {
  stop(sprintf("Requested TCGA projects were not found in recount3: %s", paste(missing_projects, collapse = ", ")), call. = FALSE)
}
tcga_projects <- tcga_projects[match(projects, tcga_projects$project), , drop = FALSE]

# Load TCGA counts and metadata for each project
message("2/7 Downloading/loading LUAD/LUSC counts and metadata")
project_data <- lapply(seq_len(nrow(tcga_projects)), function(i) {
  project_info <- tcga_projects[i, , drop = FALSE]
  project_id <- project_info$project[[1]]
  message(sprintf("Loading %s", project_id))

  project_cache <- recount3_cache(file.path(args[["cache-dir"]], project_id))
  loaded <- tryCatch({
    rse <- create_rse(
      project_info,
      annotation = args[["annotation"]],
      bfc = project_cache,
      recount3_url = args[["recount3-url"]]
    )
    list(
      counts = round(transform_counts(rse)),
      metadata = as.data.table(as.data.frame(colData(rse)), keep.rownames = "sample_id")
    )
  }, error = function(error) {
    warning(sprintf("Using cached recount3 files for %s after create_rse failed: %s", project_id, conditionMessage(error)), call. = FALSE)
    read_cached_tcga_project(args[["cache-dir"]], project_id)
  })

  counts <- loaded$counts
  metadata <- loaded$metadata
  rownames(counts) <- rownames(loaded$counts)
  if (!("sample_id" %in% names(metadata)) || !all(metadata$sample_id %in% colnames(counts))) {
    metadata[, sample_id := colnames(counts)]
  }
  metadata[, project := project_id]

  sex_column <- find_sex_column(metadata)
  metadata[, sex := normalize_sex(get(sex_column))]

  sample_type_column <- find_sample_type_column(metadata)
  if (!identical(tolower(args[["sample-type"]]), "all") && !is.na(sample_type_column)) {
    metadata <- metadata[grepl(args[["sample-type"]], get(sample_type_column), ignore.case = TRUE)]
  }

  metadata <- metadata[!is.na(sex) & sample_id %in% colnames(counts)]
  counts <- counts[, metadata$sample_id, drop = FALSE]

  list(project = project_id, counts = counts, metadata = metadata)
})

# Combine projects on shared gene IDs
common_genes <- Reduce(intersect, lapply(project_data, function(item) rownames(item$counts)))
if (length(common_genes) == 0) {
  stop("No common genes remained across selected TCGA projects.", call. = FALSE)
}

count_matrix <- do.call(cbind, lapply(project_data, function(item) item$counts[common_genes, , drop = FALSE]))
metadata <- rbindlist(lapply(project_data, `[[`, "metadata"), fill = TRUE)
metadata <- metadata[match(colnames(count_matrix), sample_id)]
age_column <- intersect(
  c(
    "tcga.gdc_cases.diagnoses.age_at_diagnosis",
    "gdc_cases.diagnoses.age_at_diagnosis",
    "tcga.cgc_case_age_at_diagnosis",
    "tcga.xml_primary_pathology_age_at_initial_pathologic_diagnosis",
    "tcga.xml_age_at_initial_pathologic_diagnosis"
  ),
  names(metadata)
)[1]
if (is.na(age_column)) {
  stop("TCGA metadata does not contain a recognized age-at-diagnosis column.", call. = FALSE)
}
metadata[, age_at_diagnosis := as.numeric(get(age_column)) / 365.25]
design_columns <- c("sex", "project", "age_at_diagnosis")
metadata <- metadata[complete_design_rows(metadata, design_columns)]
count_matrix <- count_matrix[, metadata$sample_id, drop = FALSE]

if (nrow(metadata) < 6 || length(unique(metadata$sex)) < 2) {
  stop("Need at least two sex groups and enough TCGA samples for the contrast.", call. = FALSE)
}

message(sprintf(
  "Using %s TCGA samples: %s Female, %s Male",
  nrow(metadata),
  sum(metadata$sex == "Female"),
  sum(metadata$sex == "Male")
))

# Run DESeq2 sex contrast with project adjustment
message("3/7 Running DESeq2 Female vs Male contrast")
coldata <- data.frame(
  row.names = metadata$sample_id,
  sex = factor(metadata$sex, levels = c("Male", "Female")),
  project = factor(metadata$project),
  age_at_diagnosis = metadata$age_at_diagnosis
)

design <- if (length(unique(coldata$project)) > 1) ~ project + age_at_diagnosis + sex else ~ age_at_diagnosis + sex
message(sprintf("Using DESeq2 design %s", design_formula))
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

run_summary <- data.table(
  metric = c("samples", "female_samples", "male_samples", "ranked_genes", "design_formula", "projects", "sample_type", "annotation", "min_count"),
  value = c(
    nrow(metadata),
    sum(metadata$sex == "Female"),
    sum(metadata$sex == "Male"),
    nrow(ranks),
    design_formula,
    paste(projects, collapse = ","),
    args[["sample-type"]],
    args[["annotation"]],
    min_count
  )
)
fwrite(run_summary, summary_path, sep = "\t")
}

# Run enrichment for each ALIEN GMT collection
message("5/7 Running multilevel rsfgsea")
rsfgsea_bin <- resolve_rsfgsea(args[["rsfgsea"]])
enrichment_runs <- list()
plot_theme_specs <- list(
  pathways = list(
    theme_spec("immune and antigen presentation", "IMMUN|ANTIGEN|MHC|B_CELL|B CELL|T_CELL|T CELL|LYMPHOID|LEUKOCYTE|CYTOKINE"),
    theme_spec("ribosomal and mitochondrial programs", "RRNA|RIBOSOM|MITOCHON|TRANSLATION|OXIDATIVE|RESPIRATORY_CHAIN|RESPIRATORY CHAIN")
  ),
  cancer_dependency = list(
    theme_spec("lung cancer cell lines", "LUNG|NCIH|NCI_H|LUDLU|RERFLC|DMS454")
  )
)
for (gmt_path in gmt_paths) {
  label <- gmt_label(gmt_path)
  enrichment_path <- file.path(tables_dir, sprintf("tcga_lung_%s_female_vs_male_rsfgsea.tsv", label))
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
  title = "TCGA Lung Cancer Sex Contrast",
  summary_lines = c(
    "TCGA LUAD and LUSC primary tumors from recount3 were contrasted by sex using DESeq2 on raw counts, adjusting for cancer project and age at diagnosis. Genes were ranked by the Wald statistic for Female vs Male and tested against all four ALIEN TCGA GMT collections with multilevel rsfgsea.",
    "",
    sprintf("- Projects: %s.", paste(projects, collapse = ", ")),
    sprintf("- Sample type: %s.", args[["sample-type"]]),
    sprintf(
      "- Samples: %s total (%s Female, %s Male).",
      summary_value(run_summary, "samples"),
      summary_value(run_summary, "female_samples"),
      summary_value(run_summary, "male_samples")
    ),
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
  collection_note = "The figure below gives a collection-level view of significant enrichment. Detailed term-level results are shown in the tables; figures are limited to the clearest pathway and cancer-dependency examples."
)

message(sprintf("Done: %s", report_path))
}


main()
