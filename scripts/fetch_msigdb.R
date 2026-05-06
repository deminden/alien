#!/usr/bin/env Rscript

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: fetch_msigdb.R <output_dir> [include_c4_cm] [min_msigdbr_version]")
}
output_dir <- args[[1]]
include_c4_cm <- ifelse(length(args) >= 2, toupper(args[[2]]) %in% c("TRUE", "1", "YES"), FALSE)
min_version <- ifelse(length(args) >= 3, args[[3]], "26.1.0")

# Set output paths
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Check msigdbr version
if (!requireNamespace("msigdbr", quietly = TRUE)) {
  stop("The R package msigdbr is required. Install it with install.packages('msigdbr').")
}
installed_version <- packageVersion("msigdbr")
if (installed_version < package_version(min_version)) {
  stop(sprintf("msigdbr >= %s is required; installed version is %s.", min_version, installed_version))
}

sources <- data.frame(
  source_tag = c("REACTOME", "WIKIPATHWAYS", "KEGG_MEDICUS", "GOBP", "GOMF", "GOCC", "HPO", "C6", "C9", "C4_3CA", "C4_CGN", "C4_CM"),
  collection = c("C2", "C2", "C2", "C5", "C5", "C5", "C5", "C6", "C9", "C4", "C4", "C4"),
  subcollection = c("CP:REACTOME", "CP:WIKIPATHWAYS", "CP:KEGG_MEDICUS", "GO:BP", "GO:MF", "GO:CC", "HPO", "", "", "3CA", "CGN", "CM"),
  stringsAsFactors = FALSE
)
if (!include_c4_cm) {
  sources <- sources[sources$source_tag != "C4_CM", ]
}

fetch_one <- function(collection, subcollection) {
  formals_names <- names(formals(msigdbr::msigdbr))
  base_args <- list(db_species = "HS", species = "Homo sapiens")

  if ("collection" %in% formals_names) {
    base_args$collection <- collection
  } else if ("category" %in% formals_names) {
    base_args$category <- collection
  }

  if (nzchar(subcollection)) {
    if ("subcollection" %in% formals_names) {
      base_args$subcollection <- subcollection
    } else if ("subcategory" %in% formals_names) {
      base_args$subcategory <- subcollection
    }
  }

  result <- tryCatch(do.call(msigdbr::msigdbr, base_args), error = identity)
  if (inherits(result, "error")) {
    stop(sprintf("Could not fetch MSigDB %s %s: %s", collection, subcollection, result$message))
  }
  if (nrow(result) == 0) {
    stop(sprintf("MSigDB %s %s returned no rows.", collection, subcollection))
  }
  result
}

# Fetch MSigDB collections
db_versions <- character()
for (i in seq_len(nrow(sources))) {
  source_tag <- sources$source_tag[[i]]
  message(sprintf("Fetching %s", source_tag))
  tbl <- fetch_one(sources$collection[[i]], sources$subcollection[[i]])
  if ("db_version" %in% names(tbl)) {
    db_versions <- unique(c(db_versions, tbl$db_version))
  }
  output_file <- file.path(output_dir, sprintf("msigdbr_%s.tsv.gz", source_tag))
  utils::write.table(
    tbl,
    gzfile(output_file),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

# Write session info
sink(file.path(output_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

# Write database version
writeLines(unique(db_versions), file.path(output_dir, "db_version.txt"))
