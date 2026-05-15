# Shared helpers for sex-contrast enrichment reports
# Source this file from GTEx and TCGA workflow scripts

# Command-line parsing
parse_args <- function(defaults) {
  raw <- commandArgs(trailingOnly = TRUE)
  parsed <- defaults
  i <- 1
  while (i <= length(raw)) {
    if (!startsWith(raw[[i]], "--")) {
      stop(sprintf("Unknown positional argument: %s", raw[[i]]), call. = FALSE)
    }
    if (grepl("^--[^=]+=", raw[[i]])) {
      key <- sub("^--", "", sub("=.*$", "", raw[[i]]))
      parsed[[key]] <- sub("^--[^=]+=", "", raw[[i]])
      i <- i + 1
    } else {
      key <- sub("^--", "", raw[[i]])
      if (i == length(raw) || startsWith(raw[[i + 1]], "--")) {
        parsed[[key]] <- TRUE
        i <- i + 1
      } else {
        parsed[[key]] <- raw[[i + 1]]
        i <- i + 2
      }
    }
  }
  parsed
}

as_integer <- function(value, label, minimum = 1) {
  value <- suppressWarnings(as.integer(value))
  if (is.na(value) || value < minimum) {
    stop(sprintf("%s must be an integer >= %s.", label, minimum), call. = FALSE)
  }
  value
}

summary_value <- function(summary_table, key) {
  value <- summary_table[["value"]][summary_table[["metric"]] == key][1]
  if (is.na(value)) "unknown" else as.character(value)
}

summary_matches <- function(path, expected) {
  if (!file.exists(path)) return(FALSE)
  summary_table <- fread(path)
  all(vapply(names(expected), function(key) {
    identical(summary_value(summary_table, key), as.character(expected[[key]]))
  }, logical(1), USE.NAMES = FALSE))
}

complete_design_rows <- function(table, columns) {
  Reduce(`&`, lapply(columns, function(col) {
    value <- table[[col]]
    !is.na(value) & nzchar(trimws(as.character(value)))
  }))
}

assert_full_rank_design <- function(design, coldata) {
  model_matrix <- model.matrix(design, coldata)
  if (qr(model_matrix)$rank < ncol(model_matrix)) {
    stop("DESeq2 design matrix is not full rank for the selected samples.", call. = FALSE)
  }
  invisible(TRUE)
}

download_if_missing <- function(url, path) {
  if (!file.exists(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    download.file(url, destfile = path, mode = "wb", quiet = FALSE)
  }
}

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s does not exist: %s", label, path), call. = FALSE)
  }
}

ensure_gmt_files <- function(paths, build_configs, workers, alien_bin) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) == 0) return(invisible(TRUE))
  if (length(build_configs) == 0) {
    stop(sprintf("GMT file(s) missing: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }

  message("Missing GMT file(s), building ALIEN example config(s):")
  for (path in missing) message(sprintf("  - %s", path))

  for (config_path in build_configs) {
    stop_if_missing(config_path, "ALIEN config")
    message(sprintf("Running alien build --config %s --workers %s", config_path, workers))
    status <- tryCatch(
      system2(alien_bin, args = c("build", "--config", config_path, "--workers", workers)),
      warning = function(warning) {
        message(conditionMessage(warning))
        127L
      },
      error = function(error) {
        message(conditionMessage(error))
        127L
      }
    )
    if (!identical(status, 0L)) {
      stop(sprintf("alien build failed for %s with exit status %s.", config_path, status), call. = FALSE)
    }
  }

  still_missing <- paths[!file.exists(paths)]
  if (length(still_missing) > 0) {
    stop(sprintf("GMT file(s) still missing after alien build: %s", paste(still_missing, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

ensure_recount3_g029_output_genes <- function(
  gtf_path = "data/recount3/human.gene_sums.G029.gtf.gz",
  output_genes_path = "data/recount3/tcga_gencode_v29_output_genes.tsv.gz",
  gtf_url = "https://recount-opendata.s3.amazonaws.com/recount3/release/human/annotations/gene_sums/human.gene_sums.G029.gtf.gz"
) {
  download_if_missing(gtf_url, gtf_path)
  if (file.exists(output_genes_path)) return(invisible(TRUE))

  message(sprintf("Writing recount3 G029 output-gene list: %s", output_genes_path))
  con <- gzfile(gtf_path, open = "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  lines <- lines[!startsWith(lines, "#")]
  with_gene_id <- grepl('gene_id "', lines, fixed = TRUE)
  gene_ids <- sub('.*gene_id "([^"]+)".*', "\\1", lines[with_gene_id])
  gene_ids <- unique(gene_ids[nzchar(gene_ids)])
  if (length(gene_ids) == 0) {
    stop(sprintf("Could not extract gene_id attributes from %s.", gtf_path), call. = FALSE)
  }
  dir.create(dirname(output_genes_path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(data.table::data.table(Ensembl_gene_ID = gene_ids), output_genes_path, sep = "\t", compress = "gzip")
  invisible(TRUE)
}

parse_paths <- function(value) {
  paths <- trimws(unlist(strsplit(value, "[,;]")))
  paths[nzchar(paths)]
}

common_path_prefix <- function(paths) {
  if (length(paths) == 0) return("")
  split_paths <- strsplit(paths, .Platform$file.sep, fixed = TRUE)
  common <- split_paths[[1]]
  for (parts in split_paths[-1]) {
    limit <- min(length(common), length(parts))
    keep <- 0L
    while (keep < limit && identical(common[[keep + 1L]], parts[[keep + 1L]])) {
      keep <- keep + 1L
    }
    common <- if (keep > 0) common[seq_len(keep)] else character()
  }
  paste(common, collapse = .Platform$file.sep)
}

format_gmt_summary <- function(paths) {
  if (length(paths) == 1) {
    return(sprintf("- GMT file: `%s`.", paths[[1]]))
  }
  base <- common_path_prefix(paths)
  if (!nzchar(base) || identical(base, ".")) {
    return(sprintf("- GMT files: `%s`.", paste(paths, collapse = "`, `")))
  }
  rel <- vapply(paths, function(path) {
    sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", base), .Platform$file.sep), "", path)
  }, character(1), USE.NAMES = FALSE)
  sprintf("- GMT files under `%s`: `%s`.", base, paste(rel, collapse = "`, `"))
}

# Gene-set and term label helpers
gmt_label <- function(path) {
  parts <- strsplit(normalizePath(path, mustWork = FALSE), .Platform$file.sep, fixed = TRUE)[[1]]
  gmt_index <- which(parts == "gmt")
  if (length(gmt_index) > 0 && gmt_index[[1]] > 1) {
    return(safe_name(parts[[gmt_index[[1]] - 1]]))
  }
  safe_name(tools::file_path_sans_ext(basename(path)))
}

pretty_term <- function(term) {
  parts <- split_term_id(term)
  paste(parts$collection_name, parts$pathway, sep = ": ")
}

plot_term_title <- function(term, width = 105) {
  parts <- split_term_id(term)[1]
  collection <- parts$collection_name
  if (!is.na(parts$source_type) && nzchar(parts$source_type)) {
    collection <- paste(c(parts$source_type, collection)[!is.na(c(parts$source_type, collection))], collapse = " ")
  }
  pieces <- c(collection, parts$pathway)
  pieces <- pieces[!is.na(pieces) & nzchar(pieces)]
  shorten_label(paste(pieces, collapse = ": "), width = width)
}

display_run_label <- function(label) {
  labels <- c(
    "pathways" = "Pathways",
    "disease_phenotype" = "Disease/phenotype",
    "function_location" = "Function/location",
    "cancer_dependency" = "Cancer dependency"
  )
  vapply(as.character(label), function(one) {
    if (one %in% names(labels)) labels[[one]] else format_phrase(one, title_case = TRUE)
  }, character(1), USE.NAMES = FALSE)
}

format_collection_name <- function(value) {
  known <- c(
    "PATHWAYS" = "Pathways",
    "DISEASE PHENOTYPE" = "Disease/phenotype",
    "DISEASE/PHENOTYPE" = "Disease/phenotype",
    "FUNCTION LOCATION" = "Function/location",
    "FUNCTION/LOCATION" = "Function/location",
    "CANCER DEPENDENCY" = "Cancer dependency",
    "BIOCARTA" = "BioCarta",
    "CANCER CELL LINE ENCYCLOPEDIA" = "Cancer Cell Line Encyclopedia",
    "CGN" = "CGN",
    "DISGENET" = "DisGeNET",
    "GOBP" = "GOBP",
    "GOCC" = "GOCC",
    "GOMF" = "GOMF",
    "GWAS CATALOG 2025" = "GWAS Catalog 2025",
    "HALLMARK" = "Hallmark",
    "HPO" = "HPO",
    "JENSEN DISEASES CURATED 2025" = "Jensen Diseases Curated 2025",
    "KEGG" = "KEGG",
    "KEGG MEDICUS" = "KEGG Medicus",
    "MIR" = "MIR",
    "PID" = "PID",
    "REACTOME" = "Reactome",
    "TFT" = "TFT",
    "WIKIPATHWAYS" = "WikiPathways"
  )
  clean <- gsub("_", " ", as.character(value))
  key <- toupper(trimws(clean))
  vapply(seq_along(clean), function(i) {
    if (is.na(clean[[i]])) return(NA_character_)
    if (key[[i]] %in% names(known)) known[[key[[i]]]] else format_phrase(clean[[i]], title_case = TRUE)
  }, character(1), USE.NAMES = FALSE)
}

format_term_name <- function(value) {
  clean <- gsub("_", " ", as.character(value))
  clean <- sub("^HP\\s+", "", clean)
  clean <- sub("^WP\\s+", "", clean)
  format_phrase(clean, title_case = FALSE)
}

format_phrase <- function(value, title_case = FALSE) {
  value <- gsub("_", " ", as.character(value))
  missing <- is.na(value)
  value <- gsub("\\s+", " ", trimws(value), perl = TRUE)
  if (title_case) {
    value <- tools::toTitleCase(tolower(value))
  } else {
    value <- tolower(value)
    value <- ifelse(nzchar(value), paste0(toupper(substr(value, 1, 1)), substr(value, 2, nchar(value))), value)
  }
  value <- restore_alphanumeric_tokens(value)
  value <- restore_known_acronyms(value)
  value[missing] <- NA_character_
  value
}

restore_alphanumeric_tokens <- function(value) {
  for (i in seq_along(value)) {
    if (is.na(value[[i]])) next
    matches <- gregexpr("\\b(?=[A-Za-z0-9-]*[A-Za-z])(?=[A-Za-z0-9-]*[0-9])[A-Za-z0-9-]+\\b", value[[i]], perl = TRUE)
    regmatches(value[[i]], matches) <- lapply(regmatches(value[[i]], matches), toupper)
  }
  value
}

restore_known_acronyms <- function(value) {
  replacements <- c(
    "abc" = "ABC",
    "adp" = "ADP",
    "akt" = "AKT",
    "amp" = "AMP",
    "atp" = "ATP",
    "bcr" = "BCR",
    "b cell" = "B cell",
    "brca1" = "BRCA1",
    "brca2" = "BRCA2",
    "celllines" = "cell lines",
    "covid-19" = "COVID-19",
    "crispr" = "CRISPR",
    "depmap" = "DepMap",
    "dna" = "DNA",
    "egfr" = "EGFR",
    "erk" = "ERK",
    "esc" = "ESC",
    "gdp" = "GDP",
    "genedependency" = "gene dependency",
    "gsu" = "GSU",
    "gtp" = "GTP",
    "hdacs" = "HDACs",
    "hla" = "HLA",
    "ifn" = "IFN",
    "jak" = "JAK",
    "lncrna" = "lncRNA",
    "mapk" = "MAPK",
    "mhc" = "MHC",
    "mirna" = "miRNA",
    "morf" = "MORF",
    "mrna" = "mRNA",
    "mtor" = "mTOR",
    "nadh" = "NADH",
    "nadph" = "NADPH",
    "nf-kb" = "NF-kB",
    "pi3k" = "PI3K",
    "rna" = "RNA",
    "rrna" = "rRNA",
    "sars-cov-2" = "SARS-CoV-2",
    "sarscov2" = "SARS-CoV-2",
    "stat" = "STAT",
    "t cell" = "T cell",
    "tca" = "TCA",
    "tcr" = "TCR",
    "tcf" = "TCF",
    "tgf" = "TGF",
    "tnf" = "TNF",
    "trna" = "tRNA",
    "vegf" = "VEGF",
    "vegf a" = "VEGFA",
    "vegfa" = "VEGFA",
    "wnt" = "WNT"
  )
  for (key in names(replacements)) {
    value <- gsub(paste0("\\b", key, "\\b"), replacements[[key]], value, ignore.case = TRUE, perl = TRUE)
  }
  value
}

split_term_id <- function(term) {
  rbindlist(lapply(term, function(one) {
    if (is.na(one)) {
      return(data.table(source_type = NA_character_, collection_name = NA_character_, pathway = NA_character_))
    }
    if (!grepl("__", one, fixed = TRUE)) {
      return(data.table(source_type = NA_character_, collection_name = NA_character_, pathway = format_term_name(one)))
    }
    left <- sub("__.*$", "", one)
    right <- sub("^.*?__", "", one)
    source_type <- fifelse(startsWith(left, "MSIGDB_"), "MSigDB", fifelse(startsWith(left, "ENRICHR_"), "Enrichr", "Other"))
    collection <- sub("^MSIGDB_", "", sub("^ENRICHR_", "", left))
    pathway <- right
    pathway <- sub(paste0("^", gsub("([\\W])", "\\\\\\1", collection), "_"), "", pathway)
    data.table(
      source_type = source_type,
      collection_name = format_collection_name(collection),
      pathway = format_term_name(pathway)
    )
  }), fill = TRUE)
}

shorten_label <- function(label, width = 90) {
  ifelse(nchar(label) > width, paste0(substr(label, 1, width - 3), "..."), label)
}

stop_if_missing_columns <- function(table, columns, label) {
  missing <- setdiff(columns, names(table))
  if (length(missing) > 0) {
    stop(sprintf("%s is missing required columns: %s", label, paste(missing, collapse = ", ")), call. = FALSE)
  }
}

# Metadata helpers
normalize_sex <- function(values) {
  values <- trimws(tolower(as.character(values)))
  fifelse(
    values %in% c("female", "f", "2", "xx"),
    "Female",
    fifelse(values %in% c("male", "m", "1", "xy"), "Male", NA_character_)
  )
}

find_sex_column <- function(table) {
  candidates <- grep("gender|sex", names(table), ignore.case = TRUE, value = TRUE)
  if (length(candidates) == 0) {
    stop("Could not find a sex/gender column in the metadata.", call. = FALSE)
  }
  candidates[order(grepl("gender", candidates, ignore.case = TRUE), nchar(candidates))][[1]]
}

find_sample_type_column <- function(table) {
  candidates <- grep("sample.*type|sample_type", names(table), ignore.case = TRUE, value = TRUE)
  if (length(candidates) == 0) return(NA_character_)
  candidates[[1]]
}

strip_ensembl_version <- function(ids) {
  sub("\\.[0-9]+$", "", ids)
}

# rsfgsea helpers
resolve_rsfgsea <- function(path_or_name) {
  if (file.exists(path_or_name) && file.access(path_or_name, mode = 1) == 0) {
    return(normalizePath(path_or_name))
  }
  located <- Sys.which(path_or_name)
  if (nzchar(located)) {
    return(unname(located))
  }
  cargo_bin <- file.path(Sys.getenv("HOME"), ".cargo", "bin", "rsfgsea")
  if (identical(path_or_name, "rsfgsea") && file.exists(cargo_bin) && file.access(cargo_bin, mode = 1) == 0) {
    return(normalizePath(cargo_bin))
  }
  stop("Could not find rsfgsea. Put it on PATH, set RSFGSEA_BIN, or pass --rsfgsea /path/to/rsfgsea.", call. = FALSE)
}

run_rsfgsea <- function(rsfgsea_bin, ranks_path, gmt_path, output_path, min_size, max_size, nperm_simple, eps, seed, workers) {
  status <- system2(
    rsfgsea_bin,
    args = c(
      "--ranks", ranks_path,
      "--gmt", gmt_path,
      "--output", output_path,
      "--mode", "multilevel",
      "--minSize", min_size,
      "--maxSize", max_size,
      "--nPermSimple", nperm_simple,
      "--eps", eps,
      "--seed", seed,
      "--nproc", workers
    )
  )
  if (!identical(status, 0L)) {
    stop(sprintf("rsfgsea failed with exit status %s.", status), call. = FALSE)
  }
}

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  terms <- list()
  for (line in lines) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) >= 3) {
      terms[[fields[[1]]]] <- unique(fields[-c(1, 2)])
    }
  }
  terms
}

read_enrichment_table <- function(path) {
  table <- fread(path)
  setnames(table, names(table), tolower(names(table)))
  if (!("pathway" %in% names(table))) {
    stop("rsfgsea output did not contain a pathway column.", call. = FALSE)
  }
  if ("nes" %in% names(table)) table[, nes := as.numeric(nes)]
  if ("padj" %in% names(table)) table[, padj := as.numeric(padj)]
  table
}

theme_spec <- function(name, pattern, direction = NULL) {
  list(name = name, pattern = pattern, direction = direction)
}

plot_representative_running_scores <- function(enrichment, ranks, gmt_terms, plots_dir, prefix, theme_specs = list()) {
  chosen <- select_biological_terms(enrichment, theme_specs, max_terms = max(length(theme_specs), 2))
  if (nrow(chosen) == 0) {
    candidates <- copy(enrichment)[is.finite(nes)]
    if (nrow(candidates) == 0) return(character())
    candidates[, padj_rank := if ("padj" %in% names(candidates)) fifelse(is.na(padj), 1, padj) else 1]
    chosen <- rbind(
      candidates[nes > 0][order(padj_rank, -nes)][1],
      candidates[nes < 0][order(padj_rank, nes)][1],
      fill = TRUE
    )
    chosen[, theme := "strongest available term"]
  }
  chosen <- chosen[!is.na(pathway)]

  paths <- character()
  used_names <- character()
  for (i in seq_len(nrow(chosen))) {
    term <- chosen$pathway[[i]]
    genes <- gmt_terms[[term]]
    if (is.null(genes)) next
    stem <- plot_stem(prefix, chosen$theme[[i]], "running_score")
    stem <- unique_plot_stem(stem, used_names)
    used_names <- c(used_names, stem)
    plot_path <- file.path(plots_dir, paste0(stem, ".png"))
    if (plot_running_score(ranks, genes, term, plot_path)) {
      paths <- c(paths, plot_path)
    }
  }
  paths
}

select_biological_terms <- function(enrichment, theme_specs, max_terms) {
  if (length(theme_specs) == 0 || !all(c("pathway", "nes") %in% names(enrichment))) return(data.table())
  candidates <- copy(enrichment)[is.finite(nes) & !is.na(pathway)]
  if (nrow(candidates) == 0) return(data.table())
  candidates[, padj_rank := if ("padj" %in% names(candidates)) fifelse(is.na(padj), 1, padj) else 1]
  candidates[, abs_nes := abs(nes)]

  selected <- list()
  used <- character()
  for (spec in theme_specs) {
    matches <- candidates[grepl(spec$pattern, pathway, ignore.case = TRUE)]
    if (!is.null(spec$direction)) {
      matches <- if (identical(spec$direction, "positive")) matches[nes > 0] else matches[nes < 0]
    }
    matches <- matches[!(pathway %in% used)]
    if (nrow(matches) == 0) next
    matches <- matches[order(padj_rank, -abs_nes)]
    chosen <- matches[1]
    chosen[, theme := spec$name]
    selected[[length(selected) + 1L]] <- chosen
    used <- c(used, chosen$pathway)
    if (length(selected) >= max_terms) break
  }
  if (length(selected) == 0) return(data.table())
  rbindlist(selected, fill = TRUE)
}

cleanup_running_score_plots <- function(plots_dir, prefix) {
  old <- Sys.glob(file.path(plots_dir, paste0(prefix, "_*_running_score.png")))
  if (length(old) > 0) unlink(old)
}

plot_running_score <- function(ranks, genes, term, path) {
  ranks <- ranks[order(-score)]
  hits <- ranks$gene %in% genes
  n_hits <- sum(hits)
  n_total <- nrow(ranks)
  if (n_hits == 0 || n_hits == n_total) return(FALSE)
  weights <- abs(ranks$score)
  running <- cumsum(ifelse(hits, weights / sum(weights[hits]), -1 / (n_total - n_hits)))
  curve <- data.table(position = seq_len(n_total), running = running, hit = hits)

  plot <- ggplot(curve, aes(x = position, y = running)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_line(color = "#2f5aa8", linewidth = 0.7) +
    geom_rug(data = curve[hit == TRUE], aes(x = position), sides = "b", inherit.aes = FALSE, alpha = 0.25) +
    labs(title = plot_term_title(term), x = "Ranked genes", y = "Running enrichment score") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 11.5, lineheight = 1.05, margin = margin(b = 7))
    )

  ggsave(path, plot, width = 8, height = 4.5, dpi = 150)
  TRUE
}

# Plot collection-level significant term counts
plot_collection_significance <- function(enrichment_runs, path) {
  rows <- rbindlist(lapply(names(enrichment_runs), function(label) {
    enrichment <- copy(enrichment_runs[[label]]$enrichment)
    if (!all(c("nes", "padj") %in% names(enrichment))) return(NULL)
    data.table(
      collection = label,
      direction = c("Female > Male", "Male > Female"),
      terms = c(
        sum(enrichment$padj <= 0.05 & enrichment$nes > 0, na.rm = TRUE),
        sum(enrichment$padj <= 0.05 & enrichment$nes < 0, na.rm = TRUE)
      )
    )
  }), fill = TRUE)
  if (nrow(rows) == 0) return(invisible(FALSE))
  rows[, collection_label := display_run_label(collection)]
  rows[, collection_label := factor(collection_label, levels = display_run_label(names(enrichment_runs)))]

  plot <- ggplot(rows, aes(x = collection_label, y = terms, fill = direction)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.64) +
    scale_fill_manual(values = c("Female > Male" = "#b04a7a", "Male > Female" = "#2f5aa8")) +
    labs(x = NULL, y = "FDR <= 0.05 terms", fill = NULL) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "top", panel.grid.major.x = element_blank())

  ggsave(path, plot, width = 7.2, height = 4.2, dpi = 150)
  invisible(TRUE)
}

safe_name <- function(label) {
  substr(gsub("[^A-Za-z0-9_\\-]+", "_", label), 1, 90)
}

plot_slug <- function(label, max_chars = 90) {
  slug <- tolower(gsub("[^A-Za-z0-9]+", "_", as.character(label)))
  slug <- gsub("_+", "_", slug)
  slug <- gsub("^_|_$", "", slug)
  slug <- ifelse(nzchar(slug), slug, "plot")
  substr(slug, 1, max_chars)
}

plot_stem <- function(..., max_chars = 120) {
  parts <- unlist(list(...), use.names = FALSE)
  parts <- parts[!is.na(parts) & nzchar(as.character(parts))]
  stem <- paste(vapply(parts, plot_slug, character(1), max_chars = 55), collapse = "_")
  substr(stem, 1, max_chars)
}

unique_plot_stem <- function(stem, used) {
  if (!(stem %in% used)) return(stem)
  i <- 2L
  candidate <- paste0(stem, "_", i)
  while (candidate %in% used) {
    i <- i + 1L
    candidate <- paste0(stem, "_", i)
  }
  candidate
}

# Relative path helper for Markdown links
relative_path <- function(path, base_dir) {
  target <- normalizePath(path, mustWork = FALSE)
  base <- normalizePath(base_dir, mustWork = FALSE)
  target_parts <- strsplit(target, .Platform$file.sep, fixed = TRUE)[[1]]
  base_parts <- strsplit(base, .Platform$file.sep, fixed = TRUE)[[1]]
  common <- 0L
  limit <- min(length(target_parts), length(base_parts))
  while (common < limit && identical(target_parts[[common + 1L]], base_parts[[common + 1L]])) {
    common <- common + 1L
  }
  if (common == 0L) return(target)
  up <- rep("..", length(base_parts) - common)
  down <- if (common < length(target_parts)) target_parts[seq.int(common + 1L, length(target_parts))] else character()
  pieces <- c(up, down)
  if (length(pieces) == 0) return(".")
  paste(pieces, collapse = "/")
}

# Write report body
write_report <- function(path, title, summary_lines, enrichment_runs, deseq_path, ranks_path, collection_plot_path, outdir, collection_note) {
  relative <- function(target) relative_path(target, outdir)
  body <- c(
    paste0("# ", title),
    "",
    summary_lines,
    "",
    "## Main Result",
    "",
    overall_result_summary(enrichment_runs),
    "",
    "## Data Products",
    "",
    sprintf("- DESeq2 table: [%s](%s)", basename(deseq_path), relative(deseq_path)),
    sprintf("- Rank file: [%s](%s)", basename(ranks_path), relative(ranks_path)),
    "",
    "## Collection Summary",
    "",
    collection_note,
    "",
    sprintf("![%s](%s)", basename(collection_plot_path), relative(collection_plot_path)),
    "",
    collection_overview(enrichment_runs),
    "",
    "## Enrichment Results",
    ""
  )

  for (label in names(enrichment_runs)) {
    run <- enrichment_runs[[label]]
    enrichment_sorted <- copy(run$enrichment)
    if (all(c("padj", "nes") %in% names(enrichment_sorted))) {
      enrichment_sorted <- enrichment_sorted[order(padj, -abs(nes))]
    }
    table_columns <- intersect(c("pathway", "nes", "pval", "padj", "size"), names(enrichment_sorted))
    body <- c(
      body,
      sprintf("### %s", display_run_label(label)),
      "",
      sprintf("- GMT: `%s`.", run$gmt),
      sprintf("- rsfgsea table: [%s](%s)", basename(run$enrichment_path), relative(run$enrichment_path)),
      describe_enrichment(run$enrichment),
      ""
    )
    selected_figures <- run$figure_paths[file.exists(run$figure_paths)]
    body <- c(body, "Top enrichment rows:", "", markdown_table(enrichment_sorted, table_columns, n = 8), "")
    if (length(selected_figures) > 0) {
      figure_heading <- if (length(selected_figures) == 1) "Selected figure:" else "Selected figures:"
      body <- c(body, figure_heading, "")
      for (figure in selected_figures) {
        body <- c(body, sprintf("![%s](%s)", basename(figure), relative(figure)), "")
      }
    }
  }
  writeLines(body, path)
}

collection_overview <- function(enrichment_runs) {
  rows <- rbindlist(lapply(names(enrichment_runs), function(label) {
    enrichment <- enrichment_runs[[label]]$enrichment
    data.table(
      collection = label,
      tested_terms = nrow(enrichment),
      fdr_0_05 = if ("padj" %in% names(enrichment)) sum(enrichment$padj <= 0.05, na.rm = TRUE) else NA_integer_,
      strongest_female = describe_direction(enrichment, "positive"),
      strongest_male = describe_direction(enrichment, "negative")
    )
  }), fill = TRUE)
  markdown_table(rows, names(rows), n = nrow(rows))
}

overall_result_summary <- function(enrichment_runs) {
  rows <- rbindlist(lapply(names(enrichment_runs), function(label) {
    enrichment <- enrichment_runs[[label]]$enrichment
    if (!all(c("nes", "padj") %in% names(enrichment))) return(NULL)
    data.table(
      collection = label,
      female_terms = sum(enrichment$padj <= 0.05 & enrichment$nes > 0, na.rm = TRUE),
      male_terms = sum(enrichment$padj <= 0.05 & enrichment$nes < 0, na.rm = TRUE)
    )
  }), fill = TRUE)
  if (nrow(rows) == 0 || sum(rows$female_terms + rows$male_terms) == 0) {
    return("No collection produced FDR-significant enrichment at 0.05. The report still shows the strongest nominal trends for review.")
  }
  strongest <- rows[order(-(female_terms + male_terms))][1]
  dominant <- if (strongest$female_terms > strongest$male_terms) "Female > Male" else "Male > Female"
  sprintf(
    "FDR-significant enrichment was present in %s of %s collections. The broadest signal was in `%s` (%s Female > Male terms, %s Male > Female terms), with the larger directional count in %s sets.",
    sum(rows$female_terms + rows$male_terms > 0),
    nrow(rows),
    display_run_label(strongest$collection),
    strongest$female_terms,
    strongest$male_terms,
    dominant
  )
}

describe_enrichment <- function(enrichment) {
  total <- nrow(enrichment)
  significant <- if ("padj" %in% names(enrichment)) sum(enrichment$padj <= 0.05, na.rm = TRUE) else NA_integer_
  c(
    sprintf("- Tested %s terms; %s had FDR <= 0.05.", total, ifelse(is.na(significant), "unknown number", significant)),
    sprintf("- Strongest Female > Male signal: %s.", describe_direction(enrichment, "positive")),
    sprintf("- Strongest Male > Female signal: %s.", describe_direction(enrichment, "negative"))
  )
}

describe_direction <- function(enrichment, direction) {
  table <- copy(enrichment)[is.finite(nes)]
  table <- if (identical(direction, "positive")) table[nes > 0] else table[nes < 0]
  if (nrow(table) == 0) return("none")
  table[, padj_rank := if ("padj" %in% names(table)) fifelse(is.na(padj), 1, padj) else 1]
  table <- if (identical(direction, "positive")) table[order(padj_rank, -nes)] else table[order(padj_rank, nes)]
  best <- table[1]
  fdr <- if ("padj" %in% names(best)) sprintf(", FDR %s", signif(best$padj, 3)) else ""
  sprintf("%s (NES %s%s)", shorten_label(pretty_term(best$pathway), 80), signif(best$nes, 3), fdr)
}

markdown_table <- function(table, columns, n = 12) {
  if (nrow(table) == 0 || length(columns) == 0) return("_No enrichment rows were produced._")
  shown <- copy(table)[seq_len(min(nrow(table), n)), columns, with = FALSE]
  if ("pathway" %in% names(shown)) {
    term_parts <- split_term_id(shown$pathway)
    shown[, pathway := NULL]
    shown <- cbind(term_parts, shown)
    setnames(shown, c("source_type", "collection_name"), c("source", "collection"))
  }
  if ("collection" %in% names(shown)) {
    shown[, collection := format_collection_name(collection)]
  }
  if ("term" %in% names(shown)) {
    shown[, term := format_term_name(term)]
  }
  for (column in names(shown)) {
    if (is.numeric(shown[[column]])) {
      values <- shown[[column]]
      if (all(is.na(values) | values == round(values))) {
        shown[[column]] <- as.character(as.integer(values))
      } else {
        shown[[column]] <- as.character(signif(values, 4))
      }
    }
    shown[[column]] <- gsub("\\|", "\\\\|", as.character(shown[[column]]))
  }
  paste(
    c(
      paste(markdown_column_labels(names(shown)), collapse = " | "),
      paste(rep("---", ncol(shown)), collapse = " | "),
      apply(shown, 1, paste, collapse = " | ")
    ),
    collapse = "\n"
  )
}

markdown_column_labels <- function(columns) {
  labels <- c(
    "collection" = "Collection",
    "direction" = "Direction",
    "fdr_0_05" = "FDR <= 0.05",
    "nes" = "NES",
    "padj" = "FDR",
    "pathway" = "Term",
    "pval" = "P value",
    "size" = "Genes",
    "source" = "Source",
    "strongest_female" = "Strongest Female > Male",
    "strongest_male" = "Strongest Male > Female",
    "term" = "Term",
    "tested_terms" = "Tested terms"
  )
  vapply(columns, function(column) {
    if (column %in% names(labels)) labels[[column]] else format_phrase(column, title_case = TRUE)
  }, character(1), USE.NAMES = FALSE)
}
