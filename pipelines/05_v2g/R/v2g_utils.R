suppressPackageStartupMessages({
  library(data.table)
})

missing_values <- c("", "NA", "N/A", ".", "-", "NULL", "null", "NaN", "nan")

current_script_file <- function() {
  script_arg <- commandArgs(FALSE)
  file_arg <- script_arg[grepl("^--file=", script_arg)]
  if (length(file_arg) == 0) {
    stop("Could not determine the current script path from commandArgs().")
  }
  sub("^--file=", "", file_arg[[1]])
}

clean_text <- function(x) {
  x <- as.character(x)
  x[trimws(x) %in% missing_values] <- NA_character_
  x
}

num <- function(x) suppressWarnings(as.numeric(x))

strip_ensembl_version <- function(x) {
  x <- clean_text(x)
  sub("\\.[0-9]+$", "", x)
}

normalise_chr <- function(x) {
  x <- clean_text(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x[x == "23"] <- "X"
  x[x == "24"] <- "Y"
  x[x %in% c("M", "25")] <- "MT"
  x
}

read_config <- function(config_file) {
  if (!file.exists(config_file)) {
    stop("Config file does not exist: ", config_file)
  }

  cfg <- fread(config_file, sep = "\t", header = TRUE, fill = TRUE)
  required_cols <- c("key", "value")
  missing_cols <- setdiff(required_cols, names(cfg))
  if (length(missing_cols) > 0) {
    stop("Config file is missing column(s): ", paste(missing_cols, collapse = ", "))
  }

  cfg[, key := clean_text(key)]
  cfg[, value := clean_text(value)]
  cfg <- cfg[!is.na(key)]

  out <- as.list(cfg$value)
  names(out) <- cfg$key
  attr(out, "table") <- cfg
  out
}

config_value <- function(cfg, key, default = NA_character_) {
  value <- cfg[[key]]
  if (is.null(value) || is.na(value) || value == "") default else value
}

pipeline_output_dir <- function(cfg) {
  out_dir <- config_value(cfg, "output_dir")
  if (is.na(out_dir)) stop("Config key output_dir is required.")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(out_dir, mustWork = FALSE)
}

write_tsv <- function(dt, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(dt, path, sep = "\t", quote = FALSE, na = "")
}

first_existing_col <- function(dt, candidates) {
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

pull_col <- function(dt, candidates, default = NA_character_) {
  col <- first_existing_col(dt, candidates)
  if (is.na(col)) rep(default, nrow(dt)) else dt[[col]]
}

pull_first_non_missing <- function(dt, candidates, default = NA_character_) {
  cols <- candidates[candidates %in% names(dt)]
  if (length(cols) == 0) return(rep(default, nrow(dt)))

  out <- rep(default, nrow(dt))
  for (col in cols) {
    value <- clean_text(dt[[col]])
    fill <- is.na(clean_text(out)) & !is.na(value)
    out[fill] <- value[fill]
  }
  out
}

make_variant_key <- function(chr, pos, a1, a2) {
  chr <- normalise_chr(chr)
  pos <- suppressWarnings(as.integer(pos))
  a1 <- toupper(clean_text(a1))
  a2 <- toupper(clean_text(a2))
  out <- rep(NA_character_, length(chr))
  ok <- !is.na(chr) & !is.na(pos) & !is.na(a1) & !is.na(a2)
  out[ok] <- paste(chr[ok], pos[ok], pmin(a1[ok], a2[ok]), pmax(a1[ok], a2[ok]), sep = "_")
  out
}

source_file_md5 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

file_manifest_row <- function(label, path) {
  exists <- file.exists(path)
  info <- if (exists) file.info(path) else NULL
  data.table(
    input_label = label,
    path = path,
    exists = exists,
    size_bytes = if (exists) info$size else NA_real_,
    modified_time = if (exists) as.character(info$mtime) else NA_character_,
    md5 = source_file_md5(path)
  )
}

ordered_unique <- function(x) {
  x <- clean_text(x)
  unique(x[!is.na(x)])
}

join_unique <- function(x, sep = ";") {
  x <- ordered_unique(x)
  if (length(x) == 0) NA_character_ else paste(sort(x), collapse = sep)
}

minimum_numeric <- function(x) {
  x <- num(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else min(x)
}

maximum_numeric <- function(x) {
  x <- num(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else max(x)
}

standard_evidence_columns <- function() {
  c(
    "evidence_line",
    "evidence_source",
    "evidence_class",
    "source_file",
    "sentinel_snp",
    "variant_id",
    "rsid",
    "chromosome_b37",
    "position_b37",
    "chromosome_b38",
    "position_b38",
    "effect_allele",
    "other_allele",
    "trait",
    "traits",
    "gwas_p_value",
    "gwas_z_score",
    "gwas_beta",
    "gwas_se",
    "posterior_probability",
    "gene_id",
    "gene_symbol",
    "gene_or_protein",
    "molecular_trait_type",
    "molecular_trait_id",
    "tissue_or_protein_panel",
    "p_value",
    "fdr",
    "q_value",
    "threshold_p_value",
    "effect_estimate",
    "standard_error",
    "z_score",
    "distance_from_sentinel_bp",
    "consequence",
    "impact",
    "sift",
    "polyphen",
    "cadd_phred",
    "evidence_detail"
  )
}

numeric_evidence_columns <- function() {
  c(
    "position_b37",
    "position_b38",
    "gwas_p_value",
    "gwas_z_score",
    "gwas_beta",
    "gwas_se",
    "posterior_probability",
    "p_value",
    "fdr",
    "q_value",
    "threshold_p_value",
    "effect_estimate",
    "standard_error",
    "z_score",
    "distance_from_sentinel_bp",
    "cadd_phred"
  )
}

add_missing_standard_cols <- function(dt) {
  for (col in setdiff(standard_evidence_columns(), names(dt))) {
    dt[, (col) := NA]
  }
  setcolorder(dt, standard_evidence_columns())
  invisible(dt)
}

standardise_source_table <- function(dt, mapping, constants = list(), source_file) {
  standard_cols <- standard_evidence_columns()
  source_file_name <- basename(source_file)
  out <- as.data.table(
    setNames(
      replicate(length(standard_cols), rep(NA_character_, nrow(dt)), simplify = FALSE),
      standard_cols
    )
  )

  for (col in standard_cols) {
    if (!is.null(constants[[col]])) {
      out[, (col) := constants[[col]]]
    } else if (!is.null(mapping[[col]])) {
      out[, (col) := pull_first_non_missing(dt, mapping[[col]])]
    }
  }

  out[, source_file := source_file_name]
  out[, gene_id := strip_ensembl_version(gene_id)]
  out[, gene_symbol := clean_text(gene_symbol)]
  out[, gene_or_protein := clean_text(gene_or_protein)]
  out[is.na(gene_or_protein) & !is.na(gene_symbol), gene_or_protein := gene_symbol]
  out[is.na(gene_or_protein) & !is.na(gene_id), gene_or_protein := gene_id]
  out[is.na(variant_id) & !is.na(sentinel_snp), variant_id := sentinel_snp]

  for (col in intersect(numeric_evidence_columns(), names(out))) {
    out[, (col) := num(get(col))]
  }

  out <- out[!is.na(sentinel_snp) & !is.na(gene_or_protein)]
  add_missing_standard_cols(out)
  out
}

standard_evidence_dir <- function(cfg) {
  out_dir <- file.path(pipeline_output_dir(cfg), "standardised_line_evidence")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir
}

write_standard_evidence <- function(dt, cfg, file_name) {
  add_missing_standard_cols(dt)
  path <- file.path(standard_evidence_dir(cfg), file_name)
  write_tsv(dt, path)
  message("Wrote standardised evidence: ", path)
  invisible(path)
}

read_required_table <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    stop("Missing required file: ", path)
  }
  fread(path, check.names = FALSE)
}
