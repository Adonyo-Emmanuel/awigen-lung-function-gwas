suppressPackageStartupMessages({
  library(data.table)
})

missing_values <- c("", "NA", "N/A", ".", "-", "NULL", "null", "NaN", "nan")

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
