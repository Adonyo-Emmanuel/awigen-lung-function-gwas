#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript check_v2g_config.R <config.tsv>")
}

script_arg <- commandArgs(FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
source(file.path(script_dir, "v2g_utils.R"))

cfg <- read_config(args[[1]])
cfg_table <- attr(cfg, "table")
out_dir <- pipeline_output_dir(cfg)

if (!"required" %in% names(cfg_table)) {
  cfg_table[, required := "no"]
}

cfg_table[, required := tolower(clean_text(required))]
cfg_table[, is_file_key := grepl("(_file$|_chain$)", key)]
cfg_table[, exists := NA]
cfg_table[is_file_key & !is.na(value), exists := file.exists(value)]
cfg_table[, status := "ok"]
cfg_table[required %in% c("yes", "true", "1") & is.na(value), status := "missing required value"]
cfg_table[is_file_key & !is.na(value) & exists == FALSE, status := "file not found"]

report_file <- file.path(out_dir, "config_check.tsv")
write_tsv(cfg_table, report_file)

failed <- cfg_table[status != "ok"]
if (nrow(failed) > 0) {
  stop(
    "Config check failed. See: ", report_file, "\n",
    paste(failed$key, failed$status, sep = ": ", collapse = "\n")
  )
}

message("Config check passed: ", report_file)
