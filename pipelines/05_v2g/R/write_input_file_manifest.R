#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript write_input_file_manifest.R <config.tsv>")
}

script_arg <- commandArgs(FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
source(file.path(script_dir, "v2g_utils.R"))

cfg <- read_config(args[[1]])
cfg_table <- attr(cfg, "table")
out_dir <- pipeline_output_dir(cfg)

file_rows <- cfg_table[grepl("(_file$|_chain$)", key) & !is.na(value)]
manifest <- rbindlist(
  Map(file_manifest_row, file_rows$key, file_rows$value),
  use.names = TRUE,
  fill = TRUE
)

write_tsv(manifest, file.path(out_dir, "input_file_manifest.tsv"))
write_tsv(cfg_table, file.path(out_dir, "config_snapshot.tsv"))

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "R_session_info.txt")
)

git_status <- tryCatch(
  system2("git", c("status", "--short"), stdout = TRUE, stderr = TRUE),
  error = function(e) "git status unavailable"
)
writeLines(git_status, file.path(out_dir, "git_status.txt"))

message("Input file manifest written to: ", out_dir)
