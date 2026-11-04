#!/usr/bin/env Rscript

analysis_script <- "analyze_gazelle_pooled_predator_overlap_open_spaces_quartiles.R"
if (!file.exists(analysis_script)) {
  stop("Missing shared analysis script: ", analysis_script)
}

analysis_args <- c(
  analysis_script,
  "--site_covariate=human.abundance",
  "--stratum_label=human activity",
  "--stratum_slug=human-activity",
  "--output_prefix=publication_gazelle_pooled_predator_overlap_by_human_activity_quartile",
  commandArgs(trailingOnly = TRUE)
)

status <- system2(
  file.path(R.home("bin"), "Rscript.exe"),
  args = shQuote(analysis_args)
)
if (!identical(status, 0L)) {
  stop("Human-activity overlap analysis failed with status ", status)
}
