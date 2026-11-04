#!/usr/bin/env Rscript

source_files <- c(
  "manuscript_spabundance_methods_results.md",
  "manuscript_fawn_female_beta_binomial_methods_results.md",
  "manuscript_diel_activity_methods_results.md"
)
output_file <- "conbio_methods_results_raw_package.md"

missing_files <- source_files[!file.exists(source_files)]
if (length(missing_files)) {
  stop("Missing manuscript source file(s): ", paste(missing_files, collapse = ", "))
}

header <- c(
  "# Final Methods and Results Analysis Package",
  "",
  paste0(
    "Generated from the retained corrected-data analyses on ", Sys.Date(),
    ". Early/Late and daylight diel models, seasonal N-mixture models, ",
    "gazelle-predator temporal-overlap models, PAMMs, and generated maps ",
    "are excluded."
  ),
  ""
)

sections <- unlist(lapply(source_files, function(path) {
  c(readLines(path, warn = FALSE, encoding = "UTF-8"), "")
}), use.names = FALSE)

writeLines(c(header, sections), output_file, useBytes = TRUE)
cat("Assembled final manuscript package:", output_file, "\n")
