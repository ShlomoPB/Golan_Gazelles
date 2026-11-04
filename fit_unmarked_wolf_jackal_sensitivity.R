#!/usr/bin/env Rscript

# Refit only the matched unmarked sensitivity models used to assess whether
# the primary predator gradients depend on the N-mixture implementation.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(unmarked)
})

args <- commandArgs(trailingOnly = TRUE)
arg_map <- list(
  site_cov_file = "SiteCovs_060324.csv",
  obs_cov_file = "ObsCovs_050524.csv",
  output_prefix = "publication_unmarked_wolf_jackal_sensitivity",
  seed = "20260623"
)
for (arg in args) {
  parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    key <- sub("^--", "", parts[1])
    if (!key %in% names(arg_map)) stop("Unknown argument: ", key)
    arg_map[[key]] <- parts[2]
  }
}

site_cov_file <- arg_map$site_cov_file
obs_cov_file <- arg_map$obs_cov_file
output_prefix <- arg_map$output_prefix
seed <- as.integer(arg_map$seed)

species_info <- list(
  list(
    label = "Fawn",
    code = "fawn",
    count_file = "gazelle_juvenile_abundance_60x7.csv"
  ),
  list(
    label = "Female",
    code = "female",
    count_file = "gazelle_female_abundance_60x7.csv"
  ),
  list(
    label = "Male",
    code = "male",
    count_file = "gazelle_male_abundance_60x7.csv"
  )
)

abundance_formula <- ~ wolf_z * jackal_z
detection_formula <-
  ~ Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera

z_scale <- function(x) {
  scale_sd <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(scale_sd) || scale_sd == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / scale_sd
}

read_count_matrix <- function(path) {
  counts <- read.csv(path, na.strings = c("NA", ""), check.names = FALSE)
  y <- as.matrix(counts[paste0("X", 1:7)])
  storage.mode(y) <- "numeric"
  rownames(y) <- counts$camera_name
  colnames(y) <- paste0("X", 1:7)
  y
}

to_site_by_round <- function(data, column, site_ids) {
  matrix(
    data[[column]],
    nrow = length(site_ids),
    ncol = 7,
    byrow = TRUE,
    dimnames = list(site = site_ids, replicate = paste0("X", 1:7))
  )
}

make_observation_covariates <- function(obs_covs, site_ids) {
  list(
    Path.Type = to_site_by_round(obs_covs, "Path.Type", site_ids),
    Camera.Model = to_site_by_round(obs_covs, "Camera.Model", site_ids),
    Vegetation.Structure =
      to_site_by_round(obs_covs, "Vegetation.Structure", site_ids),
    Area.Behind.Camera =
      to_site_by_round(obs_covs, "Area.Behind.Camera", site_ids)
  )
}

tidy_component <- function(fit, species, component, type) {
  estimates <- coef(fit, type = type)
  std_errors <- unmarked::SE(fit, type = type)
  z_values <- estimates / std_errors
  p_values <- 2 * stats::pnorm(abs(z_values), lower.tail = FALSE)
  term_names <- sub("^[^(]+\\(", "", names(estimates))
  term_names <- sub("\\)$", "", term_names)
  term_names[term_names == "Int"] <- "(Intercept)"
  tibble(
    species,
    component,
    term = term_names,
    estimate = unname(estimates),
    std_error = unname(std_errors),
    z_value = unname(z_values),
    p_value = unname(p_values),
    conf_low = estimate - stats::qnorm(0.975) * std_error,
    conf_high = estimate + stats::qnorm(0.975) * std_error
  )
}

site_covs_raw <- read.csv(
  site_cov_file,
  na.strings = c("NA", ""),
  check.names = FALSE
)
site_covs_raw$Camera.ID <- as.integer(site_covs_raw$Camera.ID)

obs_covs_raw <- read.csv(
  obs_cov_file,
  na.strings = c("NA", ""),
  check.names = FALSE
)
obs_covs_raw <- obs_covs_raw[
  order(obs_covs_raw$`Camera Number`, obs_covs_raw$`Round Number`),
]

fits <- list()
coefficient_rows <- list()
diagnostic_rows <- list()

for (species_i in species_info) {
  y <- read_count_matrix(species_i$count_file)
  site_ids <- rownames(y)

  expected <- expand.grid(
    `Camera Number` = as.integer(site_ids),
    `Round Number` = 1:7,
    check.names = FALSE
  )
  expected <- expected[order(expected$`Camera Number`, expected$`Round Number`), ]
  stopifnot(identical(obs_covs_raw$`Camera Number`, expected$`Camera Number`))
  stopifnot(identical(obs_covs_raw$`Round Number`, expected$`Round Number`))

  site_covs <- site_covs_raw[
    match(as.integer(site_ids), site_covs_raw$Camera.ID),
  ] %>%
    mutate(
      wolf_z = z_scale(Null.Wolf.Abundance),
      jackal_z = z_scale(null.jackal.abundance)
    )

  frame <- unmarkedFramePCount(
    y = y,
    siteCovs = data.frame(
      wolf_z = site_covs$wolf_z,
      jackal_z = site_covs$jackal_z,
      row.names = site_ids
    ),
    obsCovs = make_observation_covariates(obs_covs_raw, site_ids)
  )

  model_formula <-
    ~ Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~
      wolf_z * jackal_z
  model_k <- max(150, max(y, na.rm = TRUE) + 50)
  set.seed(seed)
  fit <- pcount(
    model_formula,
    data = frame,
    mixture = "NB",
    K = model_k
  )

  fits[[species_i$code]] <- fit
  coefficient_rows[[length(coefficient_rows) + 1]] <- bind_rows(
    tidy_component(fit, species_i$label, "abundance", "state"),
    tidy_component(fit, species_i$label, "detection", "det")
  )

  optimizer <- fit@opt
  diagnostic_rows[[length(diagnostic_rows) + 1]] <- tibble(
    species = species_i$label,
    n_sites = nrow(y),
    n_replicates = ncol(y),
    observed_count_total = sum(y, na.rm = TRUE),
    missing_count_cells = sum(is.na(y)),
    K = model_k,
    AIC = fit@AIC,
    convergence_code = if (!is.null(optimizer$convergence)) {
      optimizer$convergence
    } else {
      NA_integer_
    },
    optimizer_message = if (!is.null(optimizer$message)) {
      as.character(optimizer$message)
    } else {
      NA_character_
    }
  )
}

coefficient_table <- bind_rows(coefficient_rows)
diagnostic_table <- bind_rows(diagnostic_rows)
ecological_table <- coefficient_table %>%
  filter(
    component == "abundance",
    term %in% c("wolf_z", "jackal_z", "wolf_z:jackal_z")
  )

saveRDS(fits, paste0(output_prefix, "_fits.rds"))
write_csv(coefficient_table, paste0(output_prefix, "_all_coefficients.csv"))
write_csv(ecological_table, paste0(output_prefix, "_ecological_coefficients.csv"))
write_csv(diagnostic_table, paste0(output_prefix, "_diagnostics.csv"))

metadata <- c(
  paste0("Generated: ", format(Sys.time(), tz = "UTC"), " UTC"),
  paste0("Seed: ", seed),
  paste0("Abundance formula: ", paste(deparse(abundance_formula), collapse = " ")),
  paste0("Detection formula: ", paste(deparse(detection_formula), collapse = " ")),
  "Family: negative binomial",
  "K rule: max(150, maximum observed count + 50)",
  "",
  capture.output(sessionInfo())
)
writeLines(metadata, paste0(output_prefix, "_metadata.txt"))

cat("Saved matched unmarked sensitivity outputs with prefix:", output_prefix, "\n")
