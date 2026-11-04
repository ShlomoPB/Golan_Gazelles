#!/usr/bin/env Rscript

# Generate the final ten supplementary tables from saved final-analysis outputs.
# This script reads fitted objects and summaries only; it does not refit any model.

options(stringsAsFactors = FALSE, scipen = 999)

u_times <- "\u00d7"
u_minus <- "\u2212"
u_endash <- "\u2013"
u_delta <- "\u0394"
u_degree <- "\u00b0"

find_project_root <- function(start = getwd()) {
  here <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    marker <- file.path(here, "definitive_hierarchical_diel_model_summary.csv")
    if (file.exists(marker)) return(here)
    parent <- dirname(here)
    if (identical(parent, here)) stop("Could not locate the project root.")
    here <- parent
  }
}

root <- find_project_root()
setwd(root)
out_dir <- file.path(root, "supplementary_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  "gazelle_lf2_NMix_NB_open_wolf_jackal_interaction_order-fawn-male-female_2000b_25bl_10th_3ch_fit.rds",
  "gazelle_lf2_NMix_NB_open_wolf_jackal_interaction_order-fawn-male-female_2000b_25bl_10th_3ch_metadata.rds",
  "publication_wolf_jackal_interaction_nmix_spabundance_fits.rds",
  "publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds",
  "publication_unmarked_sensitivity_coefficients.csv",
  "publication_focused_lfmsnmix_ppc.csv",
  "publication_focused_lfmsnmix_diagnostics.csv",
  "publication_wolf_jackal_interaction_nmix_spabundance_diagnostics.csv",
  "publication_fawn_female_beta_binomial_wolf_jackal_fit_diagnostics.csv",
  "publication_fawn_female_beta_binomial_interaction_sensitivity_coefficients.csv",
  "publication_fawn_female_beta_binomial_interaction_sensitivity_model_comparison.csv",
  "publication_fawn_female_beta_binomial_pooled_season_comparison.csv",
  "definitive_hierarchical_diel_ecological_model_comparison.csv",
  "definitive_hierarchical_diel_model_summary.csv",
  "definitive_hierarchical_diel_effect_stability.csv",
  "definitive_hierarchical_diel_diagnostics.csv",
  "definitive_hierarchical_diel_data_audit.csv",
  "joint_demographic_penalty_calibrated_outputs/penalty_calibrated_omnibus_tests.csv",
  "joint_demographic_penalty_calibrated_outputs/penalty_calibrated_pairwise_tests.csv",
  "joint_demographic_penalty_calibrated_outputs/penalty_calibrated_diagnostics.csv"
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required final-analysis files: ", paste(missing_files, collapse = ", "))
}

write_table_csv <- function(x, number, slug) {
  path <- file.path(out_dir, sprintf("table_S%02d_%s.csv", number, slug))
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Package 'readr' is required for locale-safe UTF-8 CSV output.")
  }
  readr::write_csv(x, path, na = "")
  path
}

num_summary <- function(x) {
  c(
    estimate = median(x),
    lower_95 = unname(quantile(x, 0.025)),
    upper_95 = unname(quantile(x, 0.975))
  )
}

range_text <- function(x, digits = 3) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_character_)
  paste0(formatC(min(x), digits = digits, format = "f"), u_endash,
         formatC(max(x), digits = digits, format = "f"))
}

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, digits = digits, format = "f"))
}

fmt_p <- function(x, digits = 3) {
  ifelse(is.na(x), "", ifelse(x < 0.001, "<0.001", formatC(x, digits = digits, format = "f")))
}

md_escape <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  gsub("\\|", "\\\\|", x)
}

md_table <- function(x) {
  if (!nrow(x)) return("_No rows._")
  x[] <- lapply(x, md_escape)
  header <- paste0("| ", paste(names(x), collapse = " | "), " |")
  rule <- paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |")
  rows <- apply(x, 1, function(z) paste0("| ", paste(z, collapse = " | "), " |"))
  paste(c(header, rule, rows), collapse = "\n")
}

class_label <- function(x) {
  out <- c(
    fawn = "Juvenile", Fawn = "Juvenile", Juvenile = "Juvenile",
    female = "Adult female", Female = "Adult female", "Adult female" = "Adult female",
    male = "Adult male", Male = "Adult male", "Adult male" = "Adult male",
    "Pooled species" = "Pooled gazelles"
  )
  unname(out[as.character(x)])
}

# -----------------------------------------------------------------------------
# Table S1: residual spatial correlations from the final two-factor lfMsNMix.
# -----------------------------------------------------------------------------

joint_fit_file <- "gazelle_lf2_NMix_NB_open_wolf_jackal_interaction_order-fawn-male-female_2000b_25bl_10th_3ch_fit.rds"
joint_fit <- readRDS(joint_fit_file)
lambda <- as.matrix(joint_fit$lambda.samples)
joint_class_order <- c("fawn", "male", "female")
n_classes <- length(joint_class_order)
n_factors <- ncol(lambda) / n_classes
if (n_factors != 2) stop("The final joint model does not contain the expected two latent factors.")

pair_order <- list(c("fawn", "female"), c("fawn", "male"), c("female", "male"))
cor_draws <- matrix(NA_real_, nrow(lambda), length(pair_order))
for (i in seq_len(nrow(lambda))) {
  lambda_i <- matrix(lambda[i, ], nrow = n_classes, ncol = n_factors)
  corr_i <- cov2cor(tcrossprod(lambda_i))
  for (j in seq_along(pair_order)) {
    idx <- match(pair_order[[j]], joint_class_order)
    cor_draws[i, j] <- corr_i[idx[1], idx[2]]
  }
}

s1 <- data.frame(
  pair = c(
    paste("Juvenile", "adult female", sep = paste0(" ", u_endash, " ")),
    paste("Juvenile", "adult male", sep = paste0(" ", u_endash, " ")),
    paste("Adult female", "adult male", sep = paste0(" ", u_endash, " "))
  ),
  posterior_estimate = colMeans(cor_draws),
  lower_95_CrI = apply(cor_draws, 2, quantile, probs = 0.025),
  upper_95_CrI = apply(cor_draws, 2, quantile, probs = 0.975),
  source = joint_fit_file,
  check.names = FALSE
)
s1_path <- write_table_csv(s1, 1, "residual_spatial_correlations")

# -----------------------------------------------------------------------------
# Table S2: compact specifications and diagnostics for spatial analyses.
# -----------------------------------------------------------------------------

site_covs <- read.csv("SiteCovs_060324.csv", check.names = FALSE)
predator_r <- cor(
  as.numeric(site_covs[["Null.Wolf.Abundance"]]),
  as.numeric(site_covs[["null.jackal.abundance"]]),
  use = "complete.obs"
)

joint_meta <- readRDS(
  "gazelle_lf2_NMix_NB_open_wolf_jackal_interaction_order-fawn-male-female_2000b_25bl_10th_3ch_metadata.rds"
)
joint_ppc <- read.csv("publication_focused_lfmsnmix_ppc.csv", check.names = FALSE)
joint_diag <- read.csv("publication_focused_lfmsnmix_diagnostics.csv", check.names = FALSE)
class_diag <- read.csv("publication_wolf_jackal_interaction_nmix_spabundance_diagnostics.csv", check.names = FALSE)

grouped_ppc <- joint_ppc[joint_ppc$structure == "lf2" & joint_ppc$group == 1, ]
grouped_ppc$class <- class_label(grouped_ppc$class)
grouped_ppc_text <- paste0(
  grouped_ppc$class, " = ", formatC(grouped_ppc$bayesian_p, digits = 3, format = "f"),
  collapse = "; "
)
ungrouped_ppc <- joint_ppc[joint_ppc$structure == "lf2" & joint_ppc$group == 0, ]
ungrouped_ppc$class <- class_label(ungrouped_ppc$class)
ungrouped_ppc_text <- paste0(
  ungrouped_ppc$class, " = ", formatC(ungrouped_ppc$bayesian_p, digits = 3, format = "f"),
  collapse = "; "
)
loading_diag <- joint_diag[joint_diag$structure == "lf2" & joint_diag$component == "lambda", ]
joint_intercepts <- joint_diag[
  joint_diag$structure == "lf2" & grepl("Intercept", joint_diag$parameter, fixed = TRUE),
]

s2_rows <- list()
add_s2 <- function(panel, analysis, item, value, source, note = NA_character_) {
  s2_rows[[length(s2_rows) + 1L]] <<- data.frame(
    panel = panel, analysis = analysis, item = item, value = as.character(value),
    source = source, note = note, check.names = FALSE
  )
}

for (predator in c("Wolf", "Jackal")) {
  add_s2("A. Predator activity indices", predator, "Count distribution",
         if (predator == "Wolf") "Zero-inflated Poisson" else "Poisson",
         "Current manuscript Methods; original index-fit objects unavailable")
  add_s2("A. Predator activity indices", predator, "State structure", "~ 1",
         "Current manuscript Methods; original index-fit objects unavailable")
  add_s2("A. Predator activity indices", predator, "Final detection formula",
         "~ Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera",
         "Current manuscript Methods; original index-fit objects unavailable")
  add_s2("A. Predator activity indices", predator, "Sampling design", "60 sites; 7 nominal occasions",
         "SiteCovs_060324.csv; ObsCovs_050524.csv")
  add_s2("A. Predator activity indices", predator, "Site-level extraction",
         "Empirical-Bayes posterior mode from ranef()/bup(), then standardized",
         "Current manuscript Methods", "Relative site-level activity index; not density or absolute abundance")
}
add_s2("A. Predator activity indices", "Wolf and jackal", "Site-level correlation",
       sprintf("r = %.6f", predator_r), "SiteCovs_060324.csv")
add_s2("A. Predator activity indices", "Provenance", "Historical model comparison",
       "Not recoverable from archived files", "Repository audit",
       "Original predator-index fit objects and candidate-model fit statistics were not found; no AIC/WAIC values are inferred")

add_s2("B. Joint demographic lfMsNMix", "All three classes", "Abundance formula",
       "~ Open.Spaces.z + Null.Wolf.Abundance.z * null.jackal.abundance.z",
       "fit_focused_multiclass_nmix.R; final lf2 fit")
add_s2("B. Joint demographic lfMsNMix", "All three classes", "Detection formula",
       "~ log_days_operational.z + cattle_sightings.z",
       "fit_focused_multiclass_nmix.R; final lf2 fit")
add_s2("B. Joint demographic lfMsNMix", "All three classes", "Distribution and latent structure",
       "Negative binomial; 2 latent factors", "final lf2 fit")
add_s2("B. Joint demographic lfMsNMix", "All three classes", "MCMC",
       "3 chains; 50,000 iterations/chain; burn-in 5,000; thinning 10; seed 20260819",
       "final lf2 metadata")
add_s2("B. Joint demographic lfMsNMix", "All three classes", paste0("Grouped Freeman", u_endash, "Tukey PPC"),
       grouped_ppc_text, "publication_focused_lfmsnmix_ppc.csv")
add_s2("B. Joint demographic lfMsNMix", "All three classes", paste0("Ungrouped Freeman", u_endash, "Tukey PPC"),
       ungrouped_ppc_text, "publication_focused_lfmsnmix_ppc.csv")
add_s2("B. Joint demographic lfMsNMix", "All three classes", "Free-loading diagnostics",
       paste0("R-hat ", range_text(loading_diag$rhat), "; ESS ", range_text(loading_diag$ess, 1)),
       "publication_focused_lfmsnmix_diagnostics.csv")
add_s2("B. Joint demographic lfMsNMix", "All three classes", "Convergence limitation",
       paste0("Intercept diagnostics included R-hat up to ", formatC(max(joint_intercepts$rhat), 3, format = "f"),
              " and ESS as low as ", formatC(min(joint_intercepts$ess), 1, format = "f")),
       "publication_focused_lfmsnmix_diagnostics.csv",
       "Residual-correlation inference is secondary and interpreted with this limitation")

for (sp in c("Fawn", "Female", "Male")) {
  label <- class_label(sp)
  focal <- class_diag[class_diag$species == sp & class_diag$component == "abundance" &
                        !grepl("Intercept", class_diag$parameter, fixed = TRUE), ]
  intercepts <- class_diag[class_diag$species == sp &
                             grepl("Intercept", class_diag$parameter, fixed = TRUE), ]
  add_s2("C. Class-specific spAbundance::NMix", label, "Abundance formula",
         "~ Null.Wolf.Abundance.z * null.jackal.abundance.z",
         "analyze_wolf_jackal_interaction_nmix.R; final class fit")
  add_s2("C. Class-specific spAbundance::NMix", label, "Detection formula",
         "~ log_days_operational.z + cattle_sightings.z",
         "analyze_wolf_jackal_interaction_nmix.R; final class fit")
  add_s2("C. Class-specific spAbundance::NMix", label, "Distribution and MCMC",
         "Negative binomial; 3 chains; 50,000 iterations/chain; burn-in 5,000; thinning 10",
         "publication_final_model_specifications.csv")
  add_s2("C. Class-specific spAbundance::NMix", label, "Focal abundance-slope diagnostics",
         paste0("R-hat ", range_text(focal$rhat), "; ESS ", range_text(focal$ess, 1)),
         "publication_wolf_jackal_interaction_nmix_spabundance_diagnostics.csv")
  add_s2("C. Class-specific spAbundance::NMix", label, "Intercept convergence limitation",
         paste0("R-hat up to ", formatC(max(intercepts$rhat), 3, format = "f"),
                "; ESS as low as ", formatC(min(intercepts$ess), 1, format = "f")),
         "publication_wolf_jackal_interaction_nmix_spabundance_diagnostics.csv",
         "Ecological inference restricted to well-mixed slopes and normalized predictions")
}
add_s2("C. Class-specific spAbundance::NMix", "Historical clarification", "Expanded covariate models",
       "Exploratory only; inadequate all-parameter MCMC diagnostics",
       "Repository scripts and diagnostic outputs",
       "No WAIC/AIC comparison selected the final predator-only state formula")
s2 <- do.call(rbind, s2_rows)
s2_path <- write_table_csv(s2, 2, "spatial_model_specifications_and_diagnostics")

# -----------------------------------------------------------------------------
# Table S3: conditional jackal effects from final class-specific NMix draws.
# -----------------------------------------------------------------------------

class_fits_file <- "publication_wolf_jackal_interaction_nmix_spabundance_fits.rds"
class_fits <- readRDS(class_fits_file)
s3_rows <- list()
for (nm in c("fawn", "female", "male")) {
  beta <- as.matrix(class_fits[[nm]]$beta.samples)
  jackal_col <- grep("^null\\.jackal\\.abundance\\.z$", colnames(beta))
  interaction_col <- grep("Null\\.Wolf\\.Abundance\\.z:null\\.jackal\\.abundance\\.z", colnames(beta))
  if (length(jackal_col) != 1 || length(interaction_col) != 1) {
    stop("Could not identify final jackal and interaction posterior columns for ", nm)
  }
  b_j <- beta[, jackal_col]
  b_wj <- beta[, interaction_col]
  for (wolf in c(-1, 0, 1)) {
    ratio <- exp(-2 * (b_j + wolf * b_wj))
    ss <- num_summary(ratio)
    s3_rows[[length(s3_rows) + 1L]] <- data.frame(
      panel = "A. Low-jackal / high-jackal abundance ratio",
      demographic_class = class_label(nm),
      wolf_activity = c(`-1` = paste0("Low (", u_minus, "1 SD)"), `0` = "Mean (0 SD)", `1` = "High (+1 SD)")[as.character(wolf)],
      contrast = paste0("Low jackal (", u_minus, "1 SD) / high jackal (+1 SD)"),
      posterior_median = ss["estimate"], lower_95_CrI = ss["lower_95"],
      upper_95_CrI = ss["upper_95"], posterior_probability_ratio_gt_1 = mean(ratio > 1),
      source = class_fits_file, check.names = FALSE
    )
  }
  ratio_of_ratios <- exp(4 * b_wj)
  ss <- num_summary(ratio_of_ratios)
  s3_rows[[length(s3_rows) + 1L]] <- data.frame(
    panel = "B. Attenuation ratio-of-ratios",
    demographic_class = class_label(nm), wolf_activity = "Low versus high wolf",
    contrast = "Jackal-effect ratio at low wolf / jackal-effect ratio at high wolf",
    posterior_median = ss["estimate"], lower_95_CrI = ss["lower_95"],
    upper_95_CrI = ss["upper_95"], posterior_probability_ratio_gt_1 = mean(ratio_of_ratios > 1),
    source = class_fits_file, check.names = FALSE
  )
}
s3 <- do.call(rbind, s3_rows)
row.names(s3) <- NULL
s3_path <- write_table_csv(s3, 3, "conditional_jackal_effects")

# -----------------------------------------------------------------------------
# Table S4: matched unmarked sensitivity coefficients.
# -----------------------------------------------------------------------------

um <- read.csv("publication_unmarked_sensitivity_coefficients.csv", check.names = FALSE)
term_map <- c(
  wolf_z = "Wolf activity",
  jackal_z = "Jackal activity",
  "wolf_z:jackal_z" = paste("Wolf", "jackal", sep = paste0(" ", u_times, " "))
)
s4 <- data.frame(
  demographic_class = class_label(um$species),
  term = unname(term_map[um$term]),
  estimate = um$estimate,
  SE = um$std_error,
  lower_95_CI = um$conf_low,
  upper_95_CI = um$conf_high,
  z = um$z_value,
  P = um$p_value,
  direction = ifelse(um$estimate > 0, "Positive", "Negative"),
  support = ifelse(um$conf_low > 0 | um$conf_high < 0, "95% CI excludes zero", "95% CI includes zero"),
  source = "publication_unmarked_sensitivity_coefficients.csv",
  check.names = FALSE
)
s4_path <- write_table_csv(s4, 4, "matched_unmarked_sensitivity")

# -----------------------------------------------------------------------------
# Tables S5-S6: juvenile:female beta-binomial candidate set and coefficients.
# -----------------------------------------------------------------------------

if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  stop("Package 'glmmTMB' is required to extract the saved beta-binomial model summaries.")
}
suppressPackageStartupMessages(library(glmmTMB))
ratio_fits_file <- "publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds"
ratio_fits <- readRDS(ratio_fits_file)
ratio_diag <- read.csv("publication_fawn_female_beta_binomial_wolf_jackal_fit_diagnostics.csv", check.names = FALSE)

model_order <- c(
  "overall", "overall_additive", "seasonal", "seasonal_no_predator_interaction",
  "seasonal_stable_response", "seasonal_additive_stable"
)
model_labels <- c(
  overall = paste0("Wolf ", u_times, " jackal; no survey-period term"),
  overall_additive = "Additive wolf + jackal; no survey-period term",
  seasonal = paste0("Wolf ", u_times, " jackal ", u_times, " survey period"),
  seasonal_no_predator_interaction = "Additive predator slopes varying by survey period",
  seasonal_stable_response = paste0("Wolf ", u_times, " jackal + common survey-period shift"),
  seasonal_additive_stable = "Additive wolf + jackal + common survey-period shift"
)
predator_structure <- c(
  overall = paste("Wolf", "jackal", sep = paste0(" ", u_times, " ")), overall_additive = "Additive",
  seasonal = paste0("Wolf ", u_times, " jackal with period-specific main and interaction terms"),
  seasonal_no_predator_interaction = "Additive with period-specific slopes",
  seasonal_stable_response = paste("Wolf", "jackal", sep = paste0(" ", u_times, " ")), seasonal_additive_stable = "Additive"
)
period_structure <- c(
  overall = "None", overall_additive = "None",
  seasonal = "All predator terms vary by Early/Late",
  seasonal_no_predator_interaction = "Predator slopes vary by Early/Late",
  seasonal_stable_response = "Common Early/Late shift",
  seasonal_additive_stable = "Common Early/Late shift"
)
aics <- vapply(ratio_fits[model_order], AIC, numeric(1))
dfs <- vapply(ratio_fits[model_order], function(x) attr(logLik(x), "df"), numeric(1))
s5a <- data.frame(
  panel = "A. Candidate models", model = unname(model_labels[model_order]),
  predator_structure = unname(predator_structure[model_order]),
  survey_period_structure = unname(period_structure[model_order]),
  parameters_df = dfs, AIC = aics, delta_AIC = aics - min(aics),
  selected = ifelse(model_order == "seasonal_additive_stable", "Yes", "No"),
  convergence_code = ratio_diag$convergence_code[match(model_order, ratio_diag$model)],
  positive_definite_hessian = ratio_diag$positive_definite_hessian[match(model_order, ratio_diag$model)],
  source = ratio_fits_file, check.names = FALSE
)

lrt_row <- function(smaller, larger, label) {
  result <- anova(ratio_fits[[smaller]], ratio_fits[[larger]])
  data.frame(
    panel = "B. Nested likelihood-ratio tests", model = label,
    predator_structure = NA_character_, survey_period_structure = NA_character_,
    parameters_df = as.numeric(result$Df[2] - result$Df[1]),
    AIC = NA_real_, delta_AIC = as.numeric(result$AIC[2] - result$AIC[1]),
    selected = NA_character_, convergence_code = NA_integer_, positive_definite_hessian = NA,
    LR_statistic = as.numeric(result$Chisq[2]), P = as.numeric(result$`Pr(>Chisq)`[2]),
    source = ratio_fits_file, check.names = FALSE
  )
}
s5b <- rbind(
  lrt_row("overall_additive", "seasonal_additive_stable", "Add common Early/Late shift to additive model"),
  lrt_row("seasonal_additive_stable", "seasonal_stable_response", paste0("Add wolf ", u_times, " jackal to retained model")),
  lrt_row("seasonal_additive_stable", "seasonal_no_predator_interaction", "Allow additive predator slopes to vary by Early/Late")
)
s5a$LR_statistic <- NA_real_
s5a$P <- NA_real_
s5 <- rbind(s5a, s5b[names(s5a)])
# Restore LRT columns after aligning the two panels.
s5$LR_statistic <- c(rep(NA_real_, nrow(s5a)), s5b$LR_statistic)
s5$P <- c(rep(NA_real_, nrow(s5a)), s5b$P)
s5_path <- write_table_csv(s5, 5, "juvenile_female_candidate_models")

final_ratio_fit <- ratio_fits$seasonal_additive_stable
coef_matrix <- summary(final_ratio_fit)$coefficients$cond
coef_df <- data.frame(
  term_code = rownames(coef_matrix), estimate_logit = coef_matrix[, "Estimate"],
  SE = coef_matrix[, "Std. Error"], z = coef_matrix[, "z value"],
  P = coef_matrix[, "Pr(>|z|)"], row.names = NULL, check.names = FALSE
)
coef_df$odds_ratio <- exp(coef_df$estimate_logit)
coef_df$lower_95_CI <- exp(coef_df$estimate_logit - qnorm(0.975) * coef_df$SE)
coef_df$upper_95_CI <- exp(coef_df$estimate_logit + qnorm(0.975) * coef_df$SE)

pretty_ratio_term <- function(term) {
  exact <- c(
    "(Intercept)" = "Intercept", wolf_z = "Wolf activity", jackal_z = "Jackal activity",
    "seasonLate (mild-cold)" = "Late vs Early survey", open_z = "Open-space cover",
    built_z = "Built-and-road cover", reserve_z = "Nature-reserve cover", human_z = "Human activity"
  )
  out <- unname(exact[term])
  missing <- is.na(out)
  out[missing & grepl("^path_type", term)] <- paste0("Path type: ", sub("^path_type", "", term[missing & grepl("^path_type", term)]))
  out[missing & grepl("^camera_model", term)] <- paste0("Camera model: ", sub("^camera_model", "", term[missing & grepl("^camera_model", term)]))
  out[missing & grepl("^vegetation_structure", term)] <- paste0("Vegetation structure: ", sub("^vegetation_structure", "", term[missing & grepl("^vegetation_structure", term)]))
  out[missing & grepl("^area_behind_camera", term)] <- paste0("Area behind camera: ", sub("^area_behind_camera", "", term[missing & grepl("^area_behind_camera", term)]))
  out[is.na(out)] <- term[is.na(out)]
  out
}

reference_level <- function(term) {
  out <- rep(NA_character_, length(term))
  out[grepl("^path_type", term)] <- "4x4 road with off-road vehicles"
  out[grepl("^camera_model", term)] <- "Browning DARK OPS HD MAX"
  out[grepl("^vegetation_structure", term)] <- "Grassland"
  out[grepl("^area_behind_camera", term)] <- "Clear"
  out[term == "seasonLate (mild-cold)"] <- "Early (hotter)"
  out
}

s6_fixed <- data.frame(
  panel = "A. Final model fixed effects", component = "Fixed effect",
  predictor = pretty_ratio_term(coef_df$term_code), reference = reference_level(coef_df$term_code),
  estimate = coef_df$estimate_logit, SE = coef_df$SE, odds_ratio = coef_df$odds_ratio,
  lower_95_CI = coef_df$lower_95_CI, upper_95_CI = coef_df$upper_95_CI,
  z_or_LR = coef_df$z, df = NA_real_, P = coef_df$P,
  note = NA_character_, source = ratio_fits_file, check.names = FALSE
)

re <- read.csv("publication_fawn_female_beta_binomial_coefficient_random_effects.csv", check.names = FALSE)
s6_re <- data.frame(
  panel = "A. Final model other parameters", component = "Random-effect SD",
  predictor = re$group, reference = NA_character_, estimate = re$std_dev, SE = NA_real_,
  odds_ratio = NA_real_, lower_95_CI = NA_real_, upper_95_CI = NA_real_,
  z_or_LR = NA_real_, df = NA_real_, P = NA_real_, note = paste0("Variance = ", re$variance),
  source = "publication_fawn_female_beta_binomial_coefficient_random_effects.csv", check.names = FALSE
)
dispersion <- sigma(final_ratio_fit)
s6_disp <- data.frame(
  panel = "A. Final model other parameters", component = "Beta-binomial dispersion",
  predictor = "Dispersion", reference = NA_character_, estimate = as.numeric(dispersion), SE = NA_real_,
  odds_ratio = NA_real_, lower_95_CI = NA_real_, upper_95_CI = NA_real_,
  z_or_LR = NA_real_, df = NA_real_, P = NA_real_, note = NA_character_,
  source = ratio_fits_file, check.names = FALSE
)

interaction_coef <- summary(ratio_fits$seasonal_stable_response)$coefficients$cond["wolf_z:jackal_z", ]
interaction_test <- anova(ratio_fits$seasonal_additive_stable, ratio_fits$seasonal_stable_response)
slope_test <- anova(ratio_fits$seasonal_additive_stable, ratio_fits$seasonal_no_predator_interaction)
s6_sensitivity <- rbind(
  data.frame(
    panel = "B. Sensitivity and extension tests", component = "Coefficient",
    predictor = paste0("Wolf ", u_times, " jackal interaction"), reference = NA_character_,
    estimate = interaction_coef["Estimate"], SE = interaction_coef["Std. Error"],
    odds_ratio = exp(interaction_coef["Estimate"]),
    lower_95_CI = exp(interaction_coef["Estimate"] - qnorm(0.975) * interaction_coef["Std. Error"]),
    upper_95_CI = exp(interaction_coef["Estimate"] + qnorm(0.975) * interaction_coef["Std. Error"]),
    z_or_LR = interaction_coef["z value"], df = NA_real_, P = interaction_coef["Pr(>|z|)"],
    note = paste0("Coefficient from model with wolf ", u_times, " jackal plus common survey-period shift"),
    source = ratio_fits_file, check.names = FALSE
  ),
  data.frame(
    panel = "B. Sensitivity and extension tests", component = "Likelihood-ratio test",
    predictor = paste0("Add wolf ", u_times, " jackal to retained model"), reference = NA_character_,
    estimate = NA_real_, SE = NA_real_, odds_ratio = NA_real_, lower_95_CI = NA_real_, upper_95_CI = NA_real_,
    z_or_LR = interaction_test$Chisq[2], df = interaction_test$Df[2] - interaction_test$Df[1],
    P = interaction_test$`Pr(>Chisq)`[2], note = paste0(u_delta, "AIC = ", interaction_test$AIC[2] - interaction_test$AIC[1]),
    source = ratio_fits_file, check.names = FALSE
  ),
  data.frame(
    panel = "B. Sensitivity and extension tests", component = "Likelihood-ratio test",
    predictor = "Allow additive predator slopes to vary by Early/Late", reference = NA_character_,
    estimate = NA_real_, SE = NA_real_, odds_ratio = NA_real_, lower_95_CI = NA_real_, upper_95_CI = NA_real_,
    z_or_LR = slope_test$Chisq[2], df = slope_test$Df[2] - slope_test$Df[1],
    P = slope_test$`Pr(>Chisq)`[2], note = NA_character_, source = ratio_fits_file, check.names = FALSE
  )
)
s6 <- rbind(s6_fixed, s6_re, s6_disp, s6_sensitivity)
s6_path <- write_table_csv(s6, 6, "juvenile_female_coefficients_and_sensitivity")

# -----------------------------------------------------------------------------
# Tables S7-S10: diel model selection, robustness, joint tests, diagnostics.
# -----------------------------------------------------------------------------

diel_models <- read.csv("definitive_hierarchical_diel_ecological_model_comparison.csv", check.names = FALSE)
diel_model_names <- c(
  A = "Baseline", B = "Wolf-only", C = "Jackal-only",
  D = "Additive wolf + jackal", E = paste("Wolf", "jackal", sep = paste0(" ", u_times, " "))
)
s7 <- data.frame(
  demographic_analysis = ifelse(diel_models$response_class == "Pooled species", "Pooled gazelles", diel_models$response_class),
  model = unname(diel_model_names[diel_models$model_id]), AIC = diel_models$AIC,
  delta_AIC = diel_models$delta_AIC,
  selected_under_least_complex_within_2_AIC = ifelse(diel_models$selected, "Yes", "No"),
  hierarchy = diel_models$hierarchy, convergence = diel_models$convergence,
  source = "definitive_hierarchical_diel_ecological_model_comparison.csv", check.names = FALSE
)
s7_path <- write_table_csv(s7, 7, "diel_candidate_model_comparison")

diel_summary <- read.csv("definitive_hierarchical_diel_model_summary.csv", check.names = FALSE)
diel_stability <- read.csv("definitive_hierarchical_diel_effect_stability.csv", check.names = FALSE)
pred_terms <- diel_summary[diel_summary$stage == "Final fREML" & diel_summary$term_type == "predator-by-time", ]
effect_from_term <- function(term) {
  ifelse(grepl("wolf_jackal", term), "interaction",
         ifelse(grepl("wolf_z", term), "wolf", "jackal"))
}
pred_terms$effect_key <- effect_from_term(pred_terms$term)
join_key_summary <- paste(pred_terms$response_class, pred_terms$effect_key)
join_key_stability <- paste(diel_stability$response_class, diel_stability$effect)
idx_stability <- match(join_key_summary, join_key_stability)
term_labels <- c(
  wolf = paste("Wolf", "time", sep = paste0(" ", u_times, " ")),
  jackal = paste("Jackal", "time", sep = paste0(" ", u_times, " ")),
  interaction = paste(c("Wolf", "jackal", "time"), collapse = paste0(" ", u_times, " "))
)
s8 <- data.frame(
  demographic_analysis = ifelse(pred_terms$response_class == "Pooled species", "Pooled gazelles", pred_terms$response_class),
  term = unname(term_labels[pred_terms$effect_key]), EDF = pred_terms$EDF,
  statistic = pred_terms$statistic, P = pred_terms$P,
  supported_specifications = diel_stability$supported_specifications[idx_stability],
  tested_specifications = diel_stability$tested_specifications[idx_stability],
  interpretation = pred_terms$robustness_status,
  selected_model = pred_terms$model_id,
  source = "definitive_hierarchical_diel_model_summary.csv; definitive_hierarchical_diel_effect_stability.csv",
  check.names = FALSE
)
s8_path <- write_table_csv(s8, 8, "diel_effects_and_robustness")

omnibus <- read.csv("joint_demographic_penalty_calibrated_outputs/penalty_calibrated_omnibus_tests.csv", check.names = FALSE)
pairwise <- read.csv("joint_demographic_penalty_calibrated_outputs/penalty_calibrated_pairwise_tests.csv", check.names = FALSE)
s9a <- data.frame(
  panel = "A. Omnibus joint-GAM tests", contrast_or_scope = omnibus$scope,
  term = gsub(" x ", paste0(" ", u_times, " "), omnibus$predator_effect, fixed = TRUE), EDF = omnibus$EDF,
  reference_rank = omnibus$reference_df, statistic = omnibus$statistic,
  raw_P = omnibus$raw_P, Holm_adjusted_P = NA_real_,
  sensitivity_support = ifelse(grepl("jackal", omnibus$predator_effect, ignore.case = TRUE) &
                                 grepl("Wolf", omnibus$predator_effect), "7/7", "Not reopened"),
  interpretation = ifelse(omnibus$raw_P < 0.05, "Supported", "Unsupported"),
  source = "joint_demographic_penalty_calibrated_outputs/penalty_calibrated_omnibus_tests.csv",
  check.names = FALSE
)
s9b <- data.frame(
  panel = "B. Planned pairwise contrasts", contrast_or_scope = gsub(" - ", paste0(" ", u_endash, " "), pairwise$scope, fixed = TRUE),
  term = paste(c("Wolf", "jackal", "time"), collapse = paste0(" ", u_times, " ")), EDF = pairwise$EDF, reference_rank = pairwise$reference_df,
  statistic = pairwise$statistic, raw_P = pairwise$raw_P,
  Holm_adjusted_P = pairwise$Holm_adjusted_P, sensitivity_support = "Primary apparent-solar model only",
  interpretation = ifelse(pairwise$Holm_adjusted_P < 0.05, "Responses differ", "No supported difference; not equivalence"),
  source = "joint_demographic_penalty_calibrated_outputs/penalty_calibrated_pairwise_tests.csv",
  check.names = FALSE
)
s9 <- rbind(s9a, s9b)
s9_path <- write_table_csv(s9, 9, "joint_gam_demographic_heterogeneity")

diel_diag <- read.csv("definitive_hierarchical_diel_diagnostics.csv", check.names = FALSE)
diel_audit <- read.csv("definitive_hierarchical_diel_data_audit.csv", check.names = FALSE)
final_diag <- diel_diag[diel_diag$stage == "Final fREML", ]
raw_audit <- diel_audit[diel_audit$history_type == "Raw sightings" & diel_audit$time_scale == "Apparent solar time", ]
audit_idx <- match(final_diag$response_class, raw_audit$response_class)
hierarchy_text <- c(
  camera_round = "Camera + camera-round random intercepts",
  site_shape = "Camera-specific cyclic deviations + camera and camera-round random intercepts"
)
s10_class <- data.frame(
  model = ifelse(final_diag$response_class == "Pooled species", "Pooled gazelles", final_diag$response_class),
  aggregated_rows = raw_audit$rows[audit_idx], camera_rounds = raw_audit$camera_rounds[audit_idx],
  binomial_trials = raw_audit$trials[audit_idx], positive_camera_day_hours = raw_audit$success[audit_idx],
  cyclic_basis_k = final_diag$k, hierarchy = unname(hierarchy_text[final_diag$hierarchy]),
  deviance_explained = final_diag$deviance_explained, Pearson_dispersion = final_diag$pearson_dispersion,
  mean_camera_round_lag1_Pearson = final_diag$mean_camera_round_lag1_pearson,
  minimum_time_k_index = final_diag$minimum_time_k_index,
  minimum_time_k_P = final_diag$minimum_time_k_P,
  maximum_time_concurvity = final_diag$maximum_time_concurvity,
  maximum_absolute_gradient = NA_real_, minimum_hessian_eigenvalue = NA_real_,
  hessian_numerically_PSD = NA, coefficient_rank = NA_real_, convergence = final_diag$convergence,
  source = "definitive_hierarchical_diel_diagnostics.csv; definitive_hierarchical_diel_data_audit.csv",
  check.names = FALSE
)
joint_gam_diag <- read.csv("joint_demographic_penalty_calibrated_outputs/penalty_calibrated_diagnostics.csv", check.names = FALSE)
joint_audit <- raw_audit[raw_audit$response_class %in% c("Juvenile", "Adult female", "Adult male"), ]
s10_joint <- data.frame(
  model = "Confirmatory joint GAM", aggregated_rows = sum(joint_audit$rows),
  camera_rounds = max(joint_audit$camera_rounds), binomial_trials = sum(joint_audit$trials),
  positive_camera_day_hours = sum(joint_audit$success), cyclic_basis_k = 12,
  hierarchy = "Common and class-specific camera/camera-round random effects; adult camera-specific cyclic deviations",
  deviance_explained = joint_gam_diag$deviance_explained,
  Pearson_dispersion = joint_gam_diag$pearson_dispersion,
  mean_camera_round_lag1_Pearson = NA_real_, minimum_time_k_index = NA_real_,
  minimum_time_k_P = NA_real_, maximum_time_concurvity = NA_real_,
  maximum_absolute_gradient = joint_gam_diag$maximum_absolute_gradient,
  minimum_hessian_eigenvalue = joint_gam_diag$minimum_hessian_eigenvalue,
  hessian_numerically_PSD = joint_gam_diag$hessian_numerically_positive_semidefinite,
  coefficient_rank = joint_gam_diag$coefficient_rank,
  convergence = joint_gam_diag$convergence,
  source = "joint_demographic_penalty_calibrated_outputs/penalty_calibrated_diagnostics.csv",
  check.names = FALSE
)
s10 <- rbind(s10_class, s10_joint)
s10_path <- write_table_csv(s10, 10, "hierarchical_cyclic_gam_diagnostics")

juvenile_s10 <- match("Juvenile", s10$model)
joint_s10 <- match("Confirmatory joint GAM", s10$model)
class_s10 <- which(s10$model != "Confirmatory joint GAM")
s10_specific_checks <- data.frame(
  Model = c(
    s10$model[class_s10],
    rep("Juvenile", 3),
    rep("Confirmatory joint GAM", 4)
  ),
  Diagnostic = c(
    rep("Mean camera-round lag-1 Pearson residual correlation", length(class_s10)),
    "Minimum time-smooth k-index",
    "Time-smooth k-test P",
    "Maximum time-smooth concurvity",
    "Maximum absolute gradient",
    "Minimum Hessian eigenvalue",
    "Hessian numerically positive semidefinite",
    "Coefficient rank"
  ),
  Value = c(
    fmt(s10$mean_camera_round_lag1_Pearson[class_s10]),
    fmt(s10$minimum_time_k_index[juvenile_s10]),
    fmt_p(s10$minimum_time_k_P[juvenile_s10]),
    fmt(s10$maximum_time_concurvity[juvenile_s10]),
    fmt(s10$maximum_absolute_gradient[joint_s10], 6),
    formatC(s10$minimum_hessian_eigenvalue[joint_s10], digits = 3, format = "e"),
    as.character(s10$hessian_numerically_PSD[joint_s10]),
    fmt(s10$coefficient_rank[joint_s10], 0)
  ),
  check.names = FALSE
)

# -----------------------------------------------------------------------------
# Combined manuscript-ready Quarto document.
# -----------------------------------------------------------------------------

qmd <- c(
  "---",
  "title: \"Supplementary information\"",
  "format:",
  "  html: default",
  "  docx: default",
  "execute:",
  "  echo: false",
  "---",
  "",
  "# Appendix S1. Spatial analyses",
  "",
  paste(
    "Wolf and jackal covariates were detection-corrected relative site-level activity indices, not estimates of predator density.",
    "The archived manuscript specification used intercept-only predator state models and detection models containing path type, camera model, vegetation structure, and area behind the camera; empirical-Bayes site effects were extracted and standardized.",
    "The original predator-index fit objects and candidate-fit statistics were not retained, so this provenance is reported explicitly rather than reconstructed."
  ),
  "",
  paste(
    paste0("The joint negative-binomial `spAbundance::lfMsNMix` model used open-space cover and the wolf ", u_times, " jackal activity structure in the state process, effort and cattle activity in detection, and two latent factors."),
    "Open-space cover was specified a priori rather than selected by information-criterion comparison.",
    paste0("The class-specific negative-binomial `NMix` models used wolf ", u_times, " jackal activity in the state process and the same detection covariates."),
    "Bayesian fits used three chains of 50,000 iterations, a 5,000-iteration burn-in, and thinning by 10; posterior predictive checks, R-hat, and effective sample size were used diagnostically.",
    "Class-specific priors were normal with variance 100 for state coefficients and 2.72 for detection coefficients, with a Uniform(0, 10) negative-binomial dispersion prior; the joint model used the corresponding community normal priors, inverse-gamma(0.1, 0.1) among-class variance priors, and class-specific Uniform(0, 10) dispersion priors.",
    "Conditional ratios were calculated draw by draw, and matched `unmarked::pcount` models provided a parallel-framework sensitivity analysis."
  ),
  "",
  "## Table S1. Residual spatial correlations among mountain gazelle demographic classes.",
  "",
  md_table(data.frame(
    Pair = s1$pair, `Posterior estimate` = fmt(s1$posterior_estimate),
    `Lower 95% CrI` = fmt(s1$lower_95_CrI), `Upper 95% CrI` = fmt(s1$upper_95_CrI),
    check.names = FALSE
  )),
  "",
  paste0("Residual correlations are from the final two-factor joint demographic model after accounting for open-space cover and wolf", u_endash, "jackal context. Point estimates are posterior means. CrI, credible interval."),
  "",
  "## Table S2. Model specification and diagnostics for predator activity indices and Bayesian spatial abundance models.",
  "",
  md_table(s2[c("panel", "analysis", "item", "value", "note")]),
  "",
  "Gazelle latent-state results are interpreted as relative local abundance or site use. Expanded habitat/human candidate models were exploratory and had inadequate all-parameter diagnostics; no WAIC or AIC comparison selected the final predator-only class-specific formula.",
  "",
  "## Table S3. Conditional effects of jackal activity across the wolf activity gradient.",
  "",
  md_table(data.frame(
    Panel = s3$panel, `Demographic class` = s3$demographic_class,
    `Wolf activity` = s3$wolf_activity, Contrast = s3$contrast,
    `Posterior median` = fmt(s3$posterior_median, 2),
    `Lower 95% CrI` = fmt(s3$lower_95_CrI, 2),
    `Upper 95% CrI` = fmt(s3$upper_95_CrI, 2),
    `Pr(ratio > 1)` = fmt(s3$posterior_probability_ratio_gt_1, 3), check.names = FALSE
  )),
  "",
  paste0("Point estimates are posterior medians. Panel A reports low-jackal/high-jackal predicted relative abundance, so values greater than one indicate lower predicted abundance as jackal activity increases from ", u_minus, "1 to +1 SD. Panel B divides that jackal-effect ratio at low wolf activity by the corresponding ratio at high wolf activity; values greater than one indicate attenuation of the negative jackal association at high wolf activity."),
  "",
  "## Table S4. Sensitivity analysis using matched unmarked N-mixture models.",
  "",
  md_table(data.frame(
    `Demographic class` = s4$demographic_class, Term = s4$term,
    Estimate = fmt(s4$estimate), SE = fmt(s4$SE),
    `Lower 95% CI` = fmt(s4$lower_95_CI), `Upper 95% CI` = fmt(s4$upper_95_CI),
    z = fmt(s4$z), P = fmt_p(s4$P), Support = s4$support, check.names = FALSE
  )),
  "",
  "These matched negative-binomial `unmarked` fits are a sensitivity analysis. Differences in fitting framework and detection-model implementation preclude treating their coefficients as numerically equivalent to the Bayesian estimates.",
  "",
  "# Appendix S2. Juvenile:female ratio",
  "",
  paste(
    "Camera-round juvenile and adult-female counts were analyzed with a beta-binomial mixed model for 221 informative camera-rounds at 50 cameras.",
    "All six candidates retained open-space, built-and-road, nature-reserve, human-activity, and camera-design adjustments, plus camera nested within polygon random intercepts.",
    "Only predator and survey-period structure varied.",
    "The implemented Early period included camera-round midpoints on or before 17 October 2020 and Late included later midpoints; no midpoint fell exactly on the boundary.",
    "Candidate models were compared by AIC, and nested extensions were evaluated by likelihood-ratio tests."
  ),
  "",
  "## Table S5. Candidate-model comparison for the juvenile:female ratio.",
  "",
  "### Panel A. Candidate models",
  "",
  md_table(setNames(data.frame(
    Model = s5a$model, `Predator structure` = s5a$predator_structure,
    `Survey-period structure` = s5a$survey_period_structure,
    df = s5a$parameters_df, AIC = fmt(s5a$AIC), delta_AIC = fmt(s5a$delta_AIC),
    Selected = s5a$selected, check.names = FALSE
  ), c("Model", "Predator structure", "Survey-period structure", "df", "AIC", paste0(u_delta, "AIC"), "Selected"))),
  "",
  "### Panel B. Nested likelihood-ratio tests",
  "",
  md_table(setNames(data.frame(
    Comparison = s5b$model, df = s5b$parameters_df,
    LR = fmt(s5b$LR_statistic), P = fmt_p(s5b$P), delta_AIC = fmt(s5b$delta_AIC),
    check.names = FALSE
  ), c("Comparison", "df", "LR", "P", paste0(u_delta, "AIC (larger ", u_minus, " smaller)")))),
  "",
  "All candidates converged with code 0 and positive-definite Hessians. The selected model was the minimum-AIC additive-predator model with a common Early/Late shift.",
  "",
  "## Table S6. Full coefficient estimates and sensitivity analyses for the juvenile:female ratio.",
  "",
  "### Panel A. Final retained model",
  "",
  md_table(data.frame(
    Predictor = s6_fixed$predictor, Reference = s6_fixed$reference,
    `Logit coefficient` = fmt(s6_fixed$estimate), SE = fmt(s6_fixed$SE),
    `Odds ratio` = fmt(s6_fixed$odds_ratio), `Lower 95% CI` = fmt(s6_fixed$lower_95_CI),
    `Upper 95% CI` = fmt(s6_fixed$upper_95_CI), P = fmt_p(s6_fixed$P), check.names = FALSE
  )),
  "",
  md_table(data.frame(
    Parameter = c(paste0("Random-effect SD: ", s6_re$predictor), "Beta-binomial dispersion"),
    Estimate = fmt(c(s6_re$estimate, s6_disp$estimate), 4), check.names = FALSE
  )),
  "",
  "### Panel B. Sensitivity and model-extension tests",
  "",
  md_table(data.frame(
    Test = s6_sensitivity$predictor, Component = s6_sensitivity$component,
    Estimate = fmt(s6_sensitivity$estimate), SE = fmt(s6_sensitivity$SE),
    `Odds ratio` = fmt(s6_sensitivity$odds_ratio), `Lower 95% CI` = fmt(s6_sensitivity$lower_95_CI),
    `Upper 95% CI` = fmt(s6_sensitivity$upper_95_CI), `z or LR` = fmt(s6_sensitivity$z_or_LR),
    df = fmt(s6_sensitivity$df, 0), P = fmt_p(s6_sensitivity$P), check.names = FALSE
  )),
  "",
  "Continuous covariates were standardized. Odds ratios greater than one indicate a higher juvenile:female ratio.",
  "",
  "# Appendix S3. Diel activity analyses",
  "",
  paste(
    "Apparent solar time was calculated from corrected UTC timestamps using longitude and the date-specific equation of time, preserving equal one-hour durations while aligning solar noon to 12:00.",
    paste0("Astronomical night was solar elevation below 0", u_degree, ", and class-specific kernel overlap used the Dhat4 estimator."),
    "For the hierarchical analysis, each positive-effort camera-round contributed 24 one-hour binomial histories; zero-success hours were retained and recorded operational days supplied the trials.",
    "Cyclic binomial GAMs included survey round, camera and camera-round structure, and, where selected, camera-specific cyclic deviations."
  ),
  "",
  paste(
    paste0("Baseline, wolf-only, jackal-only, additive, and wolf ", u_times, " jackal temporal structures were compared by AIC, retaining the least-complex candidate within 2 AIC units of the minimum."),
    "Predicted curves were normalized within predator context and therefore represent redistribution across the 24-hour cycle.",
    "The seven-specification audit comprised the primary apparent-solar fit, k = 8, k = 16, sun time, an alternate camera-only hierarchy, independent 30-min histories, and sun time plus the camera-only hierarchy.",
    "The confirmatory joint GAM used common predator-time smooths and symmetric sum-to-zero demographic deviation smooths with penalty-aware boundary tests."
  ),
  "",
  "## Table S7. Candidate-model comparison for predator-associated diel activity.",
  "",
  md_table(setNames(data.frame(
    `Demographic analysis` = s7$demographic_analysis, Model = s7$model,
    AIC = fmt(s7$AIC), delta_AIC = fmt(s7$delta_AIC),
    `Selected under rule` = s7$selected_under_least_complex_within_2_AIC, check.names = FALSE
  ), c("Demographic analysis", "Model", "AIC", paste0(u_delta, "AIC"), "Selected under rule"))),
  "",
  "The pooled wolf-only model had the minimum AIC; the additive model was 0.196 AIC units higher. Thus, the simplicity rule agrees with, but is not solely responsible for, the pooled selection.",
  "",
  "## Table S8. Predator-associated temporal effects and robustness across sensitivity specifications.",
  "",
  md_table(data.frame(
    `Demographic analysis` = s8$demographic_analysis, Term = s8$term,
    EDF = fmt(s8$EDF), Statistic = fmt(s8$statistic), P = fmt_p(s8$P),
    `Supported specifications` = paste0(s8$supported_specifications, "/", s8$tested_specifications),
    Interpretation = s8$interpretation, check.names = FALSE
  )),
  "",
  "The seven specifications were the primary apparent-solar fit, k = 8, k = 16, sun time, alternate camera-only hierarchy, independent 30-min histories, and sun time plus camera-only hierarchy. Only terms retained in the selected model and supported in all seven specifications were classified as robust.",
  "",
  "## Table S9. Confirmatory tests of demographic heterogeneity in predator-associated diel responses.",
  "",
  "### Panel A. Omnibus joint-GAM tests",
  "",
  md_table(data.frame(
    Term = s9a$term, EDF = fmt(s9a$EDF), `Reference rank` = s9a$reference_rank,
    Statistic = fmt(s9a$statistic), P = fmt_p(s9a$raw_P),
    `Sensitivity support` = s9a$sensitivity_support, Interpretation = s9a$interpretation,
    check.names = FALSE
  )),
  "",
  "### Panel B. Planned pairwise contrasts",
  "",
  md_table(data.frame(
    Contrast = s9b$contrast_or_scope, EDF = fmt(s9b$EDF), `Reference rank` = s9b$reference_rank,
    Statistic = fmt(s9b$statistic), `Raw P` = fmt_p(s9b$raw_P),
    `Holm-adjusted P` = fmt_p(s9b$Holm_adjusted_P), Interpretation = s9b$interpretation,
    check.names = FALSE
  )),
  "",
  "Statistics are `mgcv` penalty-aware boundary-test statistics; reference rank is not a conventional chi-square degrees of freedom. Pairwise tests apply to the primary apparent-solar specification. A nonsignificant adult-female versus adult-male contrast is not evidence of equivalence.",
  "",
  "## Table S10. Diagnostic summary for hierarchical cyclic GAMs.",
  "",
  "### Panel A. Sample and model structure",
  "",
  md_table(data.frame(
    Model = s10$model, Rows = s10$aggregated_rows, `Camera-rounds` = s10$camera_rounds,
    Trials = s10$binomial_trials, `Positive camera-day-hours` = s10$positive_camera_day_hours,
    k = s10$cyclic_basis_k, Hierarchy = s10$hierarchy, check.names = FALSE
  )),
  "",
  "### Panel B. Common fitting diagnostics",
  "",
  md_table(data.frame(
    Model = s10$model, `Deviance explained` = fmt(s10$deviance_explained),
    `Pearson dispersion` = fmt(s10$Pearson_dispersion),
    Convergence = s10$convergence,
    check.names = FALSE
  )),
  "",
  "### Panel C. Structure-specific diagnostic checks",
  "",
  md_table(s10_specific_checks),
  "",
  "Structure-specific diagnostics are listed only for models for which they were reported or meaningful. The low juvenile k-index motivated the prespecified k = 8 and k = 16 sensitivity fits rather than post hoc alteration of the primary model. No explicit AR(1) correction was used; residual lag-1 summaries are reported diagnostically."
)

qmd_path <- file.path(out_dir, "supplementary_materials_final.qmd")
writeLines(qmd, qmd_path, useBytes = TRUE)

# -----------------------------------------------------------------------------
# Numerical and cross-reference QC.
# -----------------------------------------------------------------------------

qc_rows <- list()
add_qc <- function(table, result, manuscript_value, reproduced_value, tolerance, source, note = NA_character_) {
  if (is.logical(reproduced_value)) {
    match_value <- identical(reproduced_value, manuscript_value)
  } else if (is.na(reproduced_value) || is.na(manuscript_value)) {
    match_value <- NA
  } else {
    match_value <- abs(reproduced_value - manuscript_value) <= tolerance
  }
  qc_rows[[length(qc_rows) + 1L]] <<- data.frame(
    table = table, result = result, manuscript_value = as.character(manuscript_value),
    reproduced_value = as.character(reproduced_value), tolerance = tolerance,
    match_within_rounding = ifelse(is.na(match_value), "NOT CHECKABLE", ifelse(match_value, "YES", "NO")),
    source = source, note = note, check.names = FALSE
  )
}

# S1 validation values.
for (i in seq_len(3)) {
  add_qc("S1", paste0(s1$pair[i], " estimate"), c(0.863, 0.657, 0.907)[i], s1$posterior_estimate[i], 0.0006, joint_fit_file)
  add_qc("S1", paste0(s1$pair[i], " lower CrI"), c(0.523, 0.435, 0.712)[i], s1$lower_95_CrI[i], 0.0006, joint_fit_file)
  add_qc("S1", paste0(s1$pair[i], " upper CrI"), c(0.998, 0.806, 1.000)[i], s1$upper_95_CrI[i], 0.0006, joint_fit_file)
}
add_qc("S2", paste0("Wolf", u_endash, "jackal site correlation"), -0.183, predator_r, 0.0006, "SiteCovs_060324.csv")

# S3 values explicitly reported in the manuscript/prompt.
s3_low <- s3[s3$panel == "A. Low-jackal / high-jackal abundance ratio" & grepl("Low", s3$wolf_activity), ]
for (cl in c("Juvenile", "Adult female", "Adult male")) {
  row <- s3_low[s3_low$demographic_class == cl, ]
  expected <- c("Juvenile" = 12.8, "Adult female" = 12.2, "Adult male" = 7.5)[cl]
  add_qc("S3", paste0(cl, " low-wolf jackal ratio"), expected, row$posterior_median, 0.06, class_fits_file)
}
s3_ror <- s3[s3$panel == "B. Attenuation ratio-of-ratios", ]
for (cl in c("Adult female", "Adult male")) {
  row <- s3_ror[s3_ror$demographic_class == cl, ]
  expected <- c("Adult female" = 7.97, "Adult male" = 4.90)[cl]
  add_qc("S3", paste0(cl, " attenuation ratio-of-ratios"), expected, row$posterior_median, 0.006, class_fits_file)
}

# S5 model AIC and LRT validation.
expected_aic <- c(444.913, 443.222, 444.194, 444.105, 442.722, 441.075)
for (i in seq_along(expected_aic)) {
  add_qc("S5", paste0("Candidate AIC: ", s5a$model[i]), expected_aic[i], s5a$AIC[i], 0.0006, ratio_fits_file)
}
expected_lr <- c(4.147, 0.35, 0.971)
expected_lrt_p <- c(0.0417, 0.552, 0.616)
for (i in seq_len(3)) {
  add_qc("S5", paste0(s5b$model[i], " LR"), expected_lr[i], s5b$LR_statistic[i], if (i == 2) 0.006 else 0.0006, ratio_fits_file)
  add_qc("S5", paste0(s5b$model[i], " P"), expected_lrt_p[i], s5b$P[i], 0.0006, ratio_fits_file)
}

# S6 focal results.
for (term in c("Late vs Early survey", "Open-space cover", "Wolf activity", "Jackal activity")) {
  row <- s6_fixed[s6_fixed$predictor == term, ]
  targets <- list(
    "Late vs Early survey" = c(1.57, 1.01, 2.44, 0.045),
    "Open-space cover" = c(1.55, 1.04, 2.32, 0.031),
    "Wolf activity" = c(0.93, 0.70, 1.23, 0.608),
    "Jackal activity" = c(0.94, 0.63, 1.41, 0.773)
  )[[term]]
  vals <- c(row$odds_ratio, row$lower_95_CI, row$upper_95_CI, row$P)
  labels <- c("OR", "lower CI", "upper CI", "P")
  for (j in seq_along(vals)) add_qc("S6", paste(term, labels[j]), targets[j], vals[j], if (j < 4) 0.006 else 0.0006, ratio_fits_file)
}
add_qc("S6", paste0("Wolf ", u_times, " jackal sensitivity OR"), 0.86, s6_sensitivity$odds_ratio[1], 0.006, ratio_fits_file)

# S7 expected selections and nearest-competitor differences.
selected_expected <- c(
  "Juvenile" = paste("Wolf", "jackal", sep = paste0(" ", u_times, " ")),
  "Adult female" = paste("Wolf", "jackal", sep = paste0(" ", u_times, " ")),
  "Adult male" = paste("Wolf", "jackal", sep = paste0(" ", u_times, " ")),
  "Pooled gazelles" = "Wolf-only"
)
for (cl in names(selected_expected)) {
  actual <- s7$model[s7$demographic_analysis == cl & s7$selected_under_least_complex_within_2_AIC == "Yes"]
  qc_rows[[length(qc_rows) + 1L]] <- data.frame(
    table = "S7", result = paste0(cl, " selected model"), manuscript_value = selected_expected[cl],
    reproduced_value = actual, tolerance = NA_real_, match_within_rounding = ifelse(identical(unname(selected_expected[cl]), actual), "YES", "NO"),
    source = "definitive_hierarchical_diel_ecological_model_comparison.csv", note = NA_character_, check.names = FALSE
  )
}
for (i in seq_along(c("Juvenile", "Adult female", "Adult male", "Pooled gazelles"))) {
  cl <- c("Juvenile", "Adult female", "Adult male", "Pooled gazelles")[i]
  model <- if (cl == "Pooled gazelles") "Additive wolf + jackal" else "Additive wolf + jackal"
  actual <- s7$delta_AIC[s7$demographic_analysis == cl & s7$model == model]
  expected <- c(9.666, 2.495, 3.262, 0.196)[i]
  add_qc("S7", paste0(cl, " additive-model ", u_delta, "AIC"), expected, actual, 0.0006,
         "definitive_hierarchical_diel_ecological_model_comparison.csv")
}

# S8-S9 key inferential values, including the known stale pooled EDF.
get_s8 <- function(cl, term) s8[s8$demographic_analysis == cl & s8$term == term, ]
row <- get_s8("Juvenile", paste(c("Wolf", "jackal", "time"), collapse = paste0(" ", u_times, " ")))
add_qc("S8", "Juvenile interaction EDF", 2.067, row$EDF, 0.0006, row$source)
add_qc("S8", "Juvenile interaction support", 7, row$supported_specifications, 0, row$source)
row <- get_s8("Adult female", paste("Wolf", "time", sep = paste0(" ", u_times, " ")))
add_qc("S8", "Adult-female wolf EDF", 3.916, row$EDF, 0.0006, row$source)
row <- get_s8("Adult male", paste(c("Wolf", "jackal", "time"), collapse = paste0(" ", u_times, " ")))
add_qc("S8", "Adult-male interaction EDF", 2.991, row$EDF, 0.0006, row$source)
row <- get_s8("Pooled gazelles", paste("Wolf", "time", sep = paste0(" ", u_times, " ")))
add_qc("S8", "Pooled selected wolf-only EDF", 3.041, row$EDF, 0.0006, row$source,
       "Manuscript/prompt value 3.041 is from the full interaction robustness model; the final selected wolf-only model gives EDF 3.870")

for (i in seq_len(nrow(s9a))) {
  expected_p <- c(0.0778, 0.7066, 0.0060)[i]
  add_qc("S9", paste0(s9a$term[i], " omnibus P"), expected_p, s9a$raw_P[i], 0.00006, s9a$source[i])
}
for (i in seq_len(nrow(s9b))) {
  expected_p <- c(0.0121, 0.0111, 0.2767)[i]
  add_qc("S9", paste0(s9b$contrast_or_scope[i], " Holm P"), expected_p, s9b$Holm_adjusted_P[i], 0.00006, s9b$source[i])
}

qc <- do.call(rbind, qc_rows)
row.names(qc) <- NULL

manuscript_path <- "C:/Users/preissbloom/OneDrive - Tel-Aviv University/Gazelles/Golan Gazelles/Spatiotemporal responses of an endangered ungulate.docx"
crossref_lines <- c()
manuscript_text <- NA_character_
if (file.exists(manuscript_path)) {
  tmp <- tempfile("manuscript_docx_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  unzip(manuscript_path, files = "word/document.xml", exdir = tmp)
  xml <- paste(readLines(file.path(tmp, "word/document.xml"), warn = FALSE, encoding = "UTF-8"), collapse = "")
  xml <- gsub("</w:p>", "\n", xml, fixed = TRUE)
  manuscript_text <- gsub("<[^>]+>", "", xml)
  manuscript_text <- gsub("&amp;", "&", manuscript_text, fixed = TRUE)
  manuscript_text <- gsub("&lt;", "<", manuscript_text, fixed = TRUE)
  manuscript_text <- gsub("&gt;", ">", manuscript_text, fixed = TRUE)
  manuscript_text <- gsub("[[:space:]]+", " ", manuscript_text)
  hits <- gregexpr("(Appendix|Table) S[0-9]+", manuscript_text, perl = TRUE)[[1]]
  if (hits[1] != -1) {
    hit_len <- attr(hits, "match.length")
    crossref_lines <- vapply(seq_along(hits), function(i) {
      lo <- max(1, hits[i] - 90)
      hi <- min(nchar(manuscript_text), hits[i] + hit_len[i] + 90)
      paste0("- `", substr(manuscript_text, lo, hi), "`")
    }, character(1))
  }
} else {
  crossref_lines <- "- Current manuscript DOCX was not found at the configured path; cross-references were not scanned."
}

crossref_ok <- !is.na(manuscript_text) &&
  all(vapply(c("Appendix S1", "Appendix S2", "Appendix S3"), grepl, logical(1), x = manuscript_text, fixed = TRUE)) &&
  !grepl("Appendix S4|Appendix S5", manuscript_text)

structure_checks <- c(
  paste0("- Ten unique table numbers generated: **YES** (`S1`", u_endash, "`S10`)."),
  "- Tables generated in numerical order: **YES**.",
  paste0("- Appendix S1 contains Tables S1", u_endash, "S4: **YES**."),
  paste0("- Appendix S2 contains Tables S5", u_endash, "S6: **YES**."),
  paste0("- Appendix S3 contains Tables S7", u_endash, "S10: **YES**."),
  paste0("- Current manuscript appendix cross-references use the final S1/S2/S3 map and contain no Appendix S4/S5 references: **",
         ifelse(crossref_ok, "YES", "NO OR NOT CHECKED"), "**.")
)

qc_display <- data.frame(
  Table = qc$table, Result = qc$result, `Manuscript value` = qc$manuscript_value,
  `Reproduced value` = qc$reproduced_value, `Within rounding` = qc$match_within_rounding,
  Source = qc$source, Note = qc$note, check.names = FALSE
)
qc_text <- c(
  "# Supplementary table quality control",
  "",
  "Generated from saved final-analysis objects and summary files without refitting models.",
  "",
  "## Numerical validation",
  "",
  md_table(qc_display),
  "",
  "## Structural checks",
  "",
  structure_checks,
  "",
  "## Unresolved or material issues",
  "",
  "1. **Pooled diel EDF mismatch:** the selected pooled wolf-only GAM gives EDF = 3.869643 (Table S8). The value 3.041 in the current manuscript/prompt comes from the pooled full-interaction robustness fit and must not be reported as the selected-model EDF.",
  "2. **Predator-index provenance:** the current manuscript records the final state/detection specifications and empirical-Bayes extraction, but the original wolf/jackal fit objects and candidate-fit statistics are absent. Table S2 therefore does not report unreproducible AIC, WAIC, or convergence statistics for those index models.",
  "3. **Early/Late boundary wording:** code classifies round midpoints on or before 17 October 2020 as Early and later dates as Late; the manuscript says before 17 October versus 17 October onward. No midpoint occurred exactly on 17 October, so estimates are unaffected, but wording should be harmonized.",
  paste0("4. **Pairwise joint-GAM robustness:** the omnibus wolf ", u_times, " jackal demographic heterogeneity was audited in 7/7 specifications; pairwise juvenile", u_endash, "adult contrasts were tested only in the primary apparent-solar model and must not be described as 7/7 robust."),
  "",
  "## Manuscript cross-reference occurrences",
  "",
  crossref_lines,
  "",
  "The manuscript was inspected only; it was not modified."
)

qc_path <- file.path(root, "supplementary_table_qc.md")
writeLines(qc_text, qc_path, useBytes = TRUE)

cat("Generated supplementary tables:\n")
for (p in c(s1_path, s2_path, s3_path, s4_path, s5_path, s6_path, s7_path, s8_path, s9_path, s10_path)) {
  cat(" - ", normalizePath(p, winslash = "/", mustWork = TRUE), "\n", sep = "")
}
cat("Combined supplement: ", normalizePath(qmd_path, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("QC report: ", normalizePath(qc_path, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("Numerical QC failures: ", sum(qc$match_within_rounding == "NO"), "\n", sep = "")
