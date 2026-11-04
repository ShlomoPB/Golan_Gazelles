#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(emmeans)
  library(ggplot2)
  library(glmmTMB)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
arg_map <- list(
  event_file = "Gazella gazella_independent_events_30min_analysis_ready.csv",
  obs_cov_file = "ObsCovs_050524.csv",
  camera_context_file = "publication_predator_context_class_diel_2x2_median_camera_categories.csv",
  output_prefix = "publication_predator_context_nocturnality_median",
  local_timezone = "Asia/Jerusalem",
  formal_minimum_events = "20",
  formal_minimum_cameras = "5"
)
for (arg in args) {
  parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    key <- sub("^--", "", parts[1])
    if (!key %in% names(arg_map)) stop("Unknown argument: ", key)
    arg_map[[key]] <- parts[2]
  }
}

event_file <- arg_map$event_file
obs_cov_file <- arg_map$obs_cov_file
camera_context_file <- arg_map$camera_context_file
output_prefix <- arg_map$output_prefix
local_timezone <- arg_map$local_timezone
formal_minimum_events <- as.integer(arg_map$formal_minimum_events)
formal_minimum_cameras <- as.integer(arg_map$formal_minimum_cameras)

required_files <- c(event_file, obs_cov_file, camera_context_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

class_levels <- c("Fawn", "Female", "Male")
context_levels <- c(
  "Low both",
  "High jackal / Low wolf",
  "Low jackal / High wolf",
  "High both"
)
class_labels <- c(
  juvenile = "Fawn",
  female = "Female",
  male = "Male"
)
context_colors <- c(
  "Low both" = "#188977",
  "High jackal / Low wolf" = "#B12B32",
  "Low jackal / High wolf" = "#173B63",
  "High both" = "#8A5A9E"
)

equation_of_time <- function(day_of_year) {
  b <- 2 * pi * (day_of_year - 81) / 364
  9.87 * sin(2 * b) - 7.53 * cos(b) - 1.5 * sin(b)
}

parse_local_camera_time <- function(x) {
  as.POSIXct(
    sub("Z$", "", x),
    format = "%Y-%m-%dT%H:%M:%S",
    tz = local_timezone
  )
}

apparent_solar_hour <- function(timestamp, longitude) {
  day_of_year <- as.integer(format(timestamp, "%j", tz = "UTC"))
  utc_hour <- as.numeric(format(timestamp, "%H", tz = "UTC")) +
    as.numeric(format(timestamp, "%M", tz = "UTC")) / 60 +
    as.numeric(format(timestamp, "%S", tz = "UTC")) / 3600
  (utc_hour + longitude / 15 + equation_of_time(day_of_year) / 60) %% 24
}

fmt <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", fmt(x, 3)))
}

camera_context <- read_csv(camera_context_file, show_col_types = FALSE) %>%
  transmute(
    camera_name = as.integer(camera_name),
    predator_context = factor(as.character(predator_context), levels = context_levels)
  )

if (nrow(camera_context) != 60L || n_distinct(camera_context$camera_name) != 60L) {
  stop("Expected one predator-context assignment for each of 60 cameras.")
}
if (any(is.na(camera_context$predator_context))) {
  stop("All cameras must have one of the four median-split predator contexts.")
}

obs <- read.csv(
  obs_cov_file,
  na.strings = c("NA", ""),
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  transmute(
    camera_name = as.integer(`Camera Number`),
    round = as.integer(`Round Number`),
    start_date = as.Date(`Start date`, format = "%d/%m/%Y"),
    end_date = as.Date(`End Date`, format = "%d/%m/%Y")
  ) %>%
  filter(!is.na(camera_name), !is.na(start_date), !is.na(end_date))

raw_events <- read_csv(
  event_file,
  col_types = cols(event_time = col_character()),
  show_col_types = FALSE
) %>%
  mutate(
    event_row_id = row_number(),
    camera_name = as.integer(camera_name),
    class = unname(class_labels[as.character(gazelle_class)]),
    timestamp = parse_local_camera_time(event_time),
    event_date = as.Date(format(timestamp, tz = local_timezone)),
    solar_hour = apparent_solar_hour(timestamp, camelot_gps_longitude),
    # Procko et al. (2023): 0 = solar noon and 12 = solar midnight.
    nocturnality_hours = abs(solar_hour - 12)
  )

active_matches <- raw_events %>%
  select(event_row_id, camera_name, event_date) %>%
  inner_join(obs, by = "camera_name", relationship = "many-to-many") %>%
  filter(event_date >= start_date, event_date <= end_date) %>%
  group_by(event_row_id) %>%
  slice_max(start_date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(event_row_id)

events <- raw_events %>%
  inner_join(active_matches, by = "event_row_id") %>%
  left_join(camera_context, by = "camera_name") %>%
  filter(
    class %in% class_levels,
    predator_context %in% context_levels,
    is.finite(nocturnality_hours)
  ) %>%
  mutate(
    class = factor(class, levels = class_levels),
    predator_context = factor(predator_context, levels = context_levels),
    camera_name = factor(camera_name)
  )

expected_counts <- c(Fawn = 180L, Female = 1477L, Male = 2268L)
observed_counts <- table(events$class)
if (
  nrow(events) != sum(expected_counts) ||
    !all(as.integer(observed_counts[names(expected_counts)]) == expected_counts)
) {
  stop("Validated event counts do not match 180 fawn, 1477 female, and 2268 male events.")
}

context_summary <- events %>%
  group_by(class, predator_context) %>%
  summarise(
    n_events = n(),
    n_cameras_with_events = n_distinct(camera_name),
    mean_nocturnality_hours = mean(nocturnality_hours),
    sd_nocturnality_hours = sd(nocturnality_hours),
    median_nocturnality_hours = median(nocturnality_hours),
    .groups = "drop"
  ) %>%
  complete(
    class = factor(class_levels, levels = class_levels),
    predator_context = factor(context_levels, levels = context_levels),
    fill = list(n_events = 0L, n_cameras_with_events = 0L)
  ) %>%
  mutate(
    formal_test_eligible =
      n_events >= formal_minimum_events &
      n_cameras_with_events >= formal_minimum_cameras
  ) %>%
  arrange(class, predator_context)

model_tests <- list()
estimated_means <- list()
pairwise_contrasts <- list()
model_diagnostics <- list()
model_objects <- list()

for (class_i in class_levels) {
  eligible_contexts <- context_summary %>%
    filter(class == class_i, formal_test_eligible) %>%
    pull(predator_context) %>%
    as.character()

  model_data <- events %>%
    filter(class == class_i, as.character(predator_context) %in% eligible_contexts) %>%
    droplevels()

  if (length(eligible_contexts) < 2L) {
    stop("Fewer than two eligible contexts for ", class_i, ".")
  }

  null_model <- glmmTMB(
    nocturnality_hours ~ 1 + (1 | camera_name),
    data = model_data,
    family = gaussian(),
    REML = FALSE
  )
  full_model <- glmmTMB(
    nocturnality_hours ~ predator_context + (1 | camera_name),
    data = model_data,
    family = gaussian(),
    REML = FALSE
  )
  model_objects[[class_i]] <- list(null = null_model, context = full_model)

  lrt <- as.data.frame(anova(null_model, full_model))
  model_tests[[class_i]] <- tibble(
    class = class_i,
    eligible_contexts = length(eligible_contexts),
    events_used = nrow(model_data),
    cameras_used = n_distinct(model_data$camera_name),
    likelihood_ratio_chisq = lrt$Chisq[2],
    df = lrt$`Chi Df`[2],
    p_value = lrt$`Pr(>Chisq)`[2]
  )

  emm <- emmeans(full_model, ~ predator_context)
  emm_df <- as.data.frame(summary(emm, infer = c(TRUE, FALSE)))
  emm_lower_name <- intersect(c("lower.CL", "asymp.LCL"), names(emm_df))[1]
  emm_upper_name <- intersect(c("upper.CL", "asymp.UCL"), names(emm_df))[1]
  estimated_means[[class_i]] <- tibble(
    class = class_i,
    predator_context = as.character(emm_df$predator_context),
    estimated_mean_nocturnality_hours = emm_df$emmean,
    standard_error = emm_df$SE,
    df = emm_df$df,
    conf_low = emm_df[[emm_lower_name]],
    conf_high = emm_df[[emm_upper_name]]
  )

  contrast_df <- as.data.frame(summary(
    pairs(emm, adjust = "none"),
    infer = c(TRUE, TRUE),
    adjust = "none"
  ))
  contrast_lower_name <- intersect(c("lower.CL", "asymp.LCL"), names(contrast_df))[1]
  contrast_upper_name <- intersect(c("upper.CL", "asymp.UCL"), names(contrast_df))[1]
  contrast_parts <- strsplit(as.character(contrast_df$contrast), " - ", fixed = TRUE)
  clean_context <- function(x) sub("^\\((.*)\\)$", "\\1", x)
  statistic_name <- intersect(c("t.ratio", "z.ratio"), names(contrast_df))[1]
  pairwise_contrasts[[class_i]] <- tibble(
    class = class_i,
    context_a = clean_context(vapply(contrast_parts, `[`, character(1), 1)),
    context_b = clean_context(vapply(contrast_parts, `[`, character(1), 2)),
    contrast = as.character(contrast_df$contrast),
    delta_hours_a_minus_b = contrast_df$estimate,
    standard_error = contrast_df$SE,
    df = contrast_df$df,
    conf_low = contrast_df[[contrast_lower_name]],
    conf_high = contrast_df[[contrast_upper_name]],
    statistic = contrast_df[[statistic_name]],
    p_value = contrast_df$p.value
  )

  camera_variance <- VarCorr(full_model)$cond$camera_name
  model_diagnostics[[class_i]] <- tibble(
    class = class_i,
    optimizer_convergence_code = full_model$fit$convergence,
    positive_definite_hessian = isTRUE(full_model$sdr$pdHess),
    camera_random_intercept_sd = unname(attr(camera_variance, "stddev")[1]),
    residual_sd = sigma(full_model),
    AIC = AIC(full_model)
  )
}

model_tests <- bind_rows(model_tests)
estimated_means <- bind_rows(estimated_means)
pairwise_contrasts <- bind_rows(pairwise_contrasts)
model_diagnostics <- bind_rows(model_diagnostics)

write_csv(context_summary, paste0(output_prefix, "_context_summary.csv"))
write_csv(model_tests, paste0(output_prefix, "_omnibus_tests.csv"))
write_csv(estimated_means, paste0(output_prefix, "_estimated_means.csv"))
write_csv(pairwise_contrasts, paste0(output_prefix, "_pairwise_contrasts.csv"))
write_csv(model_diagnostics, paste0(output_prefix, "_model_diagnostics.csv"))
saveRDS(model_objects, paste0(output_prefix, "_models.rds"))

sparse_plot_rows <- context_summary %>%
  filter(!formal_test_eligible) %>%
  transmute(
    class = as.character(class),
    predator_context = as.character(predator_context),
    estimated_mean_nocturnality_hours = mean_nocturnality_hours,
    conf_low = NA_real_,
    conf_high = NA_real_,
    model_estimated = FALSE,
    annotation = paste0("data deficient (n=", n_events, ")")
  )

plot_data <- estimated_means %>%
  transmute(
    class,
    predator_context,
    estimated_mean_nocturnality_hours,
    conf_low,
    conf_high,
    model_estimated = TRUE,
    annotation = NA_character_
  ) %>%
  bind_rows(sparse_plot_rows) %>%
  mutate(
    class = factor(class, levels = class_levels),
    predator_context = factor(predator_context, levels = rev(context_levels))
  )

nocturnality_plot <- ggplot(
  filter(plot_data, model_estimated),
  aes(
    x = estimated_mean_nocturnality_hours,
    y = predator_context,
    color = predator_context
  )
) +
  geom_segment(
    aes(x = conf_low, xend = conf_high, yend = predator_context),
    linewidth = 0.8
  ) +
  geom_point(size = 3.0) +
  geom_point(
    data = filter(plot_data, !model_estimated),
    aes(x = estimated_mean_nocturnality_hours, y = predator_context),
    shape = 4,
    stroke = 1.0,
    size = 3.2
  ) +
  geom_text(
    data = filter(plot_data, !model_estimated),
    aes(
      x = estimated_mean_nocturnality_hours,
      y = predator_context,
      label = annotation
    ),
    hjust = -0.12,
    size = 3.0,
    color = "#4E585B"
  ) +
  facet_wrap(~class, ncol = 1) +
  scale_color_manual(values = context_colors, guide = "none") +
  scale_x_continuous(
    limits = c(0, 12),
    breaks = seq(0, 12, by = 2),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = "Mean hours from apparent solar noon (95% CI)",
    y = NULL,
    title = "Gazelle nocturnality across wolf-jackal contexts",
    subtitle = "Higher values indicate detections closer to solar midnight",
    caption = paste0(
      "Points and intervals are estimated marginal means from Gaussian mixed models with camera random intercepts.\n",
      "The high-jackal/low-wolf fawn group is shown descriptively and excluded from formal inference."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(color = "#304044"),
    axis.title.x = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "#E8E1D4", color = NA),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "#4C5A5E"),
    plot.caption = element_text(color = "#5F696B", hjust = 0)
  )

ggsave(
  paste0(output_prefix, ".png"),
  nocturnality_plot,
  width = 9.5,
  height = 8.5,
  dpi = 320
)
ggsave(
  paste0(output_prefix, ".pdf"),
  nocturnality_plot,
  width = 9.5,
  height = 8.5
)

summary_rows <- vapply(seq_len(nrow(context_summary)), function(i) {
  row <- context_summary[i, ]
  paste0(
    "| ", row$class,
    " | ", row$predator_context,
    " | ", row$n_events,
    " | ", row$n_cameras_with_events,
    " | ", fmt(row$mean_nocturnality_hours),
    " | ", fmt(row$median_nocturnality_hours),
    " | ", ifelse(row$formal_test_eligible, "Yes", "No"), " |"
  )
}, character(1))

omnibus_rows <- vapply(seq_len(nrow(model_tests)), function(i) {
  row <- model_tests[i, ]
  paste0(
    "| ", row$class,
    " | ", fmt(row$likelihood_ratio_chisq),
    " | ", row$df,
    " | ", fmt_p(row$p_value), " |"
  )
}, character(1))

contrast_rows <- vapply(seq_len(nrow(pairwise_contrasts)), function(i) {
  row <- pairwise_contrasts[i, ]
  paste0(
    "| ", row$class,
    " | ", row$context_a, " - ", row$context_b,
    " | ", fmt(row$delta_hours_a_minus_b),
    " [", fmt(row$conf_low), ", ", fmt(row$conf_high), "]",
    " | ", fmt_p(row$p_value), " |"
  )
}, character(1))

supported_contrasts <- pairwise_contrasts %>%
  filter(p_value < 0.05) %>%
  mutate(
    result = paste0(
      class, ": ", context_a, " - ", context_b,
      " = ", fmt(delta_hours_a_minus_b), " h (95% CI ",
      fmt(conf_low), " to ", fmt(conf_high), "; P = ", fmt_p(p_value), ")"
    )
  ) %>%
  pull(result)

md <- c(
  "# Gazelle nocturnality across median-split wolf-jackal contexts",
  "",
  "## Methods",
  "",
  paste0(
    "Following Procko et al. (2023), we calculated nocturnality for each independent detection as the absolute decimal-hour distance from apparent solar noon (range 0-12 h; 0 = solar noon and 12 = solar midnight). We fitted a separate Gaussian mixed model for each gazelle class, with the four-level wolf-jackal context as a fixed effect and camera as a random intercept. We tested the overall context effect by likelihood-ratio comparison with an intercept-only mixed model and report unadjusted pairwise contrasts from estimated marginal means. Contexts with fewer than ",
    formal_minimum_events, " events or ", formal_minimum_cameras,
    " cameras with detections were excluded from formal inference."
  ),
  "",
  "Procko, M., R. Naidoo, V. LeMay, and A. C. Burton. 2023. Human presence and infrastructure impact wildlife nocturnality differently across an assemblage of mammalian species. PLOS ONE 18:e0286131. https://doi.org/10.1371/journal.pone.0286131",
  "",
  "## Context summaries",
  "",
  "| Class | Predator context | Events | Cameras with events | Mean nocturnality (h) | Median (h) | Formal inference |",
  "|---|---|---:|---:|---:|---:|---:|",
  summary_rows,
  "",
  "## Omnibus context tests",
  "",
  "| Class | Likelihood-ratio chi-square | df | P |",
  "|---|---:|---:|---:|",
  omnibus_rows,
  "",
  "## Unadjusted pairwise contrasts",
  "",
  "Positive differences indicate that the first context was more nocturnal than the second.",
  "",
  "| Class | Contrast | Difference [95% CI] (h) | P |",
  "|---|---|---:|---:|",
  contrast_rows,
  "",
  "## Pairwise contrasts with raw P < 0.05",
  "",
  if (length(supported_contrasts)) {
    paste0("- ", supported_contrasts)
  } else {
    "No pairwise nocturnality contrast had an unadjusted P < 0.05."
  },
  "Because these are unadjusted post hoc contrasts, they should be interpreted in relation to the class-specific omnibus test rather than as stand-alone evidence.",
  "",
  "## Interpretation constraint",
  "",
  "Nocturnality condenses each diel distribution to distance from solar noon. It can identify shifts toward night but cannot distinguish morning from evening redistribution. It should therefore be interpreted alongside the circular kernel curves and compareCkern tests. Predator contexts are observational site categories and do not establish predator-caused behavioral change."
)
writeLines(md, paste0(output_prefix, "_results.md"), useBytes = TRUE)

writeLines(
  c(
    R.version.string,
    paste0("glmmTMB ", packageVersion("glmmTMB")),
    paste0("emmeans ", packageVersion("emmeans")),
    paste0("ggplot2 ", packageVersion("ggplot2")),
    paste0("dplyr ", packageVersion("dplyr"))
  ),
  paste0(output_prefix, "_session_info.txt")
)

cat("Saved predator-context nocturnality outputs with prefix:", output_prefix, "\n")
