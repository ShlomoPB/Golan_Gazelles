if (!requireNamespace("overlap", quietly = TRUE)) {
  stop("Package 'overlap' is required.")
}
if (!requireNamespace("activity", quietly = TRUE)) {
  stop("Package 'activity' is required.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required.")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required.")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  stop("Package 'readr' is required.")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  stop("Package 'tidyr' is required.")
}

library(overlap)
library(activity)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
arg_map <- list(
  event_file = "Gazella gazella_independent_events_30min.csv",
  site_cov_file = "SiteCovs_060324.csv",
  output_prefix = "gazelle_temporal_activity",
  timezone = "Asia/Jerusalem",
  overlap_boot = "499",
  activity_boot = "499",
  min_events_group = "20"
)

for (arg in args) {
  parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    key <- sub("^--", "", parts[1])
    arg_map[[key]] <- parts[2]
  }
}

event_file <- arg_map$event_file
site_cov_file <- arg_map$site_cov_file
output_prefix <- arg_map$output_prefix
timezone_local <- arg_map$timezone
overlap_boot <- as.integer(arg_map$overlap_boot)
activity_boot <- as.integer(arg_map$activity_boot)
min_events_group <- as.integer(arg_map$min_events_group)

pressure_vars <- c(
  "human.abundance",
  "Null.Wolf.Abundance",
  "null.jackal.abundance",
  "null.boar.abundance",
  "Cropland",
  "Built.and.Roads",
  "Firing.Zones",
  "Minefields",
  "Nature.Reserves",
  "Open.Spaces",
  "Orchards"
)

class_labels <- c(
  "juvenile" = "Fawn",
  "male" = "Male",
  "female" = "Female"
)

pretty_var <- c(
  "human.abundance" = "Human abundance",
  "Null.Wolf.Abundance" = "Wolf abundance",
  "null.jackal.abundance" = "Jackal abundance",
  "null.boar.abundance" = "Boar abundance",
  "Cropland" = "Cropland",
  "Built.and.Roads" = "Built and roads",
  "Firing.Zones" = "Firing zones",
  "Minefields" = "Minefields",
  "Nature.Reserves" = "Nature reserves",
  "Open.Spaces" = "Open spaces",
  "Orchards" = "Orchards"
)

to_radians <- function(x) {
  2 * pi * x / 24
}

to_hours <- function(x) {
  24 * x / (2 * pi)
}

pick_overlap_type <- function(n1, n2) {
  if (min(n1, n2) < 75) "Dhat1" else "Dhat4"
}

fmt_num <- function(x, digits = 3) formatC(x, digits = digits, format = "f")
fmt_ci <- function(mean, lwr, upr, digits = 3) {
  paste0(fmt_num(mean, digits), " [", fmt_num(lwr, digits), ", ", fmt_num(upr, digits), "]")
}

extract_hours <- function(x, tz_local) {
  ts_utc <- as.POSIXct(x, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  ts_local <- as.POSIXct(format(ts_utc, tz = tz_local, usetz = TRUE), tz = tz_local)
  as.numeric(format(ts_local, "%H")) +
    as.numeric(format(ts_local, "%M")) / 60 +
    as.numeric(format(ts_local, "%S")) / 3600
}

extract_density_df <- function(fit, label_cols) {
  pdf_df <- as.data.frame(fit@pdf)
  data.frame(
    hour = to_hours(pdf_df$x),
    density = pdf_df$y,
    label_cols,
    stringsAsFactors = FALSE
  )
}

calc_activity_pair <- function(times1, times2, label1, label2, activity_reps, overlap_reps) {
  n1 <- length(times1)
  n2 <- length(times2)
  out <- data.frame(
    group1 = label1,
    group2 = label2,
    n1 = n1,
    n2 = n2,
    overlap_type = pick_overlap_type(n1, n2),
    overlap = NA_real_,
    overlap_lwr = NA_real_,
    overlap_upr = NA_real_,
    activity1 = NA_real_,
    activity2 = NA_real_,
    activity_diff_p = NA_real_,
    stringsAsFactors = FALSE
  )
  if (n1 < min_events_group || n2 < min_events_group) {
    return(out)
  }

  overlap_type <- out$overlap_type[1]
  t0 <- unname(overlapEst(times1, times2, type = overlap_type))
  bt <- bootEst(
    resample(times1, overlap_reps),
    resample(times2, overlap_reps),
    type = overlap_type
  )
  ci <- bootCI(t0, bt)
  ci_use <- ci["norm0", ]
  if (any(is.na(ci_use))) {
    ci_use <- ci["norm", ]
  }

  fit1 <- fitact(times1, reps = activity_reps, sample = "data", show = FALSE)
  fit2 <- fitact(times2, reps = activity_reps, sample = "data", show = FALSE)
  cmp <- compareCkern(fit1, fit2, reps = activity_reps)

  out$overlap <- t0
  out$overlap_lwr <- unname(ci_use[1])
  out$overlap_upr <- unname(ci_use[2])
  out$activity1 <- unname(fit1@act[1])
  out$activity2 <- unname(fit2@act[1])
  out$activity_diff_p <- unname(cmp["pNull"])
  attr(out, "fit1") <- fit1
  attr(out, "fit2") <- fit2
  out
}

events <- read_csv(event_file, show_col_types = FALSE) %>%
  mutate(
    camera_name = as.integer(camera_name),
    gazelle_class = as.character(gazelle_class),
    class_label = unname(class_labels[gazelle_class]),
    local_hour = extract_hours(event_time, timezone_local),
    time_rad = to_radians(local_hour)
  )

site_covs <- read_csv(site_cov_file, show_col_types = FALSE) %>%
  rename(camera_name = Camera.ID)

missing_vars <- setdiff(pressure_vars, names(site_covs))
if (length(missing_vars) > 0) {
  stop("Missing site covariates: ", paste(missing_vars, collapse = ", "))
}

events <- events %>%
  left_join(site_covs %>% select(camera_name, all_of(pressure_vars)), by = "camera_name")

pair_mat <- t(combn(sort(unique(events$gazelle_class)), 2))
overall_results <- vector("list", nrow(pair_mat))
overall_density <- list()

for (i in seq_len(nrow(pair_mat))) {
  c1 <- pair_mat[i, 1]
  c2 <- pair_mat[i, 2]
  t1 <- events$time_rad[events$gazelle_class == c1]
  t2 <- events$time_rad[events$gazelle_class == c2]
  pair_res <- calc_activity_pair(
    t1, t2,
    unname(class_labels[c1]),
    unname(class_labels[c2]),
    activity_boot,
    overlap_boot
  )
  pair_res$class1 <- c1
  pair_res$class2 <- c2
  overall_results[[i]] <- pair_res

  fit1 <- attr(pair_res, "fit1")
  fit2 <- attr(pair_res, "fit2")
  if (!is.null(fit1) && !is.null(fit2)) {
    overall_density[[length(overall_density) + 1]] <- extract_density_df(
      fit1, list(class = unname(class_labels[c1]), comparison = paste(unname(class_labels[c1]), "vs", unname(class_labels[c2])))
    )
    overall_density[[length(overall_density) + 1]] <- extract_density_df(
      fit2, list(class = unname(class_labels[c2]), comparison = paste(unname(class_labels[c1]), "vs", unname(class_labels[c2])))
    )
  }
}

overall_results <- bind_rows(overall_results) %>%
  mutate(
    comparison = paste(group1, "vs", group2),
    overlap_ci = fmt_ci(overlap, overlap_lwr, overlap_upr),
    activity1_fmt = fmt_num(activity1),
    activity2_fmt = fmt_num(activity2),
    activity_diff_p_fmt = fmt_num(activity_diff_p)
  )

overall_density_df <- bind_rows(overall_density)

high_low_results <- list()
high_low_density <- list()
idx <- 1

for (cls in sort(unique(events$gazelle_class))) {
  class_events <- events %>% filter(gazelle_class == cls)
  for (var in pressure_vars) {
    q1 <- unname(quantile(site_covs[[var]], 0.25, na.rm = TRUE))
    q3 <- unname(quantile(site_covs[[var]], 0.75, na.rm = TRUE))
    low_cams <- site_covs$camera_name[site_covs[[var]] <= q1]
    high_cams <- site_covs$camera_name[site_covs[[var]] >= q3]
    low_times <- class_events$time_rad[class_events$camera_name %in% low_cams]
    high_times <- class_events$time_rad[class_events$camera_name %in% high_cams]

    pair_res <- calc_activity_pair(
      low_times, high_times,
      "Low quartile", "High quartile",
      activity_boot,
      overlap_boot
    )
    pair_res$class <- unname(class_labels[cls])
    pair_res$class_code <- cls
    pair_res$covariate <- var
    pair_res$covariate_label <- unname(pretty_var[var])
    pair_res$q1 <- q1
    pair_res$q3 <- q3
    pair_res$n_cameras_low <- length(unique(low_cams))
    pair_res$n_cameras_high <- length(unique(high_cams))
    high_low_results[[idx]] <- pair_res
    idx <- idx + 1

    fit1 <- attr(pair_res, "fit1")
    fit2 <- attr(pair_res, "fit2")
    if (!is.null(fit1) && !is.null(fit2)) {
      high_low_density[[length(high_low_density) + 1]] <- extract_density_df(
        fit1, list(class = unname(class_labels[cls]), covariate = unname(pretty_var[var]), group = "Low quartile")
      )
      high_low_density[[length(high_low_density) + 1]] <- extract_density_df(
        fit2, list(class = unname(class_labels[cls]), covariate = unname(pretty_var[var]), group = "High quartile")
      )
    }
  }
}

high_low_results <- bind_rows(high_low_results) %>%
  mutate(
    overlap_ci = ifelse(is.na(overlap), NA_character_, fmt_ci(overlap, overlap_lwr, overlap_upr)),
    activity_low_fmt = ifelse(is.na(activity1), NA_character_, fmt_num(activity1)),
    activity_high_fmt = ifelse(is.na(activity2), NA_character_, fmt_num(activity2)),
    activity_diff_p_fmt = ifelse(is.na(activity_diff_p), NA_character_, fmt_num(activity_diff_p)),
    q1_fmt = fmt_num(q1),
    q3_fmt = fmt_num(q3)
  )

high_low_density_df <- bind_rows(high_low_density)

write_csv(overall_results, paste0(output_prefix, "_overall_class_comparisons.csv"))
write_csv(high_low_results, paste0(output_prefix, "_high_low_pressure_comparisons.csv"))

overall_plot <- ggplot(overall_density_df, aes(x = hour, y = density, color = class)) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous(breaks = seq(0, 24, by = 4), limits = c(0, 24)) +
  labs(
    x = paste0("Local time (", timezone_local, ")"),
    y = "Kernel density",
    color = NULL,
    title = "Temporal activity curves by gazelle class"
  ) +
  facet_wrap(~ comparison, ncol = 1) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(paste0(output_prefix, "_overall_class_activity_curves.png"), overall_plot, width = 8.5, height = 8, dpi = 300)
ggsave(paste0(output_prefix, "_overall_class_activity_curves.pdf"), overall_plot, width = 8.5, height = 8)

if (nrow(high_low_density_df) > 0) {
  hl_plot <- ggplot(high_low_density_df, aes(x = hour, y = density, color = group)) +
    geom_line(linewidth = 0.8) +
    scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(0, 24)) +
    scale_color_manual(values = c("Low quartile" = "#1b9e77", "High quartile" = "#d95f02")) +
    labs(
      x = paste0("Local time (", timezone_local, ")"),
      y = "Kernel density",
      color = NULL,
      title = "Temporal activity by class across high- and low-pressure camera groups"
    ) +
    facet_grid(class ~ covariate) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
  ggsave(paste0(output_prefix, "_high_low_pressure_activity_curves.png"), hl_plot, width = 16, height = 8.5, dpi = 300)
  ggsave(paste0(output_prefix, "_high_low_pressure_activity_curves.pdf"), hl_plot, width = 16, height = 8.5)
}

md_lines <- c(
  "# Gazelle Temporal Activity Analysis",
  "",
  "## Overview",
  "",
  paste0("- Event file: `", event_file, "`"),
  paste0("- Site covariates: `", site_cov_file, "`"),
  paste0("- Local timezone used for clock time conversion: `", timezone_local, "`"),
  "- Each independent event timestamp was treated as one detection; `event_count` was not used as a weight.",
  paste0("- High and low pressure groups were defined as the top and bottom quartiles of each site covariate across cameras."),
  paste0("- Minimum event threshold per group for formal comparison: `", min_events_group, "`"),
  "",
  "## Overall class comparisons",
  "",
  "| Comparison | n1 | n2 | Overlap estimator | Overlap [95% CI] | Activity1 | Activity2 | Kernel comparison p |",
  "|---|---:|---:|---|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(overall_results))) {
  md_lines <- c(
    md_lines,
    paste0(
      "| `", overall_results$comparison[i], "` | `", overall_results$n1[i], "` | `", overall_results$n2[i],
      "` | `", overall_results$overlap_type[i], "` | `", overall_results$overlap_ci[i],
      "` | `", overall_results$activity1_fmt[i], "` | `", overall_results$activity2_fmt[i],
      "` | `", overall_results$activity_diff_p_fmt[i], "` |"
    )
  )
}

md_lines <- c(
  md_lines,
  "",
  "## High vs low pressure comparisons",
  "",
  "| Class | Covariate | Low cameras | High cameras | Low events | High events | Overlap [95% CI] | Kernel comparison p |",
  "|---|---|---:|---:|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(high_low_results))) {
  md_lines <- c(
    md_lines,
    paste0(
      "| `", high_low_results$class[i], "` | `", high_low_results$covariate_label[i], "` | `",
      high_low_results$n_cameras_low[i], "` | `", high_low_results$n_cameras_high[i], "` | `",
      high_low_results$n1[i], "` | `", high_low_results$n2[i], "` | `",
      ifelse(is.na(high_low_results$overlap_ci[i]), "NA", high_low_results$overlap_ci[i]),
      "` | `", ifelse(is.na(high_low_results$activity_diff_p_fmt[i]), "NA", high_low_results$activity_diff_p_fmt[i]), "` |"
    )
  )
}

writeLines(md_lines, paste0(output_prefix, "_results.md"))

cat("Saved temporal analysis outputs with prefix:", output_prefix, "\n")
