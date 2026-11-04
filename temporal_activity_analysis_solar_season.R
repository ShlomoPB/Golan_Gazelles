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

library(overlap)
library(activity)
library(ggplot2)
library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
arg_map <- list(
  event_file = "Gazella gazella_independent_events_30min.csv",
  site_cov_file = "SiteCovs_060324.csv",
  output_prefix = "gazelle_temporal_activity_solar_season",
  overlap_boot = "149",
  activity_boot = "149",
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

fmt_num <- function(x, digits = 3) formatC(x, digits = digits, format = "f")
fmt_ci <- function(mean, lwr, upr, digits = 3) {
  paste0(fmt_num(mean, digits), " [", fmt_num(lwr, digits), ", ", fmt_num(upr, digits), "]")
}
fmt_p_label <- function(p) {
  ifelse(
    is.na(p),
    "p = NA",
    ifelse(p < 0.001, "p < 0.001", paste0("p = ", fmt_num(p)))
  )
}
fmt_sig_label <- function(p) {
  ifelse(
    is.na(p),
    "NA",
    ifelse(
      p < 0.001,
      "***",
      ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))
    )
  )
}

equation_of_time <- function(doy) {
  b <- 2 * pi * (doy - 81) / 364
  9.87 * sin(2 * b) - 7.53 * cos(b) - 1.5 * sin(b)
}

to_solar_hour <- function(ts_utc, longitude) {
  doy <- as.integer(format(ts_utc, "%j"))
  utc_hour <- as.numeric(format(ts_utc, "%H")) +
    as.numeric(format(ts_utc, "%M")) / 60 +
    as.numeric(format(ts_utc, "%S")) / 3600
  solar_hour <- utc_hour + longitude / 15 + equation_of_time(doy) / 60
  solar_hour %% 24
}

to_radians <- function(x) 2 * pi * x / 24
to_hours <- function(x) 24 * x / (2 * pi)

pick_overlap_type <- function(n1, n2) {
  if (min(n1, n2) < 75) "Dhat1" else "Dhat4"
}

extract_density_df <- function(fit, extra) {
  pdf_df <- as.data.frame(fit@pdf)
  data.frame(
    hour = to_hours(pdf_df$x),
    density = pdf_df$y,
    extra,
    stringsAsFactors = FALSE
  )
}

calc_pair <- function(times1, times2, label1, label2, min_events, overlap_reps, activity_reps) {
  n1 <- length(times1)
  n2 <- length(times2)
  res <- data.frame(
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
  if (n1 < min_events || n2 < min_events) {
    return(res)
  }

  overlap_type <- res$overlap_type[1]
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

  res$overlap <- t0
  res$overlap_lwr <- unname(ci_use[1])
  res$overlap_upr <- unname(ci_use[2])
  res$activity1 <- unname(fit1@act[1])
  res$activity2 <- unname(fit2@act[1])
  res$activity_diff_p <- unname(cmp["pNull"])
  attr(res, "fit1") <- fit1
  attr(res, "fit2") <- fit2
  res
}

shade_df <- data.frame(
  xmin = c(0, 5, 7, 17, 19),
  xmax = c(5, 7, 17, 19, 24),
  light_phase = c("Night", "Dawn", "Day", "Dusk", "Night")
)

events <- read_csv(event_file, show_col_types = FALSE) %>%
  mutate(
    camera_name = as.integer(camera_name),
    gazelle_class = as.character(gazelle_class),
    class_label = unname(class_labels[gazelle_class]),
    ts_utc = as.POSIXct(event_time, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    event_date = as.Date(ts_utc, tz = "UTC"),
    solar_hour = to_solar_hour(ts_utc, camelot_gps_longitude),
    time_rad = to_radians(solar_hour)
  )

site_covs <- read_csv(site_cov_file, show_col_types = FALSE) %>%
  rename(camera_name = Camera.ID)

missing_vars <- setdiff(pressure_vars, names(site_covs))
if (length(missing_vars) > 0) {
  stop("Missing site covariates: ", paste(missing_vars, collapse = ", "))
}

events <- events %>%
  left_join(site_covs %>% select(camera_name, all_of(pressure_vars)), by = "camera_name")

split_time <- stats::median(events$ts_utc)
events <- events %>%
  mutate(
    season_period = ifelse(ts_utc <= split_time, "Early (hotter)", "Late (mild-cold)")
  )

## Overall class comparisons by season and pooled --------------------------------
pair_mat <- t(combn(sort(unique(events$gazelle_class)), 2))
overall_list <- list()
overall_density <- list()
idx <- 1

for (period in c("All detections", "Early (hotter)", "Late (mild-cold)")) {
  dat <- if (period == "All detections") events else filter(events, season_period == period)
  for (i in seq_len(nrow(pair_mat))) {
    c1 <- pair_mat[i, 1]
    c2 <- pair_mat[i, 2]
    t1 <- dat$time_rad[dat$gazelle_class == c1]
    t2 <- dat$time_rad[dat$gazelle_class == c2]
    res <- calc_pair(
      t1, t2,
      unname(class_labels[c1]),
      unname(class_labels[c2]),
      min_events_group,
      overlap_boot,
      activity_boot
    )
    res$class1 <- c1
    res$class2 <- c2
    res$period <- period
    overall_list[[idx]] <- res
    idx <- idx + 1

    fit1 <- attr(res, "fit1")
    fit2 <- attr(res, "fit2")
    if (!is.null(fit1) && !is.null(fit2)) {
      comp_label <- paste(unname(class_labels[c1]), "vs", unname(class_labels[c2]))
      overall_density[[length(overall_density) + 1]] <- extract_density_df(
        fit1, list(class = unname(class_labels[c1]), comparison = comp_label, period = period)
      )
      overall_density[[length(overall_density) + 1]] <- extract_density_df(
        fit2, list(class = unname(class_labels[c2]), comparison = comp_label, period = period)
      )
    }
  }
}

overall_results <- bind_rows(overall_list) %>%
  mutate(
    comparison = paste(group1, "vs", group2),
    overlap_ci = ifelse(is.na(overlap), NA_character_, fmt_ci(overlap, overlap_lwr, overlap_upr)),
    activity1_fmt = ifelse(is.na(activity1), NA_character_, fmt_num(activity1)),
    activity2_fmt = ifelse(is.na(activity2), NA_character_, fmt_num(activity2)),
    activity_diff_p_fmt = ifelse(is.na(activity_diff_p), NA_character_, fmt_num(activity_diff_p))
  )

overall_density_df <- bind_rows(overall_density)
overall_panel_df <- overall_density_df %>%
  distinct(period, comparison)
overall_shade_df <- merge(overall_panel_df, shade_df, all = TRUE)
overall_annot_df <- overall_results %>%
  transmute(
    period,
    comparison,
    x = 23.8,
    y = Inf,
    label = paste0(
      "Overlap = ", overlap_ci,
      "\n", fmt_p_label(activity_diff_p), " (", fmt_sig_label(activity_diff_p), ")"
    )
  )

## High vs low by season -------------------------------------------------------
hl_list <- list()
hl_density <- list()
idx <- 1

for (period in c("Early (hotter)", "Late (mild-cold)")) {
  period_events <- filter(events, season_period == period)
  for (cls in sort(unique(events$gazelle_class))) {
    class_events <- filter(period_events, gazelle_class == cls)
    for (var in pressure_vars) {
      q1 <- unname(quantile(site_covs[[var]], 0.25, na.rm = TRUE))
      q3 <- unname(quantile(site_covs[[var]], 0.75, na.rm = TRUE))
      low_cams <- site_covs$camera_name[site_covs[[var]] <= q1]
      high_cams <- site_covs$camera_name[site_covs[[var]] >= q3]
      low_times <- class_events$time_rad[class_events$camera_name %in% low_cams]
      high_times <- class_events$time_rad[class_events$camera_name %in% high_cams]
      res <- calc_pair(
        low_times, high_times,
        "Low quartile", "High quartile",
        min_events_group,
        overlap_boot,
        activity_boot
      )
      res$class <- unname(class_labels[cls])
      res$class_code <- cls
      res$covariate <- var
      res$covariate_label <- unname(pretty_var[var])
      res$period <- period
      res$q1 <- q1
      res$q3 <- q3
      res$n_cameras_low <- length(unique(low_cams))
      res$n_cameras_high <- length(unique(high_cams))
      hl_list[[idx]] <- res
      idx <- idx + 1

      fit1 <- attr(res, "fit1")
      fit2 <- attr(res, "fit2")
      if (!is.null(fit1) && !is.null(fit2)) {
        panel_label <- paste(unname(class_labels[cls]), ":", unname(pretty_var[var]), ":", period)
        hl_density[[length(hl_density) + 1]] <- extract_density_df(
          fit1, list(panel = panel_label, group = "Low quartile")
        )
        hl_density[[length(hl_density) + 1]] <- extract_density_df(
          fit2, list(panel = panel_label, group = "High quartile")
        )
      }
    }
  }
}

hl_results <- bind_rows(hl_list) %>%
  mutate(
    overlap_ci = ifelse(is.na(overlap), NA_character_, fmt_ci(overlap, overlap_lwr, overlap_upr)),
    activity_diff_p_fmt = ifelse(is.na(activity_diff_p), NA_character_, fmt_num(activity_diff_p))
  )

hl_density_df <- bind_rows(hl_density)

write_csv(overall_results, paste0(output_prefix, "_overall_class_comparisons.csv"))
write_csv(hl_results, paste0(output_prefix, "_high_low_pressure_comparisons.csv"))

## Plots ----------------------------------------------------------------------
p_overall <- ggplot(overall_density_df, aes(x = hour, y = density, color = class)) +
  geom_rect(
    data = overall_shade_df,
    inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = light_phase),
    alpha = 0.14,
    color = NA
  ) +
  geom_line(linewidth = 0.95) +
  geom_text(
    data = overall_annot_df,
    inherit.aes = FALSE,
    aes(x = x, y = y, label = label),
    hjust = 1.02,
    vjust = 1.15,
    size = 2.8,
    lineheight = 0.95
  ) +
  facet_grid(period ~ comparison) +
  scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0, 24)) +
  scale_fill_manual(values = c("Night" = "#264653", "Dawn" = "#f4a261", "Day" = "#f6e7b4", "Dusk" = "#e76f51")) +
  scale_color_manual(values = c("Female" = "#c1121f", "Male" = "#1d3557", "Fawn" = "#2a9d8f")) +
  labs(
    x = "Apparent solar time",
    y = "Kernel density",
    color = NULL,
    fill = NULL,
    title = "Gazelle diel activity by class in apparent solar time"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(paste0(output_prefix, "_overall_class_activity_curves.png"), p_overall, width = 12, height = 7.8, dpi = 300)
ggsave(paste0(output_prefix, "_overall_class_activity_curves.pdf"), p_overall, width = 12, height = 7.8)

selected_hl <- hl_results %>%
  filter(!is.na(activity_diff_p), activity_diff_p <= 0.05) %>%
  arrange(period, class, overlap, activity_diff_p) %>%
  group_by(period, class) %>%
  slice_head(n = 2) %>%
  ungroup()

selected_panels <- paste(selected_hl$class, ":", selected_hl$covariate_label, ":", selected_hl$period)
hl_plot_df <- hl_density_df %>%
  filter(panel %in% selected_panels)
hl_plot_annot <- selected_hl %>%
  mutate(
    panel = paste(class, ":", covariate_label, ":", period),
    x = 23.8,
    y = Inf,
    label = paste0(
      "Overlap = ", overlap_ci,
      "\n", fmt_p_label(activity_diff_p), " (", fmt_sig_label(activity_diff_p), ")"
    )
  ) %>%
  select(panel, x, y, label)
hl_shade_df <- merge(
  data.frame(panel = unique(hl_plot_df$panel), stringsAsFactors = FALSE),
  shade_df,
  all = TRUE
)

if (nrow(hl_plot_df) > 0) {
  p_hl <- ggplot(hl_plot_df, aes(x = hour, y = density, color = group)) +
    geom_rect(
      data = hl_shade_df,
      inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = light_phase),
      alpha = 0.14,
      color = NA
    ) +
    geom_line(linewidth = 0.9) +
    geom_text(
      data = hl_plot_annot,
      inherit.aes = FALSE,
      aes(x = x, y = y, label = label),
      hjust = 1.02,
      vjust = 1.15,
      size = 2.8,
      lineheight = 0.95
    ) +
    facet_wrap(~ panel, ncol = 2, scales = "free_y") +
    scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0, 24)) +
    scale_fill_manual(values = c("Night" = "#264653", "Dawn" = "#f4a261", "Day" = "#f6e7b4", "Dusk" = "#e76f51"), guide = "none") +
    scale_color_manual(values = c("Low quartile" = "#1b9e77", "High quartile" = "#d95f02")) +
    labs(
      x = "Apparent solar time",
      y = "Kernel density",
      color = NULL,
      title = "Strongest seasonal temporal shifts across pressure gradients"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )

  ggsave(paste0(output_prefix, "_selected_high_low_activity_curves.png"), p_hl, width = 12, height = 8.5, dpi = 300)
  ggsave(paste0(output_prefix, "_selected_high_low_activity_curves.pdf"), p_hl, width = 12, height = 8.5)
}

## Markdown summary -----------------------------------------------------------
month_counts <- table(format(events$ts_utc, "%Y-%m"))
split_label <- format(split_time, "%Y-%m-%d %H:%M:%S UTC", tz = "UTC")

md <- c(
  "# Gazelle Temporal Activity Analysis in Apparent Solar Time",
  "",
  "## Overview",
  "",
  paste0("- Event file: `", event_file, "`"),
  "- Time transformation: UTC timestamps were converted to apparent solar time using event date and event longitude.",
  paste0("- Seasonal split date (median event timestamp): `", split_label, "`"),
  "- Seasonal labels:",
  "  - `Early (hotter)`: detections on or before the split date",
  "  - `Late (mild-cold)`: detections after the split date",
  paste0("- Minimum event threshold per group: `", min_events_group, "`"),
  "",
  "Monthly detection counts:",
  "",
  paste(names(month_counts), month_counts, sep = ": ", collapse = " | "),
  "",
  "## Overall class comparisons",
  "",
  "| Period | Comparison | n1 | n2 | Overlap [95% CI] | Activity1 | Activity2 | Kernel comparison p |",
  "|---|---|---:|---:|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(overall_results))) {
  md <- c(
    md,
    paste0(
      "| `", overall_results$period[i], "` | `", overall_results$comparison[i], "` | `",
      overall_results$n1[i], "` | `", overall_results$n2[i], "` | `",
      ifelse(is.na(overall_results$overlap_ci[i]), "NA", overall_results$overlap_ci[i]), "` | `",
      ifelse(is.na(overall_results$activity1_fmt[i]), "NA", overall_results$activity1_fmt[i]), "` | `",
      ifelse(is.na(overall_results$activity2_fmt[i]), "NA", overall_results$activity2_fmt[i]), "` | `",
      ifelse(is.na(overall_results$activity_diff_p_fmt[i]), "NA", overall_results$activity_diff_p_fmt[i]), "` |"
    )
  )
}

sig_hl <- hl_results %>%
  filter(!is.na(activity_diff_p), activity_diff_p <= 0.05) %>%
  arrange(period, class, overlap, activity_diff_p)

md <- c(
  md,
  "",
  "## Significant high-vs-low quartile temporal shifts",
  "",
  "| Period | Class | Covariate | Low events | High events | Overlap [95% CI] | Kernel comparison p |",
  "|---|---|---|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(sig_hl))) {
  md <- c(
    md,
    paste0(
      "| `", sig_hl$period[i], "` | `", sig_hl$class[i], "` | `", sig_hl$covariate_label[i], "` | `",
      sig_hl$n1[i], "` | `", sig_hl$n2[i], "` | `", sig_hl$overlap_ci[i], "` | `",
      sig_hl$activity_diff_p_fmt[i], "` |"
    )
  )
}

writeLines(md, paste0(output_prefix, "_results.md"))

cat("Saved solar-time seasonal temporal outputs with prefix:", output_prefix, "\n")
