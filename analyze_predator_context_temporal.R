if (!requireNamespace("readr", quietly = TRUE)) {
  stop("Package 'readr' is required.")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required.")
}
if (!requireNamespace("overlap", quietly = TRUE)) {
  stop("Package 'overlap' is required.")
}
if (!requireNamespace("activity", quietly = TRUE)) {
  stop("Package 'activity' is required.")
}
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required.")
}

library(readr)
library(dplyr)
library(ggplot2)
library(overlap)
library(activity)
library(readxl)

event_file <- "Gazella gazella_independent_events_30min.csv"
site_cov_file <- "SiteCovs_060324.csv"
coords_file <- "Camera Coordinates.xlsx"
output_prefix <- "publication_predator_context_q33_67"
overlap_boot <- 99
activity_boot <- 99
plot_reps <- 29
min_events_group <- 20
lower_q <- 0.33
upper_q <- 0.67

class_labels <- c(
  "juvenile" = "Fawn",
  "male" = "Male",
  "female" = "Female"
)

pretty_context <- c(
  "high_jackal_low_wolf" = "High jackal / Low wolf",
  "low_jackal_high_wolf" = "Low jackal / High wolf",
  "high_both" = "High both",
  "low_both" = "Low both",
  "mixed" = "Intermediate / mixed"
)

context_colors <- c(
  "High jackal / Low wolf" = "#9b2226",
  "Low jackal / High wolf" = "#1d3557",
  "High both" = "#6a4c93",
  "Low both" = "#2a9d8f",
  "Intermediate / mixed" = "#8d99ae"
)

contrast_specs <- data.frame(
  contrast_id = c("HJLW_vs_LJHW", "HB_vs_LB", "Extreme_vs_Mixed"),
  group1_code = c("high_jackal_low_wolf", "high_both", "extreme"),
  group2_code = c("low_jackal_high_wolf", "low_both", "mixed"),
  group1_label = c("High jackal / Low wolf", "High both", "Extreme predator contexts"),
  group2_label = c("Low jackal / High wolf", "Low both", "Intermediate / mixed"),
  stringsAsFactors = FALSE
)

fmt_num <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
fmt_ci <- function(mean, lwr, upr, digits = 3) {
  paste0(fmt_num(mean, digits), " [", fmt_num(lwr, digits), ", ", fmt_num(upr, digits), "]")
}
fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", fmt_num(x, 3)))
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

site_covs <- read_csv(site_cov_file, show_col_types = FALSE)
q_j1 <- quantile(site_covs$null.jackal.abundance, lower_q, na.rm = TRUE)
q_j3 <- quantile(site_covs$null.jackal.abundance, upper_q, na.rm = TRUE)
q_w1 <- quantile(site_covs$Null.Wolf.Abundance, lower_q, na.rm = TRUE)
q_w3 <- quantile(site_covs$Null.Wolf.Abundance, upper_q, na.rm = TRUE)

site_context <- site_covs %>%
  transmute(
    camera_name = as.integer(Camera.ID),
    Polygon.Name,
    jackal = null.jackal.abundance,
    wolf = Null.Wolf.Abundance,
    predator_context = case_when(
      jackal >= q_j3 & wolf <= q_w1 ~ "high_jackal_low_wolf",
      jackal <= q_j1 & wolf >= q_w3 ~ "low_jackal_high_wolf",
      jackal >= q_j3 & wolf >= q_w3 ~ "high_both",
      jackal <= q_j1 & wolf <= q_w1 ~ "low_both",
      TRUE ~ "mixed"
    ),
    predator_context_label = unname(pretty_context[predator_context]),
    predator_context_group = ifelse(predator_context == "mixed", "mixed", "extreme")
  )

events <- read_csv(event_file, show_col_types = FALSE) %>%
  mutate(
    camera_name = as.integer(camera_name),
    gazelle_class = as.character(gazelle_class),
    class_label = unname(class_labels[gazelle_class]),
    ts_utc = as.POSIXct(event_time, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    solar_hour = to_solar_hour(ts_utc, camelot_gps_longitude),
    time_rad = to_radians(solar_hour)
  ) %>%
  left_join(site_context, by = "camera_name")

split_time <- stats::median(events$ts_utc)
events <- events %>%
  mutate(
    season_period = ifelse(ts_utc <= split_time, "Early (hotter)", "Late (mild-cold)")
  )

coords_raw <- read_excel(coords_file)
names(coords_raw) <- c("site_label", "lat", "lon")
coords_raw$camera_name <- as.integer(gsub("[^0-9]", "", as.character(coords_raw$site_label)))
camera_map <- site_context %>%
  left_join(coords_raw[, c("camera_name", "lat", "lon")], by = "camera_name")

category_summary <- site_context %>%
  count(predator_context_label, name = "n_cameras") %>%
  arrange(match(predator_context_label, pretty_context))

event_summary <- events %>%
  count(predator_context_label, class_label, season_period, name = "n_events") %>%
  arrange(match(predator_context_label, pretty_context), class_label, season_period)

comparison_results <- list()
density_list <- list()
idx <- 1

for (period_i in c("All detections", "Early (hotter)", "Late (mild-cold)")) {
  period_events <- if (period_i == "All detections") events else filter(events, season_period == period_i)
  for (cls in sort(unique(period_events$gazelle_class))) {
    class_events <- filter(period_events, gazelle_class == cls)
    for (k in seq_len(nrow(contrast_specs))) {
      spec <- contrast_specs[k, ]
      group1_var <- if (spec$group1_code == "extreme") "predator_context_group" else "predator_context"
      group2_var <- if (spec$group2_code == "extreme") "predator_context_group" else "predator_context"
      times1 <- if (group1_var == "predator_context_group") {
        class_events$time_rad[class_events$predator_context_group == spec$group1_code]
      } else {
        class_events$time_rad[class_events$predator_context == spec$group1_code]
      }
      times2 <- if (group2_var == "predator_context_group") {
        class_events$time_rad[class_events$predator_context_group == spec$group2_code]
      } else {
        class_events$time_rad[class_events$predator_context == spec$group2_code]
      }
      res <- calc_pair(
        times1, times2,
        spec$group1_label,
        spec$group2_label,
        min_events_group,
        overlap_boot,
        activity_boot
      )
      res$class <- unname(class_labels[cls])
      res$class_code <- cls
      res$period <- period_i
      res$contrast_id <- spec$contrast_id
      comparison_results[[idx]] <- res
      idx <- idx + 1

      fit1 <- attr(res, "fit1")
      fit2 <- attr(res, "fit2")
      if (!is.null(fit1) && !is.null(fit2)) {
        panel <- paste(unname(class_labels[cls]), ":", spec$group1_label, "vs", spec$group2_label, ":", period_i)
        density_list[[length(density_list) + 1]] <- extract_density_df(
          fit1, list(panel = panel, group = spec$group1_label, class = unname(class_labels[cls]), period = period_i)
        )
        density_list[[length(density_list) + 1]] <- extract_density_df(
          fit2, list(panel = panel, group = spec$group2_label, class = unname(class_labels[cls]), period = period_i)
        )
      }
    }
  }
}

comparison_df <- bind_rows(comparison_results) %>%
  mutate(
    overlap_ci = ifelse(is.na(overlap), NA_character_, fmt_ci(overlap, overlap_lwr, overlap_upr)),
    activity_diff_p_fmt = fmt_p(activity_diff_p)
  )

density_df <- bind_rows(density_list)

write_csv(site_context, paste0(output_prefix, "_camera_categories.csv"))
write_csv(category_summary, paste0(output_prefix, "_camera_category_counts.csv"))
write_csv(event_summary, paste0(output_prefix, "_event_summary.csv"))
write_csv(comparison_df, paste0(output_prefix, "_temporal_comparisons.csv"))

map_plot <- ggplot(camera_map, aes(x = lon, y = lat, color = predator_context_label)) +
  geom_point(size = 3.1, alpha = 0.95) +
  scale_color_manual(values = context_colors, drop = FALSE) +
  coord_equal() +
  labs(
    x = "Longitude",
    y = "Latitude",
    color = NULL,
    title = "Camera predator-context categories"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(paste0(output_prefix, "_camera_map.png"), map_plot, width = 8.2, height = 5.6, dpi = 300)
ggsave(paste0(output_prefix, "_camera_map.pdf"), map_plot, width = 8.2, height = 5.6)

selected_comp <- comparison_df %>%
  filter(!is.na(activity_diff_p)) %>%
  arrange(period, class, activity_diff_p, overlap) %>%
  group_by(period, class) %>%
  slice_head(n = 1) %>%
  ungroup()

selected_panels <- paste(selected_comp$class, ":", selected_comp$group1, "vs", selected_comp$group2, ":", selected_comp$period)
selected_density_df <- density_df %>%
  filter(panel %in% selected_panels)

if (nrow(selected_density_df) > 0) {
  shade_plot_df <- merge(data.frame(panel = unique(selected_density_df$panel), stringsAsFactors = FALSE), shade_df, by = NULL)
  annot_df <- selected_comp %>%
    mutate(
      panel = paste(class, ":", group1, "vs", group2, ":", period),
      x = 23.7,
      y = Inf,
      label = paste0("Overlap = ", overlap_ci, "\np = ", activity_diff_p_fmt)
    ) %>%
    select(panel, x, y, label)

  temporal_plot <- ggplot(selected_density_df, aes(x = hour, y = density, color = group)) +
    geom_rect(
      data = shade_plot_df,
      inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = light_phase),
      alpha = 0.14,
      color = NA
    ) +
    geom_line(linewidth = 0.95) +
    geom_text(
      data = annot_df,
      inherit.aes = FALSE,
      aes(x = x, y = y, label = label),
      hjust = 1,
      vjust = 1.1,
      size = 2.8,
      lineheight = 0.95
    ) +
    facet_wrap(~ panel, ncol = 2, scales = "free_y") +
    scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0, 24)) +
    scale_fill_manual(values = c("Night" = "#264653", "Dawn" = "#f4a261", "Day" = "#f6e7b4", "Dusk" = "#e76f51"), guide = "none") +
    scale_color_manual(values = c(
      "High jackal / Low wolf" = "#9b2226",
      "Low jackal / High wolf" = "#1d3557",
      "High both" = "#6a4c93",
      "Low both" = "#2a9d8f",
      "Extreme predator contexts" = "#7f5539",
      "Intermediate / mixed" = "#8d99ae"
    )) +
    labs(
      x = "Apparent solar time",
      y = "Kernel density",
      color = NULL,
      title = "Selected class-specific temporal contrasts across predator contexts"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )

  ggsave(paste0(output_prefix, "_selected_temporal_contrasts.png"), temporal_plot, width = 12, height = 8.5, dpi = 300)
  ggsave(paste0(output_prefix, "_selected_temporal_contrasts.pdf"), temporal_plot, width = 12, height = 8.5)
}

md <- c(
  "# Predator-context temporal analysis",
  "",
  "## Camera classification",
  "",
  paste0(
    "Cameras were classified using quartiles of `null.jackal.abundance` and `Null.Wolf.Abundance`: jackal Q1 = `", fmt_num(q_j1),
    "`, jackal Q3 = `", fmt_num(q_j3),
    "`, wolf Q1 = `", fmt_num(q_w1),
    "`, and wolf Q3 = `", fmt_num(q_w3), "`."
  ),
  "",
  "| Category | Cameras |",
  "|---|---:|"
)

for (i in seq_len(nrow(category_summary))) {
  md <- c(md, paste0("| `", category_summary$predator_context_label[i], "` | `", category_summary$n_cameras[i], "` |"))
}

md <- c(
  md,
  "",
  "## Event totals by class and predator context",
  "",
  "| Category | Class | Period | Events |",
  "|---|---|---|---:|"
)
for (i in seq_len(nrow(event_summary))) {
  md <- c(
    md,
    paste0(
      "| `", event_summary$predator_context_label[i], "` | `", event_summary$class_label[i], "` | `",
      event_summary$season_period[i], "` | `", event_summary$n_events[i], "` |"
    )
  )
}

md <- c(
  md,
  "",
  "## Temporal comparisons",
  "",
  "Only comparisons with at least the minimum event threshold in both groups were tested formally.",
  "",
  "| Period | Class | Contrast | n1 | n2 | Overlap [95% CI] | Kernel p |",
  "|---|---|---|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(comparison_df))) {
  md <- c(
    md,
    paste0(
      "| `", comparison_df$period[i], "` | `", comparison_df$class[i], "` | `",
      comparison_df$group1[i], " vs ", comparison_df$group2[i], "` | `",
      comparison_df$n1[i], "` | `", comparison_df$n2[i], "` | `",
      ifelse(is.na(comparison_df$overlap_ci[i]), "NA", comparison_df$overlap_ci[i]), "` | `",
      comparison_df$activity_diff_p_fmt[i], "` |"
    )
  )
}

tested_df <- comparison_df %>% filter(!is.na(activity_diff_p))
if (nrow(tested_df) > 0) {
  sig_df <- tested_df %>% filter(activity_diff_p <= 0.05) %>% arrange(activity_diff_p, overlap)
  md <- c(md, "", "## Supported contrasts", "")
  if (nrow(sig_df) == 0) {
    md <- c(md, "No predator-context temporal contrasts were clearly supported at `p < 0.05` after applying the minimum event threshold.")
  } else {
    for (i in seq_len(nrow(sig_df))) {
      md <- c(
        md,
        paste0(
          "- In `", sig_df$period[i], "`, `", sig_df$class[i], "` detections differed between `",
          sig_df$group1[i], "` and `", sig_df$group2[i], "` cameras (overlap `",
          sig_df$overlap_ci[i], "`, `p = ", sig_df$activity_diff_p_fmt[i], "`)."
        )
      )
    }
  }
}

writeLines(md, paste0(output_prefix, "_results.md"))

cat("Saved predator-context temporal outputs with prefix:", output_prefix, "\n")
