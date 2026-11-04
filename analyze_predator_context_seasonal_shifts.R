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
if (!requireNamespace("tidyr", quietly = TRUE)) {
  stop("Package 'tidyr' is required.")
}

library(readr)
library(dplyr)
library(ggplot2)
library(overlap)
library(activity)
library(tidyr)

event_file <- "Gazella gazella_independent_events_30min.csv"
site_cov_file <- "SiteCovs_060324.csv"
output_prefix <- "publication_predator_context_q33_67_seasonal_shifts"
overlap_boot <- 99
activity_boot <- 99
plot_reps <- 49
min_events_group <- 20
lower_q <- 0.33
upper_q <- 0.67

class_labels <- c(
  "juvenile" = "Fawn",
  "male" = "Male",
  "female" = "Female"
)

context_levels <- c(
  "High jackal / Low wolf",
  "Low jackal / High wolf",
  "High both",
  "Low both",
  "Intermediate / mixed"
)

season_colors <- c(
  "Early (hotter)" = "#d62828",
  "Late (mild-cold)" = "#1d3557"
)

fmt_num <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
fmt_ci <- function(mean, lwr, upr, digits = 3) {
  paste0(fmt_num(mean, digits), " [", fmt_num(lwr, digits), ", ", fmt_num(upr, digits), "]")
}
fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", fmt_num(x, 3)))
}
fmt_sig <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(x < 0.001, "***", ifelse(x < 0.01, "**", ifelse(x < 0.05, "*", "ns")))
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

site_covs <- read_csv(site_cov_file, show_col_types = FALSE)
q_j1 <- quantile(site_covs$null.jackal.abundance, lower_q, na.rm = TRUE)
q_j3 <- quantile(site_covs$null.jackal.abundance, upper_q, na.rm = TRUE)
q_w1 <- quantile(site_covs$Null.Wolf.Abundance, lower_q, na.rm = TRUE)
q_w3 <- quantile(site_covs$Null.Wolf.Abundance, upper_q, na.rm = TRUE)

site_context <- site_covs %>%
  transmute(
    camera_name = as.integer(Camera.ID),
    predator_context_label = case_when(
      null.jackal.abundance >= q_j3 & Null.Wolf.Abundance <= q_w1 ~ "High jackal / Low wolf",
      null.jackal.abundance <= q_j1 & Null.Wolf.Abundance >= q_w3 ~ "Low jackal / High wolf",
      null.jackal.abundance >= q_j3 & Null.Wolf.Abundance >= q_w3 ~ "High both",
      null.jackal.abundance <= q_j1 & Null.Wolf.Abundance <= q_w1 ~ "Low both",
      TRUE ~ "Intermediate / mixed"
    )
  ) %>%
  mutate(predator_context_label = factor(predator_context_label, levels = context_levels))

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
    season_period = ifelse(ts_utc <= split_time, "Early (hotter)", "Late (mild-cold)"),
    season_period = factor(season_period, levels = c("Early (hotter)", "Late (mild-cold)")),
    class_label = factor(class_label, levels = c("Fawn", "Female", "Male"))
  )

comparison_list <- list()
density_list <- list()
idx <- 1

for (ctx in context_levels) {
  for (cls in c("juvenile", "female", "male")) {
    dat <- events %>%
      filter(predator_context_label == ctx, gazelle_class == cls)
    times_early <- dat$time_rad[dat$season_period == "Early (hotter)"]
    times_late <- dat$time_rad[dat$season_period == "Late (mild-cold)"]
    res <- calc_pair(
      times_early,
      times_late,
      "Early (hotter)",
      "Late (mild-cold)",
      min_events_group,
      overlap_boot,
      activity_boot
    )
    res$predator_context_label <- ctx
    res$class <- unname(class_labels[cls])
    res$class_code <- cls
    comparison_list[[idx]] <- res
    idx <- idx + 1

    fit1 <- attr(res, "fit1")
    fit2 <- attr(res, "fit2")
    if (!is.null(fit1) && !is.null(fit2)) {
      panel <- paste(unname(class_labels[cls]), ":", ctx)
      density_list[[length(density_list) + 1]] <- extract_density_df(
        fit1,
        list(panel = panel, class = unname(class_labels[cls]), predator_context_label = ctx, season_period = "Early (hotter)")
      )
      density_list[[length(density_list) + 1]] <- extract_density_df(
        fit2,
        list(panel = panel, class = unname(class_labels[cls]), predator_context_label = ctx, season_period = "Late (mild-cold)")
      )
    }
  }
}

comparison_df <- bind_rows(comparison_list) %>%
  mutate(
    overlap_ci = ifelse(is.na(overlap), NA_character_, fmt_ci(overlap, overlap_lwr, overlap_upr)),
    activity_diff_p_fmt = fmt_p(activity_diff_p),
    sig = fmt_sig(activity_diff_p)
  )

density_df <- bind_rows(density_list)

counts_df <- events %>%
  count(predator_context_label, class_label, season_period, name = "n_events") %>%
  tidyr::complete(predator_context_label, class_label, season_period, fill = list(n_events = 0)) %>%
  group_by(predator_context_label, class_label) %>%
  summarise(
    early_n = n_events[season_period == "Early (hotter)"],
    late_n = n_events[season_period == "Late (mild-cold)"],
    .groups = "drop"
  )

write_csv(comparison_df, paste0(output_prefix, "_comparisons.csv"))
write_csv(counts_df, paste0(output_prefix, "_event_counts.csv"))

if (nrow(density_df) > 0) {
  panel_levels <- comparison_df %>%
    mutate(panel = paste(class, ":", predator_context_label)) %>%
    arrange(match(class, c("Fawn", "Female", "Male")), match(predator_context_label, context_levels)) %>%
    pull(panel)

  density_df$panel <- factor(density_df$panel, levels = unique(panel_levels))
  shade_plot_df <- merge(
    data.frame(panel = unique(density_df$panel), stringsAsFactors = FALSE),
    shade_df,
    by = NULL
  )

  annot_df <- comparison_df %>%
    left_join(counts_df, by = c("predator_context_label", "class" = "class_label")) %>%
    mutate(
      panel = paste(class, ":", predator_context_label),
      x = 23.7,
      y = Inf,
      label = paste0(
        "Early n=", early_n,
        "\nLate n=", late_n,
        "\nOverlap=", ifelse(is.na(overlap_ci), "NA", overlap_ci),
        "\np=", activity_diff_p_fmt,
        " ", sig
      )
    ) %>%
    select(panel, x, y, label)

  p <- ggplot(density_df, aes(x = hour, y = density, color = season_period)) +
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
      size = 2.6,
      lineheight = 0.95
    ) +
    facet_wrap(~ panel, ncol = 3, scales = "free_y") +
    scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0, 24)) +
    scale_fill_manual(
      values = c("Night" = "#264653", "Dawn" = "#f4a261", "Day" = "#f6e7b4", "Dusk" = "#e76f51"),
      guide = "none"
    ) +
    scale_color_manual(values = season_colors) +
    labs(
      x = "Apparent solar time",
      y = "Kernel density",
      color = NULL,
      title = "Seasonal shifts in gazelle activity within predator-context categories"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )

  ggsave(paste0(output_prefix, ".png"), p, width = 14, height = 11, dpi = 300)
  ggsave(paste0(output_prefix, ".pdf"), p, width = 14, height = 11)
}

md <- c(
  "# Seasonal shifts within predator-context categories",
  "",
  paste0(
    "Predator-context categories used the 33rd and 67th percentiles of jackal and wolf abundance. ",
    "Season was split at the median gazelle detection timestamp (`", as.character(split_time), "` UTC)."
  ),
  "",
  "## Event counts",
  "",
  "| Predator context | Class | Early events | Late events |",
  "|---|---|---:|---:|"
)

for (i in seq_len(nrow(counts_df))) {
  md <- c(
    md,
    paste0(
      "| `", counts_df$predator_context_label[i], "` | `", counts_df$class_label[i],
      "` | `", counts_df$early_n[i], "` | `", counts_df$late_n[i], "` |"
    )
  )
}

md <- c(
  md,
  "",
  "## Early vs late activity comparisons",
  "",
  "| Predator context | Class | Early n | Late n | Overlap [95% CI] | Kernel p |",
  "|---|---|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(comparison_df))) {
  md <- c(
    md,
    paste0(
      "| `", comparison_df$predator_context_label[i], "` | `", comparison_df$class[i], "` | `",
      comparison_df$n1[i], "` | `", comparison_df$n2[i], "` | `",
      ifelse(is.na(comparison_df$overlap_ci[i]), "NA", comparison_df$overlap_ci[i]),
      "` | `", comparison_df$activity_diff_p_fmt[i], "` |"
    )
  )
}

sig_df <- comparison_df %>%
  filter(!is.na(activity_diff_p), activity_diff_p <= 0.05) %>%
  arrange(match(class, c("Fawn", "Female", "Male")), activity_diff_p)

md <- c(md, "", "## Supported seasonal shifts", "")
if (nrow(sig_df) == 0) {
  md <- c(md, "No within-context early versus late activity differences met the testing threshold.")
} else {
  for (i in seq_len(nrow(sig_df))) {
    md <- c(
      md,
      paste0(
        "- `", sig_df$class[i], "` in `", sig_df$predator_context_label[i],
        "` shifted between early and late periods (overlap `", sig_df$overlap_ci[i],
        "`, `p = ", sig_df$activity_diff_p_fmt[i], "`)."
      )
    )
  }
}

writeLines(md, paste0(output_prefix, "_results.md"))

cat("Saved seasonal predator-context shift outputs with prefix:", output_prefix, "\n")
