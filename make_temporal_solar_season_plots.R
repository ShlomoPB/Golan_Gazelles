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

library(activity)
library(ggplot2)
library(dplyr)
library(readr)

event_file <- "Gazella gazella_independent_events_30min.csv"
site_cov_file <- "SiteCovs_060324.csv"
overall_file <- "gazelle_temporal_activity_solar_season_overall_class_comparisons.csv"
hl_file <- "gazelle_temporal_activity_solar_season_high_low_pressure_comparisons.csv"
output_prefix <- "gazelle_temporal_activity_solar_season"
plot_reps <- 49

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

extract_density_df <- function(fit, extra) {
  pdf_df <- as.data.frame(fit@pdf)
  data.frame(
    hour = to_hours(pdf_df$x),
    density = pdf_df$y,
    extra,
    stringsAsFactors = FALSE
  )
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
    solar_hour = to_solar_hour(ts_utc, camelot_gps_longitude),
    time_rad = to_radians(solar_hour)
  )

site_covs <- read_csv(site_cov_file, show_col_types = FALSE) %>%
  rename(camera_name = Camera.ID)

split_time <- stats::median(events$ts_utc)
events <- events %>%
  mutate(season_period = ifelse(ts_utc <= split_time, "Early (hotter)", "Late (mild-cold)"))

overall_results <- read_csv(overall_file, show_col_types = FALSE)
hl_results <- read_csv(hl_file, show_col_types = FALSE)

pair_map <- data.frame(
  comparison = c("Fawn vs Male", "Female vs Fawn", "Female vs Male"),
  class1 = c("juvenile", "female", "female"),
  class2 = c("male", "juvenile", "male"),
  label1 = c("Fawn", "Female", "Female"),
  label2 = c("Male", "Fawn", "Male"),
  stringsAsFactors = FALSE
)

overall_density <- list()
for (i in seq_len(nrow(pair_map))) {
  for (period_i in c("All detections", "Early (hotter)", "Late (mild-cold)")) {
    dat <- if (period_i == "All detections") events else filter(events, season_period == period_i)
    fit1 <- fitact(dat$time_rad[dat$gazelle_class == pair_map$class1[i]], reps = plot_reps, sample = "data", show = FALSE)
    fit2 <- fitact(dat$time_rad[dat$gazelle_class == pair_map$class2[i]], reps = plot_reps, sample = "data", show = FALSE)
    overall_density[[length(overall_density) + 1]] <- extract_density_df(
      fit1,
      list(comparison = pair_map$comparison[i], period = period_i, class = pair_map$label1[i])
    )
    overall_density[[length(overall_density) + 1]] <- extract_density_df(
      fit2,
      list(comparison = pair_map$comparison[i], period = period_i, class = pair_map$label2[i])
    )
  }
}
overall_density_df <- bind_rows(overall_density)
overall_shade_df <- merge(
  unique(overall_density_df[c("period", "comparison")]),
  shade_df,
  by = NULL
)
overall_annot_df <- overall_results %>%
  transmute(
    period,
    comparison,
    x = 23.7,
    y = Inf,
    label = paste0("Overlap = ", overlap_ci, "\n", fmt_p_label(activity_diff_p), " (", fmt_sig_label(activity_diff_p), ")")
  )

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
    hjust = 1,
    vjust = 1.1,
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

selected_density <- list()
for (i in seq_len(nrow(selected_hl))) {
  cls_code <- c("Fawn" = "juvenile", "Male" = "male", "Female" = "female")[selected_hl$class[i]]
  var <- names(pretty_var)[pretty_var == selected_hl$covariate_label[i]][1]
  q1 <- quantile(site_covs[[var]], 0.25, na.rm = TRUE)
  q3 <- quantile(site_covs[[var]], 0.75, na.rm = TRUE)
  low_cams <- site_covs$camera_name[site_covs[[var]] <= q1]
  high_cams <- site_covs$camera_name[site_covs[[var]] >= q3]
  class_events <- events %>%
    filter(gazelle_class == cls_code, season_period == selected_hl$period[i])
  fit_low <- fitact(class_events$time_rad[class_events$camera_name %in% low_cams], reps = plot_reps, sample = "data", show = FALSE)
  fit_high <- fitact(class_events$time_rad[class_events$camera_name %in% high_cams], reps = plot_reps, sample = "data", show = FALSE)
  panel_label <- paste(selected_hl$class[i], ":", selected_hl$covariate_label[i], ":", selected_hl$period[i])
  selected_density[[length(selected_density) + 1]] <- extract_density_df(
    fit_low,
    list(panel = panel_label, group = "Low quartile")
  )
  selected_density[[length(selected_density) + 1]] <- extract_density_df(
    fit_high,
    list(panel = panel_label, group = "High quartile")
  )
}

selected_density_df <- bind_rows(selected_density)
hl_shade_df <- merge(
  data.frame(panel = unique(selected_density_df$panel), stringsAsFactors = FALSE),
  shade_df,
  by = NULL
)
hl_annot_df <- selected_hl %>%
  mutate(
    panel = paste(class, ":", covariate_label, ":", period),
    x = 23.7,
    y = Inf,
    label = paste0("Overlap = ", overlap_ci, "\n", fmt_p_label(activity_diff_p), " (", fmt_sig_label(activity_diff_p), ")")
  ) %>%
  select(panel, x, y, label)

p_hl <- ggplot(selected_density_df, aes(x = hour, y = density, color = group)) +
  geom_rect(
    data = hl_shade_df,
    inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = light_phase),
    alpha = 0.14,
    color = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_text(
    data = hl_annot_df,
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

cat("Saved corrected solar-season temporal plots\n")
