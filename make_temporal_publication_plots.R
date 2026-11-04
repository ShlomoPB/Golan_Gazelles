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

event_file <- "Gazella gazella_independent_events_30min.csv"
site_cov_file <- "SiteCovs_060324.csv"
timezone_local <- "Asia/Jerusalem"

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
  "Built.and.Roads" = "Built and roads",
  "Nature.Reserves" = "Nature reserves",
  "Open.Spaces" = "Open spaces",
  "Cropland" = "Cropland",
  "Orchards" = "Orchards"
)

extract_hours <- function(x, tz_local) {
  ts_utc <- as.POSIXct(x, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  ts_local <- as.POSIXct(format(ts_utc, tz = tz_local, usetz = TRUE), tz = tz_local)
  as.numeric(format(ts_local, "%H")) +
    as.numeric(format(ts_local, "%M")) / 60 +
    as.numeric(format(ts_local, "%S")) / 3600
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
  xmin = c(0, 5.5, 7.0, 18.5),
  xmax = c(5.5, 7.0, 18.5, 24),
  period = c("Night", "Crepuscular", "Day", "Crepuscular/Night")
)

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

## Overall class curves -------------------------------------------------------
overall_density <- list()
for (cls in sort(unique(events$gazelle_class))) {
  fit <- fitact(events$time_rad[events$gazelle_class == cls], reps = 199, sample = "data", show = FALSE)
  overall_density[[length(overall_density) + 1]] <- extract_density_df(
    fit, list(class = unname(class_labels[cls]))
  )
}
overall_density_df <- bind_rows(overall_density)

p_overall <- ggplot(overall_density_df, aes(x = hour, y = density, color = class)) +
  geom_rect(
    data = shade_df,
    inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = period),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0, 24)) +
  scale_fill_manual(
    values = c(
      "Night" = "#203a43",
      "Crepuscular" = "#f4a261",
      "Day" = "#f6e7b4",
      "Crepuscular/Night" = "#34515e"
    ),
    guide = "none"
  ) +
  scale_color_manual(values = c("Female" = "#c1121f", "Male" = "#1d3557", "Fawn" = "#2a9d8f")) +
  labs(
    x = paste0("Local time (", timezone_local, ")"),
    y = "Kernel density",
    color = NULL,
    title = "Diel activity patterns of gazelle classes"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave("publication_temporal_overall_activity.png", p_overall, width = 8.5, height = 4.8, dpi = 300)
ggsave("publication_temporal_overall_activity.pdf", p_overall, width = 8.5, height = 4.8)

## Selected high-low curves ---------------------------------------------------
hl <- read_csv("gazelle_temporal_activity_high_low_pressure_comparisons.csv", show_col_types = FALSE)
selected <- hl %>%
  filter(!is.na(activity_diff_p), activity_diff_p <= 0.05) %>%
  arrange(class, overlap, activity_diff_p) %>%
  group_by(class) %>%
  slice_head(n = 2) %>%
  ungroup()

selected_density <- list()

for (i in seq_len(nrow(selected))) {
  cls_code <- c("Fawn" = "juvenile", "Male" = "male", "Female" = "female")[selected$class[i]]
  var <- names(pretty_var)[pretty_var == selected$covariate_label[i]][1]
  q1 <- quantile(site_covs[[var]], 0.25, na.rm = TRUE)
  q3 <- quantile(site_covs[[var]], 0.75, na.rm = TRUE)
  low_cams <- site_covs$camera_name[site_covs[[var]] <= q1]
  high_cams <- site_covs$camera_name[site_covs[[var]] >= q3]
  class_events <- events %>% filter(gazelle_class == cls_code)
  low_times <- class_events$time_rad[class_events$camera_name %in% low_cams]
  high_times <- class_events$time_rad[class_events$camera_name %in% high_cams]
  fit_low <- fitact(low_times, reps = 199, sample = "data", show = FALSE)
  fit_high <- fitact(high_times, reps = 199, sample = "data", show = FALSE)
  facet_label <- paste(selected$class[i], ":", selected$covariate_label[i])
  selected_density[[length(selected_density) + 1]] <- extract_density_df(
    fit_low, list(panel = facet_label, group = "Low quartile")
  )
  selected_density[[length(selected_density) + 1]] <- extract_density_df(
    fit_high, list(panel = facet_label, group = "High quartile")
  )
}

selected_density_df <- bind_rows(selected_density)

p_selected <- ggplot(selected_density_df, aes(x = hour, y = density, color = group)) +
  geom_rect(
    data = shade_df,
    inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = period),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 0.95) +
  facet_wrap(~ panel, ncol = 2, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0, 24)) +
  scale_fill_manual(
    values = c(
      "Night" = "#203a43",
      "Crepuscular" = "#f4a261",
      "Day" = "#f6e7b4",
      "Crepuscular/Night" = "#34515e"
    ),
    guide = "none"
  ) +
  scale_color_manual(values = c("Low quartile" = "#1b9e77", "High quartile" = "#d95f02")) +
  labs(
    x = paste0("Local time (", timezone_local, ")"),
    y = "Kernel density",
    color = NULL,
    title = "Strongest temporal shifts across pressure gradients"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave("publication_temporal_selected_high_low_activity.png", p_selected, width = 11, height = 7.2, dpi = 300)
ggsave("publication_temporal_selected_high_low_activity.pdf", p_selected, width = 11, height = 7.2)

cat("Saved publication temporal plots\n")
