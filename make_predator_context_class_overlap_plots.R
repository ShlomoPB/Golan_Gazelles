if (!requireNamespace("readr", quietly = TRUE)) {
  stop("Package 'readr' is required.")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required.")
}
if (!requireNamespace("activity", quietly = TRUE)) {
  stop("Package 'activity' is required.")
}

library(readr)
library(dplyr)
library(ggplot2)
library(activity)

event_file <- "Gazella gazella_independent_events_30min.csv"
site_cov_file <- "SiteCovs_060324.csv"
output_prefix <- "publication_predator_context_class_overlap_q33_67"
season_output_prefix <- "publication_predator_context_class_overlap_q33_67_by_season"
plot_reps <- 49
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

context_colors <- c(
  "High jackal / Low wolf" = "#9b2226",
  "Low jackal / High wolf" = "#1d3557",
  "High both" = "#6a4c93",
  "Low both" = "#2a9d8f",
  "Intermediate / mixed" = "#8d99ae"
)

class_colors <- c(
  "Female" = "#c1121f",
  "Male" = "#1d3557",
  "Fawn" = "#2a9d8f"
)

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
  )

context_counts <- site_context %>%
  count(predator_context_label, name = "n_cameras") %>%
  mutate(
    facet_label = paste0(predator_context_label, " (n=", n_cameras, ")")
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
  left_join(site_context, by = "camera_name") %>%
  left_join(context_counts, by = "predator_context_label") %>%
  mutate(
    predator_context_label = factor(predator_context_label, levels = context_levels),
    facet_label = factor(
      facet_label,
      levels = context_counts$facet_label[match(context_levels, context_counts$predator_context_label)]
    ),
    class_label = factor(class_label, levels = c("Fawn", "Female", "Male"))
  )

split_time <- stats::median(events$ts_utc)
events <- events %>%
  mutate(
    season_period = ifelse(ts_utc <= split_time, "Early (hotter)", "Late (mild-cold)"),
    season_period = factor(season_period, levels = c("Early (hotter)", "Late (mild-cold)"))
  )

density_list <- list()
for (ctx in context_levels) {
  ctx_dat <- filter(events, predator_context_label == ctx)
  for (cls in c("Fawn", "Female", "Male")) {
    cls_times <- ctx_dat$time_rad[ctx_dat$class_label == cls]
    if (length(cls_times) < 5) {
      next
    }
    fit <- fitact(cls_times, reps = plot_reps, sample = "data", show = FALSE)
    density_list[[length(density_list) + 1]] <- extract_density_df(
      fit,
      list(
        predator_context_label = ctx,
        facet_label = context_counts$facet_label[context_counts$predator_context_label == ctx],
        class = cls
      )
    )
  }
}

density_df <- bind_rows(density_list)

counts_df <- events %>%
  count(predator_context_label, class_label, name = "n_events") %>%
  tidyr::complete(predator_context_label, class_label, fill = list(n_events = 0)) %>%
  arrange(predator_context_label, class_label) %>%
  left_join(
    context_counts,
    by = "predator_context_label"
  ) %>%
  group_by(predator_context_label, facet_label) %>%
  summarise(
    label = paste0(
      "Fawn n=", n_events[class_label == "Fawn"],
      "\nFemale n=", n_events[class_label == "Female"],
      "\nMale n=", n_events[class_label == "Male"]
    ),
    .groups = "drop"
  ) %>% 
  mutate(x = 23.7, y = Inf)

shade_plot_df <- merge(
  data.frame(
    predator_context_label = context_levels,
    facet_label = context_counts$facet_label[match(context_levels, context_counts$predator_context_label)],
    stringsAsFactors = FALSE
  ),
  shade_df,
  by = NULL
)

p <- ggplot(density_df, aes(x = hour, y = density, color = class)) +
  geom_rect(
    data = shade_plot_df,
    inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = light_phase),
    alpha = 0.14,
    color = NA
  ) +
  geom_line(linewidth = 0.95) +
  geom_text(
    data = counts_df,
    inherit.aes = FALSE,
    aes(x = x, y = y, label = label),
    hjust = 1,
    vjust = 1.1,
    size = 2.8,
    lineheight = 0.95
  ) +
  facet_wrap(~ facet_label, ncol = 2, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0, 24)) +
  scale_fill_manual(
    values = c("Night" = "#264653", "Dawn" = "#f4a261", "Day" = "#f6e7b4", "Dusk" = "#e76f51"),
    guide = "none"
  ) +
  scale_color_manual(values = class_colors) +
  labs(
    x = "Apparent solar time",
    y = "Kernel density",
    color = NULL,
    title = "Overlapping gazelle-class activity curves within predator-context categories"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(paste0(output_prefix, ".png"), p, width = 12, height = 9, dpi = 300)
ggsave(paste0(output_prefix, ".pdf"), p, width = 12, height = 9)

season_density_list <- list()
for (period_i in levels(events$season_period)) {
  period_dat <- filter(events, season_period == period_i)
  for (ctx in context_levels) {
    ctx_dat <- filter(period_dat, predator_context_label == ctx)
    for (cls in c("Fawn", "Female", "Male")) {
      cls_times <- ctx_dat$time_rad[ctx_dat$class_label == cls]
      if (length(cls_times) < 5) {
        next
      }
      fit <- fitact(cls_times, reps = plot_reps, sample = "data", show = FALSE)
      season_density_list[[length(season_density_list) + 1]] <- extract_density_df(
        fit,
        list(
          season_period = period_i,
          predator_context_label = ctx,
          facet_label = context_counts$facet_label[context_counts$predator_context_label == ctx],
          class = cls
        )
      )
    }
  }
}

season_density_df <- bind_rows(season_density_list)

season_counts_df <- events %>%
  count(season_period, predator_context_label, class_label, name = "n_events") %>%
  tidyr::complete(season_period, predator_context_label, class_label, fill = list(n_events = 0)) %>%
  arrange(season_period, predator_context_label, class_label) %>%
  left_join(context_counts, by = "predator_context_label") %>%
  group_by(season_period, predator_context_label, facet_label) %>%
  summarise(
    label = paste0(
      "Fawn n=", n_events[class_label == "Fawn"],
      "\nFemale n=", n_events[class_label == "Female"],
      "\nMale n=", n_events[class_label == "Male"]
    ),
    .groups = "drop"
  ) %>%
  mutate(x = 23.7, y = Inf)

season_shade_plot_df <- merge(
  expand.grid(
    season_period = levels(events$season_period),
    predator_context_label = context_levels,
    stringsAsFactors = FALSE
  ) %>%
    left_join(
      context_counts %>% select(predator_context_label, facet_label),
      by = "predator_context_label"
    ),
  shade_df,
  by = NULL
)

p_season <- ggplot(season_density_df, aes(x = hour, y = density, color = class)) +
  geom_rect(
    data = season_shade_plot_df,
    inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = light_phase),
    alpha = 0.14,
    color = NA
  ) +
  geom_line(linewidth = 0.95) +
  geom_text(
    data = season_counts_df,
    inherit.aes = FALSE,
    aes(x = x, y = y, label = label),
    hjust = 1,
    vjust = 1.1,
    size = 2.6,
    lineheight = 0.95
  ) +
  facet_grid(season_period ~ facet_label, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0, 24)) +
  scale_fill_manual(
    values = c("Night" = "#264653", "Dawn" = "#f4a261", "Day" = "#f6e7b4", "Dusk" = "#e76f51"),
    guide = "none"
  ) +
  scale_color_manual(values = class_colors) +
  labs(
    x = "Apparent solar time",
    y = "Kernel density",
    color = NULL,
    title = "Gazelle-class activity curves by predator context and season"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(paste0(season_output_prefix, ".png"), p_season, width = 15, height = 8.8, dpi = 300)
ggsave(paste0(season_output_prefix, ".pdf"), p_season, width = 15, height = 8.8)

cat("Saved predator-context overlapping class plots with prefixes:", output_prefix, "and", season_output_prefix, "\n")
