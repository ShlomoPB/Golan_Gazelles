#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

output_prefix <- "publication_fawn_female_count_ratio_wolf_jackal"
juvenile_file <- "gazelle_juvenile_abundance_60x7.csv"
female_file <- "gazelle_female_abundance_60x7.csv"
site_cov_file <- "SiteCovs_060324.csv"
obs_cov_file <- "ObsCovs_050524.csv"
event_file <- "Gazella gazella_independent_events_30min.csv"
lower_q <- 0.33
upper_q <- 0.67

read_count_long <- function(path, value_name) {
  read.csv(path, na.strings = c("NA", ""), check.names = FALSE) %>%
    as_tibble() %>%
    mutate(camera_id = as.integer(camera_name)) %>%
    pivot_longer(
      cols = all_of(paste0("X", 1:7)),
      names_to = "round_label",
      values_to = value_name
    ) %>%
    mutate(round = as.integer(sub("^X", "", round_label))) %>%
    select(camera_id, round, all_of(value_name))
}

predator_level <- function(x, lower, upper) {
  factor(
    case_when(
      x <= lower ~ "Low",
      x >= upper ~ "High",
      TRUE ~ "Middle"
    ),
    levels = c("Low", "Middle", "High")
  )
}

fawn_counts <- read_count_long(juvenile_file, "fawn")
female_counts <- read_count_long(female_file, "female")

counts <- full_join(fawn_counts, female_counts, by = c("camera_id", "round"))
if (any(xor(is.na(counts$fawn), is.na(counts$female)))) {
  stop("Fawn and female count histories have different missing-data patterns.")
}

site_covs <- read_csv(site_cov_file, show_col_types = FALSE)
wolf_cutoffs <- quantile(
  site_covs$Null.Wolf.Abundance,
  probs = c(lower_q, upper_q),
  na.rm = TRUE,
  names = FALSE
)
jackal_cutoffs <- quantile(
  site_covs$null.jackal.abundance,
  probs = c(lower_q, upper_q),
  na.rm = TRUE,
  names = FALSE
)

site_groups <- site_covs %>%
  transmute(
    camera_id = as.integer(Camera.ID),
    wolf_abundance = Null.Wolf.Abundance,
    jackal_abundance = null.jackal.abundance,
    wolf_level = predator_level(Null.Wolf.Abundance, wolf_cutoffs[1], wolf_cutoffs[2]),
    jackal_level = predator_level(null.jackal.abundance, jackal_cutoffs[1], jackal_cutoffs[2])
  )

# Match the season assignment used by the seasonal spAbundance analysis.
events <- read_csv(event_file, show_col_types = FALSE) %>%
  mutate(ts_utc = as.POSIXct(event_time, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
split_time <- median(events$ts_utc, na.rm = TRUE)
split_date <- as.Date(split_time, tz = "UTC")

obs_covs <- read.csv(obs_cov_file, na.strings = c("NA", ""), check.names = FALSE) %>%
  as_tibble() %>%
  mutate(
    camera_id = as.integer(`Camera Number`),
    round = as.integer(`Round Number`),
    start_date = as.Date(`Start date`, format = "%d/%m/%Y"),
    end_date = as.Date(`End Date`, format = "%d/%m/%Y"),
    mid_date = start_date + floor(as.numeric(end_date - start_date) / 2)
  )

round_mid_lookup <- obs_covs %>%
  filter(!is.na(mid_date)) %>%
  group_by(round) %>%
  summarise(mid_num = median(as.numeric(mid_date)), .groups = "drop")

season_lookup <- obs_covs %>%
  left_join(round_mid_lookup, by = "round") %>%
  mutate(
    mid_date = as.Date(
      ifelse(is.na(mid_date), mid_num, as.numeric(mid_date)),
      origin = "1970-01-01"
    ),
    period = if_else(mid_date <= split_date, "Early (hotter)", "Late (mild-cold)")
  ) %>%
  select(camera_id, round, period)

if (anyDuplicated(season_lookup[c("camera_id", "round")])) {
  stop("Observation covariates contain duplicate camera-round records.")
}

count_data <- counts %>%
  left_join(season_lookup, by = c("camera_id", "round")) %>%
  left_join(site_groups, by = "camera_id") %>%
  filter(is.finite(fawn), is.finite(female))

if (anyNA(count_data$period) || anyNA(count_data$wolf_level) || anyNA(count_data$jackal_level)) {
  stop("At least one valid count record could not be matched to season or predator covariates.")
}

analysis_data <- bind_rows(
  count_data,
  count_data %>% mutate(period = "Overall")
) %>%
  mutate(
    period = factor(
      period,
      levels = c("Overall", "Early (hotter)", "Late (mild-cold)")
    )
  )

# Ratio of pooled counts: detection-rich camera-rounds contribute more information.
pooled <- analysis_data %>%
  group_by(period, wolf_level, jackal_level) %>%
  summarise(
    n_cameras = n_distinct(camera_id),
    n_camera_sessions = n(),
    fawn_count = sum(fawn),
    female_count = sum(female),
    ratio = if_else(female_count > 0, fawn_count / female_count, NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    method = "Pooled counts",
    label = case_when(
      n_cameras == 0 ~ "No cameras",
      female_count == 0 & fawn_count > 0 ~ sprintf("Inf\n(%d/0)\nC=%d", fawn_count, n_cameras),
      TRUE ~ sprintf("%.2f\n(%d/%d)\nC=%d", ratio, fawn_count, female_count, n_cameras)
    )
  )

# Give each informative session equal weight within camera, then each camera equal
# weight within a predator cell. The mean composition is converted back to a ratio.
camera_composition <- analysis_data %>%
  mutate(
    identified_total = fawn + female,
    fawn_proportion = if_else(identified_total > 0, fawn / identified_total, NA_real_)
  ) %>%
  filter(is.finite(fawn_proportion)) %>%
  group_by(period, wolf_level, jackal_level, camera_id) %>%
  summarise(
    camera_mean_proportion = mean(fawn_proportion),
    n_informative_sessions = n(),
    .groups = "drop"
  )

balanced <- camera_composition %>%
  group_by(period, wolf_level, jackal_level) %>%
  summarise(
    mean_fawn_proportion = mean(camera_mean_proportion),
    n_informative_cameras = n(),
    n_informative_sessions = sum(n_informative_sessions),
    ratio = if_else(
      mean_fawn_proportion < 1,
      mean_fawn_proportion / (1 - mean_fawn_proportion),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    method = "Equal camera/session weight",
    label = sprintf(
      "%.2f\nC=%d, S=%d",
      ratio,
      n_informative_cameras,
      n_informative_sessions
    )
  )

grid <- expand_grid(
  period = factor(
    c("Overall", "Early (hotter)", "Late (mild-cold)"),
    levels = c("Overall", "Early (hotter)", "Late (mild-cold)")
  ),
  wolf_level = factor(c("Low", "Middle", "High"), levels = c("Low", "Middle", "High")),
  jackal_level = factor(c("Low", "Middle", "High"), levels = c("Low", "Middle", "High"))
)

pooled_plot_data <- grid %>%
  left_join(pooled, by = c("period", "wolf_level", "jackal_level")) %>%
  mutate(
    method = "Pooled counts",
    label = replace_na(label, "No cameras")
  )

balanced_plot_data <- grid %>%
  left_join(balanced, by = c("period", "wolf_level", "jackal_level")) %>%
  mutate(
    method = "Equal camera/session weight",
    label = replace_na(label, "No identified\nfawn/female\ndetections")
  )

plot_data <- bind_rows(pooled_plot_data, balanced_plot_data) %>%
  mutate(
    method = factor(
      method,
      levels = c("Pooled counts", "Equal camera/session weight")
    ),
    text_color = if_else(is.finite(ratio) & ratio >= 0.22, "white", "#263238")
  )

write_csv(
  pooled_plot_data %>%
    mutate(across(c(period, wolf_level, jackal_level), as.character)),
  paste0(output_prefix, "_pooled_counts.csv")
)
write_csv(
  balanced_plot_data %>%
    mutate(across(c(period, wolf_level, jackal_level), as.character)),
  paste0(output_prefix, "_balanced_composition.csv")
)

finite_ratios <- plot_data$ratio[is.finite(plot_data$ratio)]
fill_upper <- max(0.4, max(finite_ratios))

ratio_plot <- ggplot(
  plot_data,
  aes(x = jackal_level, y = wolf_level, fill = ratio)
) +
  geom_tile(color = "#fbfaf7", linewidth = 1.1) +
  geom_text(
    aes(label = label, color = text_color),
    size = 2.75,
    lineheight = 0.9,
    fontface = "bold"
  ) +
  facet_grid(method ~ period) +
  scale_color_identity() +
  scale_fill_gradientn(
    colors = c("#f7f0df", "#e7c283", "#d17443", "#7e231f"),
    limits = c(0, fill_upper),
    oob = scales::squish,
    na.value = "#d8d8d4",
    name = "Fawn:female\ncount ratio"
  ) +
  coord_equal() +
  labs(
    x = "Jackal abundance tertile",
    y = "Wolf abundance tertile",
    title = "Observed fawn:female ratios across wolf-jackal contexts",
    subtitle = paste0(
      "Predator groups use 33/67% thresholds; seasonal assignment matches the spAbundance analysis\n",
      "Pooled labels show ratio (fawn/female counts) and cameras; balanced labels show ratio and informative cameras/sessions"
    ),
    caption = paste0(
      "Descriptive camera-trap count indices, not detection-corrected demographic ratios. ",
      "Early/late split date: ", split_date, "."
    )
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "#283336"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "#4c5a5e", lineheight = 1.05),
    plot.caption = element_text(color = "#5f696b", hjust = 0),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "#e8e1d4", color = NA),
    legend.position = "right"
  )

ggsave(paste0(output_prefix, "_combined.png"), ratio_plot, width = 12.8, height = 8.1, dpi = 320)
ggsave(paste0(output_prefix, "_combined.pdf"), ratio_plot, width = 12.8, height = 8.1)

cutoff_table <- tibble(
  predator = c("Wolf", "Jackal"),
  lower_33_cutoff = c(wolf_cutoffs[1], jackal_cutoffs[1]),
  upper_67_cutoff = c(wolf_cutoffs[2], jackal_cutoffs[2])
)
write_csv(cutoff_table, paste0(output_prefix, "_tertile_cutoffs.csv"))

cat("Saved count-based fawn:female heatmap outputs with prefix:", output_prefix, "\n")
