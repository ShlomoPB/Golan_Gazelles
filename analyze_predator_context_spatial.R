if (!requireNamespace("readr", quietly = TRUE)) {
  stop("Package 'readr' is required.")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required.")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  stop("Package 'tidyr' is required.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required.")
}
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required.")
}

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

juvenile_file <- "gazelle_juvenile_abundance_60x7.csv"
female_file <- "gazelle_female_abundance_60x7.csv"
male_file <- "gazelle_male_abundance_60x7.csv"
category_file <- "publication_predator_context_q33_67_camera_categories.csv"
coords_file <- "Camera Coordinates.xlsx"
output_prefix <- "publication_predator_context_q33_67_spatial"
set.seed(20260623)

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

read_count_matrix <- function(path) {
  x <- read.csv(path, na.strings = c("NA", ""), check.names = FALSE)
  y <- as.matrix(x[paste0("X", 1:7)])
  storage.mode(y) <- "numeric"
  rownames(y) <- x$camera_name
  colnames(y) <- paste0("X", 1:7)
  y
}

fmt_num <- function(x, digits = 2) formatC(x, format = "f", digits = digits)

juvenile_mat <- read_count_matrix(juvenile_file)
female_mat <- read_count_matrix(female_file)
male_mat <- read_count_matrix(male_file)

stopifnot(identical(rownames(juvenile_mat), rownames(female_mat)))
stopifnot(identical(rownames(juvenile_mat), rownames(male_mat)))

camera_ids <- rownames(juvenile_mat)
active_rounds <- rowSums(is.finite(female_mat))

camera_df <- tibble(
  camera_name = as.integer(camera_ids),
  juvenile_total = rowSums(juvenile_mat, na.rm = TRUE),
  female_total = rowSums(female_mat, na.rm = TRUE),
  male_total = rowSums(male_mat, na.rm = TRUE),
  active_rounds = active_rounds,
  juvenile_per_active_round = juvenile_total / active_rounds,
  female_per_active_round = female_total / active_rounds,
  male_per_active_round = male_total / active_rounds,
  any_juvenile = juvenile_total > 0,
  juvenile_female_ratio = ifelse(female_total > 0, juvenile_total / female_total, NA_real_),
  juveniles_per_100_females = ifelse(female_total > 0, 100 * juvenile_total / female_total, NA_real_)
)

category_df <- read_csv(category_file, show_col_types = FALSE) %>%
  select(camera_name, Polygon.Name, predator_context_label)

coords_raw <- read_excel(coords_file)
names(coords_raw) <- c("site_label", "lat", "lon")
coords_raw$camera_name <- as.integer(gsub("[^0-9]", "", as.character(coords_raw$site_label)))
coords_raw <- coords_raw[, c("camera_name", "lat", "lon")]

camera_df <- camera_df %>%
  left_join(category_df, by = "camera_name") %>%
  left_join(coords_raw, by = "camera_name") %>%
  mutate(
    predator_context_label = factor(predator_context_label, levels = context_levels)
  )

category_summary <- camera_df %>%
  group_by(predator_context_label) %>%
  summarise(
    n_cameras = n(),
    juvenile_total = sum(juvenile_total),
    female_total = sum(female_total),
    male_total = sum(male_total),
    mean_juvenile_per_camera = mean(juvenile_total),
    mean_juvenile_per_active_round = mean(juvenile_per_active_round),
    median_juvenile_per_active_round = median(juvenile_per_active_round),
    cameras_with_any_juvenile = sum(any_juvenile),
    juvenile_camera_occupancy = mean(any_juvenile),
    juveniles_per_100_females = 100 * sum(juvenile_total) / sum(female_total),
    .groups = "drop"
  )

hjlw_df <- filter(camera_df, predator_context_label == "High jackal / Low wolf")
n_hjlw <- nrow(hjlw_df)
obs_total <- sum(hjlw_df$juvenile_total)
obs_mean <- mean(hjlw_df$juvenile_per_active_round)
obs_occ <- sum(hjlw_df$any_juvenile)

perm_n <- 20000
perm_total <- numeric(perm_n)
perm_mean <- numeric(perm_n)
perm_occ <- numeric(perm_n)
for (i in seq_len(perm_n)) {
  idx <- sample(seq_len(nrow(camera_df)), size = n_hjlw, replace = FALSE)
  sub <- camera_df[idx, ]
  perm_total[i] <- sum(sub$juvenile_total)
  perm_mean[i] <- mean(sub$juvenile_per_active_round)
  perm_occ[i] <- sum(sub$any_juvenile)
}

null_summary <- tibble(
  metric = c("juvenile_total", "juvenile_per_active_round_mean", "cameras_with_any_juvenile"),
  observed = c(obs_total, obs_mean, obs_occ),
  null_mean = c(mean(perm_total), mean(perm_mean), mean(perm_occ)),
  null_median = c(median(perm_total), median(perm_mean), median(perm_occ)),
  null_q025 = c(quantile(perm_total, 0.025), quantile(perm_mean, 0.025), quantile(perm_occ, 0.025)),
  null_q975 = c(quantile(perm_total, 0.975), quantile(perm_mean, 0.975), quantile(perm_occ, 0.975)),
  p_lower = c(mean(perm_total <= obs_total), mean(perm_mean <= obs_mean), mean(perm_occ <= obs_occ))
)

write_csv(camera_df, paste0(output_prefix, "_camera_summary.csv"))
write_csv(category_summary, paste0(output_prefix, "_category_summary.csv"))
write_csv(null_summary, paste0(output_prefix, "_high_jackal_low_wolf_null_summary.csv"))

map_fawn <- ggplot(camera_df, aes(x = lon, y = lat)) +
  geom_point(aes(color = predator_context_label), size = 2.4, alpha = 0.45) +
  geom_point(aes(size = juvenile_total, fill = predator_context_label), shape = 21, color = "black", stroke = 0.25, alpha = 0.9) +
  scale_color_manual(values = context_colors, drop = FALSE, guide = "none") +
  scale_fill_manual(values = context_colors, drop = FALSE) +
  scale_size_area(max_size = 12) +
  coord_equal() +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = NULL,
    size = "Total juvenile count",
    title = "Spatial distribution of juvenile detections across predator-context cameras"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(paste0(output_prefix, "_juvenile_map.png"), map_fawn, width = 8.5, height = 6, dpi = 300)
ggsave(paste0(output_prefix, "_juvenile_map.pdf"), map_fawn, width = 8.5, height = 6)

map_long <- camera_df %>%
  select(camera_name, lon, lat, predator_context_label, juvenile_total, female_total, male_total) %>%
  pivot_longer(
    cols = c(juvenile_total, female_total, male_total),
    names_to = "class",
    values_to = "count"
  ) %>%
  mutate(
    class = recode(class,
      juvenile_total = "Fawn",
      female_total = "Female",
      male_total = "Male"
    )
  )

map_classes <- ggplot(map_long, aes(x = lon, y = lat)) +
  geom_point(aes(color = predator_context_label), size = 1.8, alpha = 0.35) +
  geom_point(aes(size = count, fill = predator_context_label), shape = 21, color = "black", stroke = 0.2, alpha = 0.9) +
  scale_color_manual(values = context_colors, drop = FALSE, guide = "none") +
  scale_fill_manual(values = context_colors, drop = FALSE) +
  scale_size_area(max_size = 10) +
  facet_wrap(~ class, ncol = 3) +
  coord_equal() +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = NULL,
    size = "Total count",
    title = "Camera-level class counts across predator-context categories"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(paste0(output_prefix, "_class_maps.png"), map_classes, width = 12, height = 4.6, dpi = 300)
ggsave(paste0(output_prefix, "_class_maps.pdf"), map_classes, width = 12, height = 4.6)

bar_plot <- ggplot(category_summary, aes(x = predator_context_label, y = mean_juvenile_per_active_round, fill = predator_context_label)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = fmt_num(mean_juvenile_per_active_round, 2)), vjust = -0.35, size = 3.1) +
  scale_fill_manual(values = context_colors, guide = "none") +
  labs(
    x = NULL,
    y = "Mean juvenile count per active round",
    title = "Juvenile camera-level abundance by predator context"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 18, hjust = 1)
  )

ggsave(paste0(output_prefix, "_juvenile_bar.png"), bar_plot, width = 8.5, height = 4.8, dpi = 300)
ggsave(paste0(output_prefix, "_juvenile_bar.pdf"), bar_plot, width = 8.5, height = 4.8)

md <- c(
  "# Spatial analysis of predator-context categories",
  "",
  "## Category summaries",
  "",
  "| Category | Cameras | Juvenile total | Female total | Male total | Mean juvenile / active round | Cameras with juveniles | Juveniles per 100 females |",
  "|---|---:|---:|---:|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(category_summary))) {
  md <- c(
    md,
    paste0(
      "| `", category_summary$predator_context_label[i], "` | `", category_summary$n_cameras[i],
      "` | `", category_summary$juvenile_total[i], "` | `", category_summary$female_total[i],
      "` | `", category_summary$male_total[i], "` | `", fmt_num(category_summary$mean_juvenile_per_active_round[i], 2),
      "` | `", category_summary$cameras_with_any_juvenile[i], "/", category_summary$n_cameras[i],
      "` | `", fmt_num(category_summary$juveniles_per_100_females[i], 1), "` |"
    )
  )
}

md <- c(
  md,
  "",
  "## Null comparison for `High jackal / Low wolf` cameras",
  "",
  paste0("Observed `High jackal / Low wolf` camera count: `", n_hjlw, "`."),
  "",
  "| Metric | Observed | Null mean | Null 95% interval | Lower-tail p |",
  "|---|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(null_summary))) {
  md <- c(
    md,
    paste0(
      "| `", null_summary$metric[i], "` | `", fmt_num(null_summary$observed[i], 3),
      "` | `", fmt_num(null_summary$null_mean[i], 3),
      "` | `", fmt_num(null_summary$null_q025[i], 3), " to ", fmt_num(null_summary$null_q975[i], 3),
      "` | `", fmt_num(null_summary$p_lower[i], 3), "` |"
    )
  )
}

writeLines(md, paste0(output_prefix, "_results.md"))

cat("Saved predator-context spatial outputs with prefix:", output_prefix, "\n")
