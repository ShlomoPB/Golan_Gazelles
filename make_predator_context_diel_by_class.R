#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
arg_map <- list(
  source_prefix = "publication_predator_context_class_diel_2x2_median",
  output_prefix = "publication_predator_context_diel_by_class_median",
  event_file = "Gazella gazella_independent_events_30min_analysis_ready.csv",
  obs_cov_file = "ObsCovs_050524.csv",
  local_timezone = "Asia/Jerusalem",
  minimum_curve_events = "5"
)
for (arg in args) {
  parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    key <- sub("^--", "", parts[1])
    if (!key %in% names(arg_map)) stop("Unknown argument: ", key)
    arg_map[[key]] <- parts[2]
  }
}

source_prefix <- arg_map$source_prefix
output_prefix <- arg_map$output_prefix
event_file <- arg_map$event_file
obs_cov_file <- arg_map$obs_cov_file
local_timezone <- arg_map$local_timezone
minimum_curve_events <- as.integer(arg_map$minimum_curve_events)

density_file <- paste0(source_prefix, "_density_curves.csv")
count_file <- paste0(source_prefix, "_event_counts.csv")
camera_file <- paste0(source_prefix, "_camera_categories.csv")

required_files <- c(density_file, count_file, camera_file, event_file, obs_cov_file)
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
context_abbreviations <- c(
  "Low both" = "Low both",
  "High jackal / Low wolf" = "High jackal / low wolf",
  "Low jackal / High wolf" = "Low jackal / high wolf",
  "High both" = "High both"
)
context_colors <- c(
  "Low both" = "#188977",
  "High jackal / Low wolf" = "#B12B32",
  "Low jackal / High wolf" = "#173B63",
  "High both" = "#8A5A9E"
)
class_labels <- c(
  juvenile = "Fawn",
  female = "Female",
  male = "Male"
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

density_curves <- read_csv(density_file, show_col_types = FALSE) %>%
  transmute(
    class = factor(class, levels = class_levels),
    predator_context = factor(as.character(predator_context), levels = context_levels),
    solar_hour,
    density
  )

event_counts <- read_csv(count_file, show_col_types = FALSE) %>%
  mutate(
    class = factor(class, levels = class_levels),
    predator_context = factor(as.character(predator_context), levels = context_levels)
  ) %>%
  arrange(class, predator_context)

camera_context <- read_csv(camera_file, show_col_types = FALSE) %>%
  transmute(
    camera_name = as.integer(camera_name),
    predator_context = factor(as.character(predator_context), levels = context_levels)
  )

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
    solar_hour = apparent_solar_hour(timestamp, camelot_gps_longitude)
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
  filter(class %in% class_levels, predator_context %in% context_levels) %>%
  mutate(
    class = factor(class, levels = class_levels),
    predator_context = factor(predator_context, levels = context_levels)
  )

sparse_groups <- event_counts %>%
  filter(n_events < minimum_curve_events)
sparse_events <- events %>%
  inner_join(
    sparse_groups %>% select(class, predator_context),
    by = c("class", "predator_context")
  )

light_background <- tibble(
  xmin = c(0, 5, 7, 17, 19),
  xmax = c(5, 7, 17, 19, 24),
  light_phase = c("Night", "Dawn", "Day", "Dusk", "Night")
)
combined_background <- crossing(
  class = factor(class_levels, levels = class_levels),
  light_background
)

count_annotations <- event_counts %>%
  mutate(
    context_short = unname(context_abbreviations[as.character(predator_context)]),
    context_line = paste0(
      context_short,
      ": n=",
      n_events,
      ifelse(n_events < minimum_curve_events, " (data deficient)", "")
    )
  ) %>%
  group_by(class) %>%
  summarise(label = paste(context_line, collapse = "\n"), .groups = "drop") %>%
  mutate(x = 23.65, y = Inf)

base_layers <- list(
  scale_x_continuous(
    breaks = seq(0, 24, by = 3),
    limits = c(0, 24),
    expand = expansion(mult = c(0, 0))
  ),
  scale_color_manual(values = context_colors, drop = FALSE, name = NULL),
  scale_fill_manual(
    values = c(
      Night = "#294C60",
      Dawn = "#E9A45D",
      Day = "#F5E6B8",
      Dusk = "#CF6A4C"
    ),
    guide = "none"
  ),
  labs(
    x = "Apparent solar time",
    y = "Relative activity density",
    caption = paste0(
      "All 60 cameras classified by median wolf and jackal activity. Curves are normalized distributions, not event rates.\n",
      "Camera clocks were interpreted as Israel local civil time and converted to apparent solar time; shading indicates approximate light periods."
    )
  ),
  theme_minimal(base_size = 11),
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "#314043"),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "#E8E1D4", color = NA),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "#4C5A5E", lineheight = 1.05),
    plot.caption = element_text(color = "#5F696B", hjust = 0),
    legend.position = "top"
  )
)

combined_plot <- ggplot(
  density_curves,
  aes(x = solar_hour, y = density, color = predator_context)
) +
  geom_rect(
    data = combined_background,
    inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = light_phase),
    alpha = 0.16,
    color = NA
  ) +
  geom_line(linewidth = 1.0) +
  geom_rug(
    data = sparse_events,
    inherit.aes = FALSE,
    aes(x = solar_hour, color = predator_context),
    sides = "b",
    linewidth = 0.9,
    alpha = 0.95
  ) +
  geom_text(
    data = count_annotations,
    inherit.aes = FALSE,
    aes(x = x, y = y, label = label),
    hjust = 1,
    vjust = 1.1,
    size = 2.75,
    lineheight = 0.98,
    color = "#263238"
  ) +
  facet_wrap(~class, ncol = 1, scales = "free_y") +
  labs(
    title = "Class-specific diel activity across wolf-jackal contexts",
    subtitle = paste0(
      "Each panel compares the four median-split predator contexts within one gazelle class\n",
      "The high-jackal/low-wolf fawn sample is data deficient (n=3) and is shown as event ticks without a density curve"
    )
  ) +
  base_layers

ggsave(
  paste0(output_prefix, "_combined.png"),
  combined_plot,
  width = 10.6,
  height = 11.8,
  dpi = 320
)
ggsave(
  paste0(output_prefix, "_combined.pdf"),
  combined_plot,
  width = 10.6,
  height = 11.8
)

make_class_plot <- function(class_i) {
  class_density <- filter(density_curves, class == class_i) %>%
    left_join(
      event_counts %>%
        filter(class == class_i) %>%
        select(predator_context, n_events),
      by = "predator_context"
    ) %>%
    mutate(
      legend_label = paste0(
        unname(context_abbreviations[as.character(predator_context)]),
        " (n=", n_events, ")"
      )
    )
  label_lookup <- event_counts %>%
    filter(class == class_i) %>%
    transmute(
      predator_context,
      legend_label = paste0(
        unname(context_abbreviations[as.character(predator_context)]),
        " (n=", n_events,
        ifelse(n_events < minimum_curve_events, ", data deficient", ""),
        ")"
      )
    )
  legend_levels <- label_lookup$legend_label
  class_density <- class_density %>%
    select(-legend_label) %>%
    left_join(label_lookup, by = "predator_context") %>%
    mutate(legend_label = factor(legend_label, levels = legend_levels))
  class_sparse <- sparse_events %>%
    filter(class == class_i) %>%
    left_join(label_lookup, by = "predator_context") %>%
    mutate(legend_label = factor(legend_label, levels = legend_levels))

  legend_values <- setNames(
    unname(context_colors[as.character(label_lookup$predator_context)]),
    label_lookup$legend_label
  )

  ggplot(class_density, aes(x = solar_hour, y = density, color = legend_label)) +
    geom_rect(
      data = light_background,
      inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = light_phase),
      alpha = 0.16,
      color = NA
    ) +
    geom_line(linewidth = 1.15) +
    geom_rug(
      data = class_sparse,
      inherit.aes = FALSE,
      aes(x = solar_hour, color = legend_label),
      sides = "b",
      linewidth = 1.0,
      alpha = 0.95
    ) +
    scale_color_manual(
      values = legend_values,
      breaks = legend_levels,
      drop = FALSE,
      name = NULL
    ) +
    scale_x_continuous(
      breaks = seq(0, 24, by = 3),
      limits = c(0, 24),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_fill_manual(
      values = c(
        Night = "#294C60",
        Dawn = "#E9A45D",
        Day = "#F5E6B8",
        Dusk = "#CF6A4C"
      ),
      guide = "none"
    ) +
    labs(
      x = "Apparent solar time",
      y = "Relative activity density",
      title = paste0(class_i, " diel activity across wolf-jackal contexts"),
      subtitle = "All-camera median split; curves are normalized within predator context",
      caption = paste0(
        "Independent 30-min events within active camera rounds. Curves describe timing, not abundance.\n",
        "Associations across camera contexts should not be interpreted as causal temporal displacement."
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "#314043"),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(color = "#4C5A5E"),
      plot.caption = element_text(color = "#5F696B", hjust = 0),
      legend.position = "top",
      plot.margin = margin(12, 18, 12, 22)
    )
}

for (class_i in class_levels) {
  class_plot <- make_class_plot(class_i)
  class_slug <- tolower(class_i)
  ggsave(
    paste0(output_prefix, "_", class_slug, ".png"),
    class_plot,
    width = 9.3,
    height = 6.3,
    dpi = 320
  )
  ggsave(
    paste0(output_prefix, "_", class_slug, ".pdf"),
    class_plot,
    width = 9.3,
    height = 6.3
  )
}

write_csv(
  event_counts,
  paste0(output_prefix, "_event_counts.csv")
)

md <- c(
  "# Class-specific predator-context diel figures",
  "",
  "The figures compare normalized apparent-solar-time activity distributions among four median-split wolf-jackal contexts separately for fawns, adult females, and adult males. All 60 cameras were assigned to a context. Curves use 30-min independent events within documented active camera rounds.",
  "",
  "## Interpretation constraint",
  "",
  "Differences among curves describe associations between camera-level predator context and gazelle activity timing. They do not demonstrate that predators caused temporal displacement because habitat, camera placement, and other site differences may also vary among contexts.",
  "",
  "The high-jackal/low-wolf fawn context was data deficient, containing only three events from two cameras. Those times are shown as ticks, and no density curve was estimated.",
  "",
  "## Outputs",
  "",
  paste0("- Combined: `", output_prefix, "_combined.png` and `.pdf`."),
  paste0("- Fawn: `", output_prefix, "_fawn.png` and `.pdf`."),
  paste0("- Female: `", output_prefix, "_female.png` and `.pdf`."),
  paste0("- Male: `", output_prefix, "_male.png` and `.pdf`.")
)
writeLines(md, paste0(output_prefix, "_results.md"), useBytes = TRUE)

cat("Saved class-specific predator-context diel figures with prefix:", output_prefix, "\n")
