#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
arg_map <- list(
  input_file = "publication_predator_context_diel_pairwise_overlap_median_results.csv",
  output_prefix = "publication_predator_context_diel_pairwise_overlap_median",
  small_sample_threshold = "75",
  formal_minimum_events = "20",
  compare_reps = "999",
  camera_bootstrap_reps = "500"
)
for (arg in args) {
  parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    key <- sub("^--", "", parts[1])
    if (!key %in% names(arg_map)) stop("Unknown argument: ", key)
    arg_map[[key]] <- parts[2]
  }
}

input_file <- arg_map$input_file
output_prefix <- arg_map$output_prefix
small_sample_threshold <- as.integer(arg_map$small_sample_threshold)
formal_minimum_events <- as.integer(arg_map$formal_minimum_events)
compare_reps <- as.integer(arg_map$compare_reps)
camera_bootstrap_reps <- as.integer(arg_map$camera_bootstrap_reps)

class_levels <- c("Fawn", "Female", "Male")
context_levels <- c(
  "Low both",
  "High jackal / Low wolf",
  "Low jackal / High wolf",
  "High both"
)
context_short <- c(
  "Low both" = "Low both",
  "High jackal / Low wolf" = "High jackal / low wolf",
  "Low jackal / High wolf" = "Low jackal / high wolf",
  "High both" = "High both"
)

results <- read_csv(input_file, show_col_types = FALSE) %>%
  select(-any_of(c(
    "compareCkern_p_fdr_within_class",
    "compareCkern_p_fdr_global",
    "significance"
  ))) %>%
  mutate(
    significance = case_when(
      !formal_test_eligible ~ "Descriptive only",
      compareCkern_p < 0.001 ~ "***",
      compareCkern_p < 0.01 ~ "**",
      compareCkern_p < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(factor(class, levels = class_levels), overlap)

write_csv(results, paste0(output_prefix, "_results.csv"))

context_pairs <- combn(context_levels, 2, simplify = FALSE)
pair_levels <- unique(unlist(lapply(context_pairs, function(pair_i) {
  paste(
    unname(context_short[pair_i[1]]),
    "vs",
    unname(context_short[pair_i[2]])
  )
})))

plot_data <- results %>%
  mutate(
    class = factor(class, levels = class_levels),
    pair = factor(pair, levels = rev(pair_levels)),
    inference = case_when(
      !formal_test_eligible ~ "Descriptive only (<20 events)",
      compareCkern_p < 0.05 ~ "Raw P < 0.05",
      TRUE ~ "Not supported"
    ),
    estimate_label = paste0(sprintf("%.2f", overlap), " ", significance)
  )

inference_colors <- c(
  "Raw P < 0.05" = "#A52A2A",
  "Not supported" = "#245B78",
  "Descriptive only (<20 events)" = "#8A8F91"
)

overlap_plot <- ggplot(plot_data, aes(x = overlap, y = pair, color = inference)) +
  geom_segment(
    aes(x = conf_low, xend = conf_high, yend = pair),
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  geom_point(size = 2.8) +
  geom_text(
    aes(label = estimate_label),
    hjust = -0.12,
    size = 3.0,
    color = "#263238"
  ) +
  facet_wrap(~class, ncol = 1) +
  scale_color_manual(values = inference_colors, name = NULL) +
  scale_x_continuous(
    limits = c(0, 1.08),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Temporal overlap (Dhat)",
    y = NULL,
    title = "Within-class temporal overlap among wolf-jackal contexts",
    subtitle = paste0(
      "Intervals use camera-cluster bootstrap; tests use unadjusted compareCkern P values\n",
      "Dhat1 used when the smaller sample was <", small_sample_threshold,
      " events; otherwise Dhat4"
    ),
    caption = paste0(
      "* P<0.05; ** P<0.01; *** P<0.001. P values are unadjusted. Tests with fewer than ",
      formal_minimum_events,
      " events in either context are descriptive only. Lower overlap indicates a larger temporal shift."
    )
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(color = "#304044"),
    axis.title.x = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "#E8E1D4", color = NA),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "#4C5A5E", lineheight = 1.05),
    plot.caption = element_text(color = "#5F696B", hjust = 0),
    legend.position = "top"
  )

ggsave(paste0(output_prefix, ".png"), overlap_plot, width = 10.6, height = 10.2, dpi = 320)
ggsave(paste0(output_prefix, ".pdf"), overlap_plot, width = 10.6, height = 10.2)

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}
fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", fmt(x, 3)))
}

table_rows <- vapply(seq_len(nrow(results)), function(i) {
  row <- results[i, ]
  interval_text <- if (is.finite(row$conf_low)) {
    paste0(fmt(row$overlap), " [", fmt(row$conf_low), ", ", fmt(row$conf_high), "]")
  } else {
    paste0(fmt(row$overlap), " [NA]")
  }
  paste0(
    "| ", row$class,
    " | ", row$pair,
    " | `", row$n_a, " / ", row$n_b,
    "` | `", row$estimator,
    "` | `", interval_text,
    "` | `", fmt_p(row$compareCkern_p), "` |"
  )
}, character(1))

lowest_rows <- unlist(lapply(class_levels, function(class_i) {
  lowest <- results %>%
    filter(class == class_i, formal_test_eligible, is.finite(overlap)) %>%
    slice_min(overlap, n = 1, with_ties = FALSE)
  if (!nrow(lowest)) return(character())
  paste0(
    "- **", class_i, ":** ", lowest$pair,
    ", `Dhat = ", fmt(lowest$overlap),
    "` (95% camera-bootstrap CI `", fmt(lowest$conf_low), "-", fmt(lowest$conf_high),
    "`); unadjusted `compareCkern P = ", fmt_p(lowest$compareCkern_p), "`."
  )
}))

md <- c(
  "# Pairwise predator-context temporal overlap within gazelle classes",
  "",
  "## Methods",
  "",
  paste0(
    "All six context pairs were evaluated within each class using the all-camera median classification. `Dhat1` was used when the smaller sample contained fewer than ",
    small_sample_threshold,
    " events; otherwise `Dhat4` was used. Intervals were obtained from ",
    camera_bootstrap_reps,
    " camera-cluster bootstrap replicates. Distributional differences were tested with `activity::compareCkern` using ",
    compare_reps,
    " null replicates and a +1 Monte Carlo correction. We report unadjusted P-values."
  ),
  "",
  paste0(
    "Because `compareCkern` resamples events rather than cameras, its P-values do not fully account for camera-level clustering. Comparisons with fewer than ",
    formal_minimum_events,
    " events in either context were calculated descriptively but excluded from formal inference."
  ),
  "",
  "## Results",
  "",
  "| Class | Context pair | Events | Estimator | Dhat [camera-bootstrap 95% CI] | compareCkern P |",
  "|---|---|---:|---|---:|---:|",
  table_rows,
  "",
  "## Lowest eligible overlap by class",
  "",
  lowest_rows
)
writeLines(md, paste0(output_prefix, "_results.md"), useBytes = TRUE)

cat("Rendered raw-P overlap outputs with prefix:", output_prefix, "\n")
