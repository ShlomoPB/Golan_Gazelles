#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

fit_file <- "publication_wolf_jackal_interaction_nmix_spabundance_fits.rds"
output_prefix <- "publication_spabundance_wolf_jackal_relative"

fits <- readRDS(fit_file)
class_labels <- c(fawn = "Juvenile", female = "Female", male = "Male")
wolf_levels <- c(-1, 0, 1)
wolf_labels <- c(
  `-1` = "Low wolf activity (-1 SD)",
  `0` = "Mean wolf activity",
  `1` = "High wolf activity (+1 SD)"
)

summarize_grid <- function(fit, grid) {
  draws <- as.matrix(fit$beta.samples)
  beta <- draws[, c(
    "Null.Wolf.Abundance.z",
    "null.jackal.abundance.z",
    "Null.Wolf.Abundance.z:null.jackal.abundance.z"
  ), drop = FALSE]
  design <- cbind(
    grid$wolf_z,
    grid$jackal_z,
    grid$wolf_z * grid$jackal_z
  )
  relative_draws <- exp(beta %*% t(design))
  data.frame(
    grid,
    estimate = apply(relative_draws, 2, median),
    conf_low = apply(relative_draws, 2, quantile, probs = 0.025),
    conf_high = apply(relative_draws, 2, quantile, probs = 0.975)
  )
}

line_grid <- expand.grid(
  wolf_z = wolf_levels,
  jackal_z = seq(-1.5, 1.5, length.out = 121)
)
heat_grid <- expand.grid(wolf_z = wolf_levels, jackal_z = wolf_levels)

line_results <- bind_rows(lapply(names(fits), function(class_name) {
  summarize_grid(fits[[class_name]], line_grid) %>%
    mutate(class = unname(class_labels[class_name]))
})) %>%
  mutate(
    class = factor(class, levels = unname(class_labels)),
    wolf_level = factor(
      wolf_labels[as.character(wolf_z)],
      levels = unname(wolf_labels)
    )
  )

heat_results <- bind_rows(lapply(names(fits), function(class_name) {
  summarize_grid(fits[[class_name]], heat_grid) %>%
    mutate(class = unname(class_labels[class_name]))
})) %>%
  mutate(
    class = factor(class, levels = unname(class_labels)),
    is_reference = wolf_z == 0 & jackal_z == 0,
    uncertainty_label = if_else(
      is_reference,
      "(reference)",
      sprintf("(%.2f-%.2f)", conf_low, conf_high)
    )
  )

conditional_results <- bind_rows(lapply(names(fits), function(class_name) {
  draws <- as.matrix(fits[[class_name]]$beta.samples)
  beta_wolf <- draws[, "Null.Wolf.Abundance.z"]
  beta_jackal <- draws[, "null.jackal.abundance.z"]
  beta_interaction <- draws[, "Null.Wolf.Abundance.z:null.jackal.abundance.z"]
  bind_rows(
    lapply(wolf_levels, function(jackal_value) {
      ratio <- exp(beta_wolf + jackal_value * beta_interaction)
      data.frame(
        class = unname(class_labels[class_name]),
        focal_predator = "Wolf",
        conditioning_predator = "Jackal",
        conditioning_value_sd = jackal_value,
        mean = mean(ratio),
        median = median(ratio),
        conf_low = quantile(ratio, 0.025),
        conf_high = quantile(ratio, 0.975)
      )
    }),
    lapply(wolf_levels, function(wolf_value) {
      ratio <- exp(beta_jackal + wolf_value * beta_interaction)
      data.frame(
        class = unname(class_labels[class_name]),
        focal_predator = "Jackal",
        conditioning_predator = "Wolf",
        conditioning_value_sd = wolf_value,
        mean = mean(ratio),
        median = median(ratio),
        conf_low = quantile(ratio, 0.025),
        conf_high = quantile(ratio, 0.975)
      )
    })
  )
}))

write_csv(line_results, paste0(output_prefix, "_interaction_curves.csv"))
write_csv(heat_results, paste0(output_prefix, "_prediction_grid.csv"))
write_csv(
  conditional_results,
  paste0(output_prefix, "_conditional_effect_ratios.csv")
)

palette <- c(
  "Low wolf activity (-1 SD)" = "#B34135",
  "Mean wolf activity" = "#68747A",
  "High wolf activity (+1 SD)" = "#173F5F"
)

curve_plot <- ggplot(
  line_results,
  aes(x = jackal_z, y = estimate, color = wolf_level, fill = wolf_level)
) +
  geom_hline(yintercept = 1, color = "#AAB2B5", linewidth = 0.45) +
  geom_ribbon(
    aes(ymin = conf_low, ymax = conf_high),
    alpha = 0.10,
    color = NA
  ) +
  geom_line(linewidth = 1.05) +
  facet_wrap(~class, nrow = 1) +
  scale_y_log10() +
  scale_color_manual(values = palette, name = NULL) +
  scale_fill_manual(values = palette, name = NULL) +
  labs(
    x = "Jackal activity (standardized)",
    y = "Relative expected abundance",
    title = "Class-specific relative abundance across wolf-jackal gradients",
    subtitle = "Posterior medians and 95% credible bands; reference = mean wolf and jackal activity",
    caption = "The relative scale excludes the abundance intercept; curves describe gradients rather than absolute population size."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "top"
  )

ggsave(
  paste0(output_prefix, "_interaction_curves.png"),
  curve_plot,
  width = 9.4,
  height = 4.8,
  dpi = 320
)
ggsave(
  paste0(output_prefix, "_interaction_curves.pdf"),
  curve_plot,
  width = 9.4,
  height = 4.8
)

heat_plot <- ggplot(
  heat_results,
  aes(x = factor(jackal_z), y = factor(wolf_z), fill = log(estimate))
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = sprintf("%.2f", estimate)),
    nudge_y = 0.10,
    size = 3.85,
    fontface = "bold"
  ) +
  geom_text(
    aes(label = uncertainty_label),
    nudge_y = -0.12,
    size = 2.65
  ) +
  facet_wrap(~class, nrow = 1) +
  scale_fill_gradient2(
    low = "#C85C4A",
    mid = "#F6F3EA",
    high = "#204D69",
    midpoint = 0,
    breaks = log(c(0.25, 0.5, 1, 2, 4)),
    labels = c("0.25", "0.5", "1", "2", "4"),
    name = "Gazelle relative\nactivity"
  ) +
  labs(
    x = "Jackal activity (SD)",
    y = "Wolf activity (SD)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 15)
  )

ggsave(
  paste0(output_prefix, "_heatmap.png"),
  heat_plot,
  width = 10.34,
  height = 4.8,
  dpi = 320
)
ggsave(
  paste0(output_prefix, "_heatmap.pdf"),
  heat_plot,
  width = 10.34,
  height = 4.8
)
