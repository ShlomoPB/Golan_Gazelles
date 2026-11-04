if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required.")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required.")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  stop("Package 'readr' is required.")
}

library(ggplot2)
library(dplyr)
library(readr)

fawn_file <- "fawn_NMix_Open_Spaces_jackal_yes_wolf_yes_extra_null_boar_abundance_human_abundance_Null_Wolf_Abundance_det_effort_NB_summary.txt"
male_file <- "male_NMix_Built_and_Roads_jackal_yes_wolf_yes_extra_null_boar_abundance_human_abundance_Null_Wolf_Abundance_det_effort_NB_summary.txt"
female_file <- "female_NMix_Nature_Reserves_jackal_yes_wolf_yes_extra_null_boar_abundance_human_abundance_Null_Wolf_Abundance_det_effort_NB_summary.txt"

read_abundance_block <- function(path, class_label, habitat_term) {
  lines <- readLines(path, warn = FALSE)
  start <- grep("^Abundance \\(log scale\\):", lines)
  end <- grep("^Detection \\(logit scale\\):", lines)
  if (length(start) == 0 || length(end) == 0) {
    stop("Could not parse abundance block in ", path)
  }
  block <- lines[(start + 2):(end - 2)]
  block <- block[nzchar(trimws(block))]
  block <- gsub("^\\s+", "", block)
  con <- textConnection(paste(block, collapse = "\n"))
  on.exit(close(con), add = TRUE)
  df <- read.table(con, header = FALSE, stringsAsFactors = FALSE)
  names(df) <- c("term", "mean", "sd", "lwr", "median", "upr", "rhat", "ess")
  df$class <- class_label
  df$primary_habitat <- habitat_term
  df
}

coef_df <- bind_rows(
  read_abundance_block(fawn_file, "Fawn", "Open.Spaces.z"),
  read_abundance_block(male_file, "Male", "Built.and.Roads.z"),
  read_abundance_block(female_file, "Female", "Nature.Reserves.z")
)

term_labels <- c(
  "Open.Spaces.z" = "Open spaces",
  "Built.and.Roads.z" = "Built and roads",
  "Nature.Reserves.z" = "Nature reserves",
  "null.jackal.abundance.z" = "Jackal abundance",
  "Null.Wolf.Abundance.z" = "Wolf abundance",
  "null.boar.abundance.z" = "Boar abundance",
  "human.abundance.z" = "Human abundance"
)

term_order <- c(
  "Open.Spaces.z",
  "Built.and.Roads.z",
  "Nature.Reserves.z",
  "null.jackal.abundance.z",
  "Null.Wolf.Abundance.z",
  "null.boar.abundance.z",
  "human.abundance.z"
)

plot_df <- coef_df %>%
  filter(term %in% names(term_labels)) %>%
  mutate(
    term_label = factor(term_labels[term], levels = term_labels[term_order]),
    highlight = ifelse(term == primary_habitat | term %in% c("null.jackal.abundance.z", "Null.Wolf.Abundance.z"),
                       "Key effect", "Other effect"),
    supported = ifelse(lwr > 0 | upr < 0, "CrI excludes 0", "CrI overlaps 0")
  )

p <- ggplot(plot_df, aes(x = mean, y = term_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.5) +
  geom_errorbar(
    aes(xmin = lwr, xmax = upr, color = supported),
    width = 0.18,
    linewidth = 0.9,
    orientation = "y"
  ) +
  geom_point(aes(shape = highlight, fill = highlight), size = 3.2, color = "black") +
  facet_wrap(~ class, ncol = 3, scales = "free_y") +
  scale_color_manual(values = c("CrI excludes 0" = "#0b6e4f", "CrI overlaps 0" = "#9f9f9f")) +
  scale_fill_manual(values = c("Key effect" = "#d95f02", "Other effect" = "#fdd0a2")) +
  scale_shape_manual(values = c("Key effect" = 21, "Other effect" = 21)) +
  labs(
    x = "Posterior mean effect on log abundance",
    y = NULL,
    color = NULL,
    fill = NULL,
    shape = NULL,
    title = "Class-specific responses to habitat, predator, and anthropogenic covariates"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave("publication_class_coefficient_figure.png", p, width = 10.5, height = 5.8, dpi = 300)
ggsave("publication_class_coefficient_figure.pdf", p, width = 10.5, height = 5.8)

write_csv(
  plot_df %>%
    transmute(
      class,
      term,
      term_label = as.character(term_label),
      mean,
      lwr,
      upr,
      highlight,
      supported
    ),
  "publication_class_coefficient_figure_data.csv"
)

cat("Saved figure to publication_class_coefficient_figure.png and .pdf\n")
