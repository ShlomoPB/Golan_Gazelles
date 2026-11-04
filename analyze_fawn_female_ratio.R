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
if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  stop("Package 'glmmTMB' is required.")
}

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)

juvenile_file <- "gazelle_juvenile_abundance_60x7.csv"
female_file <- "gazelle_female_abundance_60x7.csv"
site_cov_file <- "SiteCovs_060324.csv"
output_prefix <- "publication_fawn_female_ratio"

read_count_matrix <- function(path) {
  x <- read.csv(path, na.strings = c("NA", ""), check.names = FALSE)
  y <- as.matrix(x[paste0("X", 1:7)])
  storage.mode(y) <- "numeric"
  rownames(y) <- x$camera_name
  colnames(y) <- paste0("X", 1:7)
  y
}

z_scale <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / s
}

fmt_num <- function(x, digits = 2) formatC(x, format = "f", digits = digits)
fmt_ci <- function(est, lwr, upr, digits = 2) {
  paste0(fmt_num(est, digits), " [", fmt_num(lwr, digits), ", ", fmt_num(upr, digits), "]")
}
fmt_p <- function(x) {
  ifelse(x < 0.001, "<0.001", fmt_num(x, 3))
}

juvenile_mat <- read_count_matrix(juvenile_file)
female_mat <- read_count_matrix(female_file)
stopifnot(identical(rownames(juvenile_mat), rownames(female_mat)))
stopifnot(identical(colnames(juvenile_mat), colnames(female_mat)))

keep <- is.finite(juvenile_mat) & is.finite(female_mat)
ratio_all <- sum(juvenile_mat[keep]) / sum(female_mat[keep])

round_summary <- tibble(
  round = colnames(female_mat),
  juvenile = colSums(juvenile_mat, na.rm = TRUE),
  female = colSums(female_mat, na.rm = TRUE)
) %>%
  mutate(
    juvenile_female_ratio = juvenile / female,
    juveniles_per_100_females = 100 * juvenile_female_ratio
  )

camera_round <- expand_grid(
  camera_index = seq_len(nrow(female_mat)),
  round_index = seq_len(ncol(female_mat))
) %>%
  transmute(
    camera_id = rownames(female_mat)[camera_index],
    round = factor(colnames(female_mat)[round_index], levels = paste0("X", 1:7)),
    juvenile = juvenile_mat[cbind(camera_index, round_index)],
    female = female_mat[cbind(camera_index, round_index)]
  ) %>%
  filter(is.finite(juvenile), is.finite(female)) %>%
  mutate(
    juvenile_female_ratio = ifelse(female > 0, juvenile / female, NA_real_),
    juveniles_per_100_females = 100 * juvenile_female_ratio,
    total_identified = juvenile + female
  )

site_covs <- read_csv(site_cov_file, show_col_types = FALSE) %>%
  rename(camera_id = Camera.ID) %>%
  mutate(
    camera_id = as.character(camera_id),
    Polygon.Name = factor(Polygon.Name),
    Open.Spaces.z = z_scale(Open.Spaces),
    Nature.Reserves.z = z_scale(Nature.Reserves),
    Built.and.Roads.z = z_scale(Built.and.Roads),
    Null.Wolf.Abundance.z = z_scale(Null.Wolf.Abundance),
    null.jackal.abundance.z = z_scale(null.jackal.abundance),
    null.boar.abundance.z = z_scale(null.boar.abundance),
    human.abundance.z = z_scale(human.abundance)
  ) %>%
  select(
    camera_id,
    Polygon.Name,
    Open.Spaces.z,
    Nature.Reserves.z,
    Built.and.Roads.z,
    Null.Wolf.Abundance.z,
    null.jackal.abundance.z,
    null.boar.abundance.z,
    human.abundance.z
  )

camera_round <- camera_round %>%
  left_join(site_covs, by = "camera_id") %>%
  mutate(camera_id = factor(camera_id))

model_data <- camera_round %>%
  filter(total_identified > 0)

ratio_fit <- glmmTMB(
  cbind(juvenile, female) ~
    Open.Spaces.z +
    Nature.Reserves.z +
    Built.and.Roads.z +
    Null.Wolf.Abundance.z +
    null.jackal.abundance.z +
    null.boar.abundance.z +
    human.abundance.z +
    round +
    (1 | Polygon.Name / camera_id),
  family = betabinomial(link = "logit"),
  data = model_data
)

fit_sum <- summary(ratio_fit)
coef_tab <- as.data.frame(fit_sum$coefficients$cond)
coef_tab$term <- rownames(coef_tab)
rownames(coef_tab) <- NULL
coef_tab <- coef_tab %>%
  transmute(
    term,
    estimate = Estimate,
    std_error = `Std. Error`,
    z_value = `z value`,
    p_value = `Pr(>|z|)`,
    odds_ratio = exp(estimate),
    conf_low = exp(estimate - 1.96 * std_error),
    conf_high = exp(estimate + 1.96 * std_error)
  )

round_ref <- levels(model_data$round)[1]
round_effects <- coef_tab %>%
  filter(term == "(Intercept)" | grepl("^round", term)) %>%
  mutate(
    round = case_when(
      term == "(Intercept)" ~ round_ref,
      TRUE ~ sub("^round", "", term)
    )
  ) %>%
  arrange(match(round, levels(model_data$round)))

vc <- as.data.frame(VarCorr(ratio_fit)$cond)
dispersion_par <- sigma(ratio_fit)

round_pred <- tibble(
  round = factor(levels(model_data$round), levels = levels(model_data$round)),
  Open.Spaces.z = 0,
  Nature.Reserves.z = 0,
  Built.and.Roads.z = 0,
  Null.Wolf.Abundance.z = 0,
  null.jackal.abundance.z = 0,
  null.boar.abundance.z = 0,
  human.abundance.z = 0,
  Polygon.Name = model_data$Polygon.Name[1],
  camera_id = model_data$camera_id[1]
)
round_pred$pred_prop <- predict(ratio_fit, newdata = round_pred, type = "response", re.form = NA)
round_pred <- round_pred %>%
  mutate(
    pred_ratio = pred_prop / (1 - pred_prop),
    pred_per_100 = 100 * pred_ratio
  )

write_csv(round_summary, paste0(output_prefix, "_by_round.csv"))
write_csv(camera_round, paste0(output_prefix, "_camera_round.csv"))
write_csv(coef_tab, paste0(output_prefix, "_model_coefficients.csv"))
write_csv(round_pred %>% select(round, pred_prop, pred_ratio, pred_per_100),
          paste0(output_prefix, "_model_round_predictions.csv"))

plot_round <- ggplot(round_summary, aes(x = round, y = juveniles_per_100_females, group = 1)) +
  geom_col(fill = "#6c8e5e", width = 0.72) +
  geom_text(aes(label = fmt_num(juveniles_per_100_females, 1)), vjust = -0.35, size = 3.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "Survey round",
    y = "Juveniles per 100 adult females",
    title = "Camera-based juvenile:adult-female ratio by round"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(paste0(output_prefix, "_by_round.png"), plot_round, width = 8.2, height = 4.8, dpi = 300)
ggsave(paste0(output_prefix, "_by_round.pdf"), plot_round, width = 8.2, height = 4.8)

sig_terms <- coef_tab %>%
  filter(term != "(Intercept)", p_value < 0.05) %>%
  arrange(p_value)

md <- c(
  "# Juvenile:Adult-female Ratio Analysis",
  "",
  "## Descriptive ratios",
  "",
  paste0(
    "Across all camera-rounds with finite counts for both classes, the pooled juvenile:adult-female ratio was `",
    fmt_num(ratio_all, 3),
    "`, equivalent to `",
    fmt_num(100 * ratio_all, 1),
    "` juveniles per 100 adult females."
  ),
  "",
  "| Round | Juveniles | Adult females | Ratio | Juveniles per 100 females |",
  "|---|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(round_summary))) {
  md <- c(
    md,
    paste0(
      "| `", round_summary$round[i], "` | `", round_summary$juvenile[i],
      "` | `", round_summary$female[i], "` | `",
      fmt_num(round_summary$juvenile_female_ratio[i], 3), "` | `",
      fmt_num(round_summary$juveniles_per_100_females[i], 1), "` |"
    )
  )
}

md <- c(
  md,
  "",
  "## Beta-binomial GLMM",
  "",
  "Response: `cbind(juvenile, female)` at the camera x round level, using only rows with at least one identified juvenile or adult female.",
  "",
  "Fixed effects: `Open.Spaces`, `Nature.Reserves`, `Built.and.Roads`, `Null.Wolf.Abundance`, `null.jackal.abundance`, `null.boar.abundance`, `human.abundance`, and `round`.",
  "",
  "Random effects: `(1 | Polygon.Name / camera_id)`.",
  "",
  paste0("Beta-binomial dispersion parameter (`sigma`) = `", fmt_num(dispersion_par, 3), "`."),
  "",
  "| Term | Log-odds | Odds ratio [95% CI] | p |",
  "|---|---:|---:|---:|"
)

for (i in seq_len(nrow(coef_tab))) {
  md <- c(
    md,
    paste0(
      "| `", coef_tab$term[i], "` | `", fmt_num(coef_tab$estimate[i], 3),
      "` | `", fmt_ci(coef_tab$odds_ratio[i], coef_tab$conf_low[i], coef_tab$conf_high[i], 2),
      "` | `", fmt_p(coef_tab$p_value[i]), "` |"
    )
  )
}

md <- c(md, "", "## Short interpretation", "")
if (nrow(sig_terms) == 0) {
  md <- c(md, "No fixed effects were clearly associated with the juvenile:adult-female composition index at `p < 0.05`.")
} else {
  md <- c(md, "Terms with `p < 0.05` were:")
  for (i in seq_len(nrow(sig_terms))) {
    direction <- ifelse(sig_terms$estimate[i] > 0, "higher", "lower")
    md <- c(
      md,
      paste0(
        "- `", sig_terms$term[i], "` was associated with `", direction,
        "` juvenile representation among identified juveniles and adult females (odds ratio `",
        fmt_num(sig_terms$odds_ratio[i], 2), "`, 95% CI `",
        fmt_num(sig_terms$conf_low[i], 2), " to ", fmt_num(sig_terms$conf_high[i], 2),
        "`, `p = ", fmt_p(sig_terms$p_value[i]), "`)."
      )
    )
  }
}

writeLines(md, paste0(output_prefix, "_results.md"))

cat("Saved juvenile:female ratio outputs with prefix:", output_prefix, "\n")
