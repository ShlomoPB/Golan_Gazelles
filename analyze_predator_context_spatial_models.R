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
male_file <- "gazelle_male_abundance_60x7.csv"
category_file <- "publication_predator_context_q33_67_camera_categories.csv"
output_prefix <- "publication_predator_context_q33_67_spatial_model"

context_levels <- c(
  "Intermediate / mixed",
  "Low both",
  "High both",
  "Low jackal / High wolf",
  "High jackal / Low wolf"
)

class_levels <- c("Female", "Male", "Fawn")

read_count_matrix <- function(path) {
  x <- read.csv(path, na.strings = c("NA", ""), check.names = FALSE)
  y <- as.matrix(x[paste0("X", 1:7)])
  storage.mode(y) <- "numeric"
  rownames(y) <- x$camera_name
  colnames(y) <- paste0("X", 1:7)
  y
}

fmt_num <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
fmt_ci <- function(est, lwr, upr, digits = 2) {
  paste0(fmt_num(est, digits), " [", fmt_num(lwr, digits), ", ", fmt_num(upr, digits), "]")
}
fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", fmt_num(x, 3)))
}

juvenile_mat <- read_count_matrix(juvenile_file)
female_mat <- read_count_matrix(female_file)
male_mat <- read_count_matrix(male_file)
stopifnot(identical(rownames(juvenile_mat), rownames(female_mat)))
stopifnot(identical(rownames(juvenile_mat), rownames(male_mat)))

camera_ids <- rownames(juvenile_mat)
active_rounds <- rowSums(is.finite(female_mat))

camera_counts <- tibble(
  camera_name = as.integer(camera_ids),
  active_rounds = active_rounds,
  Fawn = rowSums(juvenile_mat, na.rm = TRUE),
  Female = rowSums(female_mat, na.rm = TRUE),
  Male = rowSums(male_mat, na.rm = TRUE)
) %>%
  pivot_longer(cols = c(Fawn, Female, Male), names_to = "class", values_to = "count")

category_df <- read_csv(category_file, show_col_types = FALSE) %>%
  select(camera_name, Polygon.Name, predator_context_label)

dat <- camera_counts %>%
  left_join(category_df, by = "camera_name") %>%
  mutate(
    class = factor(class, levels = class_levels),
    predator_context_label = factor(predator_context_label, levels = context_levels),
    camera_name = factor(camera_name),
    Polygon.Name = factor(Polygon.Name)
  )

additive_fit <- glmmTMB(
  count ~ class + predator_context_label + offset(log(active_rounds)) + (1 | Polygon.Name / camera_name),
  family = nbinom2(link = "log"),
  data = dat
)

interaction_fit <- glmmTMB(
  count ~ class * predator_context_label + offset(log(active_rounds)) + (1 | Polygon.Name / camera_name),
  family = nbinom2(link = "log"),
  data = dat
)

lrt <- anova(additive_fit, interaction_fit)

coef_tab <- summary(interaction_fit)$coefficients$cond
coef_df <- as.data.frame(coef_tab)
coef_df$term <- rownames(coef_tab)
rownames(coef_df) <- NULL

newdat <- expand.grid(
  class = factor(class_levels, levels = class_levels),
  predator_context_label = factor(context_levels, levels = context_levels),
  active_rounds = 1
)

X <- model.matrix(~ class * predator_context_label, data = newdat)
beta <- fixef(interaction_fit)$cond
V <- vcov(interaction_fit)$cond
eta <- as.vector(X %*% beta)
se_eta <- sqrt(diag(X %*% V %*% t(X)))

pred_df <- newdat %>%
  mutate(
    eta = eta,
    se_eta = se_eta,
    pred_count_per_round = exp(eta),
    lwr = exp(eta - 1.96 * se_eta),
    upr = exp(eta + 1.96 * se_eta)
  )

make_contrast <- function(class_name, context_a, context_b = "Intermediate / mixed") {
  nd <- data.frame(
    class = factor(c(class_name, class_name), levels = class_levels),
    predator_context_label = factor(c(context_a, context_b), levels = context_levels)
  )
  Xc <- model.matrix(~ class * predator_context_label, data = nd)
  cvec <- Xc[1, ] - Xc[2, ]
  est <- sum(cvec * beta)
  se <- sqrt(as.numeric(t(cvec) %*% V %*% cvec))
  tibble(
    class = class_name,
    contrast = paste0(context_a, " vs ", context_b),
    log_rate_ratio = est,
    rate_ratio = exp(est),
    lwr = exp(est - 1.96 * se),
    upr = exp(est + 1.96 * se),
    z = est / se,
    p_value = 2 * pnorm(abs(est / se), lower.tail = FALSE)
  )
}

contrast_df <- bind_rows(
  lapply(class_levels, function(cls) make_contrast(cls, "High jackal / Low wolf")),
  lapply(class_levels, function(cls) make_contrast(cls, "Low jackal / High wolf")),
  lapply(class_levels, function(cls) make_contrast(cls, "High both")),
  lapply(class_levels, function(cls) make_contrast(cls, "Low both")),
  lapply(class_levels, function(cls) make_contrast(cls, "High jackal / Low wolf", "Low jackal / High wolf"))
)

write_csv(coef_df, paste0(output_prefix, "_coefficients.csv"))
write_csv(pred_df, paste0(output_prefix, "_predicted_counts.csv"))
write_csv(contrast_df, paste0(output_prefix, "_contrasts.csv"))

heatmap_plot <- ggplot(pred_df, aes(x = predator_context_label, y = class, fill = pred_count_per_round)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = fmt_num(pred_count_per_round, 2)), size = 3.2) +
  scale_fill_gradient(low = "#dceaf5", high = "#1d3557") +
  labs(
    x = NULL,
    y = NULL,
    fill = "Predicted count\nper active round",
    title = "Predicted class-specific camera-level counts across predator contexts"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 18, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(paste0(output_prefix, "_heatmap.png"), heatmap_plot, width = 9.2, height = 4.8, dpi = 300)
ggsave(paste0(output_prefix, "_heatmap.pdf"), heatmap_plot, width = 9.2, height = 4.8)

md <- c(
  "# Spatial predator-context count model",
  "",
  "Camera-level class counts were modeled with a negative-binomial GLMM using an offset for the number of active rounds.",
  "",
  "Model formula:",
  "",
  "`count ~ class * predator_context_label + offset(log(active_rounds)) + (1 | Polygon.Name / camera_name)`",
  "",
  "## Omnibus interaction test",
  "",
  paste0(
    "Likelihood-ratio test comparing additive and interaction models: `Chisq = ",
    fmt_num(lrt$Chisq[2], 3), "`, `df = ", lrt$Df[2] - lrt$Df[1], "`, `p = ", fmt_p(lrt$`Pr(>Chisq)`[2]), "`."
  ),
  "",
  "## Predicted count per active round",
  "",
  "| Class | Predator context | Predicted count [95% CI] |",
  "|---|---|---:|"
)

for (i in seq_len(nrow(pred_df))) {
  md <- c(
    md,
    paste0(
      "| `", pred_df$class[i], "` | `", pred_df$predator_context_label[i], "` | `",
      fmt_ci(pred_df$pred_count_per_round[i], pred_df$lwr[i], pred_df$upr[i], 2), "` |"
    )
  )
}

md <- c(
  md,
  "",
  "## Context contrasts",
  "",
  "| Class | Contrast | Rate ratio [95% CI] | p |",
  "|---|---|---:|---:|"
)

for (i in seq_len(nrow(contrast_df))) {
  md <- c(
    md,
    paste0(
      "| `", contrast_df$class[i], "` | `", contrast_df$contrast[i], "` | `",
      fmt_ci(contrast_df$rate_ratio[i], contrast_df$lwr[i], contrast_df$upr[i], 2),
      "` | `", fmt_p(contrast_df$p_value[i]), "` |"
    )
  )
}

writeLines(md, paste0(output_prefix, "_results.md"))

cat("Saved predator-context spatial model outputs with prefix:", output_prefix, "\n")
