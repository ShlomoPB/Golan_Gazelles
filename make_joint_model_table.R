if (!requireNamespace("readr", quietly = TRUE)) {
  stop("Package 'readr' is required.")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required.")
}

library(readr)
library(dplyr)

fmt_num <- function(x, digits = 2) formatC(x, digits = digits, format = "f")
fmt_ci <- function(mean, lwr, upr, digits = 2) {
  paste0(fmt_num(mean, digits), " [", fmt_num(lwr, digits), ", ", fmt_num(upr, digits), "]")
}

beta_tbl <- read_csv("multiclass_lfMsNMix_class_threats_beta.csv", show_col_types = FALSE)
ppc_tbl <- read_csv("multiclass_lfMsNMix_class_threats_ppc.csv", show_col_types = FALSE)
cor_tbl <- read_csv("multiclass_lfMsNMix_class_threats_latent_mu_cor.csv", show_col_types = FALSE)
kappa_tbl <- read_csv("multiclass_lfMsNMix_class_threats_kappa.csv", show_col_types = FALSE)

final_beta <- beta_tbl %>%
  filter(model == "final", term %in% c("Open.Spaces.z", "Nature.Reserves.z", "null.jackal.abundance.z", "Null.Wolf.Abundance.z")) %>%
  mutate(effect = fmt_ci(mean, lwr, upr)) %>%
  select(class, term, effect, rhat, ESS)

final_kappa <- kappa_tbl %>%
  filter(model == "final") %>%
  mutate(kappa = fmt_ci(mean, lwr, upr)) %>%
  select(class, kappa, rhat, ESS)

final_ppc <- ppc_tbl %>%
  filter(model == "final") %>%
  mutate(
    fit_y = fmt_num(fit_y),
    fit_y_rep = fmt_num(fit_y_rep),
    bayes_p = fmt_num(bayes_p, 3)
  )

final_cor <- cor_tbl %>%
  filter(model == "final") %>%
  mutate(correlation = fmt_ci(mean, lwr, upr)) %>%
  select(pair, correlation)

write_csv(final_beta, "publication_joint_model_coefficients.csv")
write_csv(final_kappa, "publication_joint_model_kappa.csv")
write_csv(final_ppc, "publication_joint_model_ppc.csv")
write_csv(final_cor, "publication_joint_model_correlations.csv")

lines <- c(
  "# Joint-model Summary Table",
  "",
  "Exploratory latent-factor negative-binomial N-mixture model with shared detection",
  "`Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera + log_days_operational.z`.",
  "",
  "## Species-level coefficients",
  "",
  "| Class | Term | Posterior mean [95% CrI] | Rhat | ESS |",
  "|---|---|---:|---:|---:|"
)

for (i in seq_len(nrow(final_beta))) {
  lines <- c(
    lines,
    paste0(
      "| `", final_beta$class[i], "` | `", final_beta$term[i], "` | `", final_beta$effect[i],
      "` | `", fmt_num(final_beta$rhat[i]), "` | `", round(final_beta$ESS[i]), "` |"
    )
  )
}

lines <- c(
  lines,
  "",
  "## Negative-binomial overdispersion",
  "",
  "| Class | kappa | Rhat | ESS |",
  "|---|---:|---:|---:|"
)

for (i in seq_len(nrow(final_kappa))) {
  lines <- c(
    lines,
    paste0(
      "| `", final_kappa$class[i], "` | `", final_kappa$kappa[i], "` | `",
      fmt_num(final_kappa$rhat[i]), "` | `", round(final_kappa$ESS[i]), "` |"
    )
  )
}

lines <- c(
  lines,
  "",
  "## Posterior predictive check",
  "",
  "| Class | Fit(y) | Fit(yrep) | Bayesian p |",
  "|---|---:|---:|---:|"
)

for (i in seq_len(nrow(final_ppc))) {
  lines <- c(
    lines,
    paste0(
      "| `", final_ppc$class[i], "` | `", final_ppc$fit_y[i], "` | `",
      final_ppc$fit_y_rep[i], "` | `", final_ppc$bayes_p[i], "` |"
    )
  )
}

lines <- c(
  lines,
  "",
  "## Residual cross-class correlations",
  "",
  "| Pair | Posterior mean [95% CrI] |",
  "|---|---:|"
)

for (i in seq_len(nrow(final_cor))) {
  lines <- c(
    lines,
    paste0("| `", final_cor$pair[i], "` | `", final_cor$correlation[i], "` |")
  )
}

writeLines(lines, "publication_joint_model_table.md")
cat("Saved joint-model table outputs\n")
