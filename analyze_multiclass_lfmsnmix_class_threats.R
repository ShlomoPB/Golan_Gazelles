args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(x) {
  out <- list(
    baseline = "gazelle_lfMsNMix_2f_NB_shared_env_1200b_25bl_5th_3ch_fit.rds",
    design = "gazelle_lfMsNMix_2f_NB_class_threats_det_design_2000b_25bl_5th_3ch_fit.rds",
    design_effort = "gazelle_lfMsNMix_2f_NB_class_threats_det_design_effort_2000b_25bl_5th_3ch_fit.rds",
    final = "gazelle_lfMsNMix_2f_NB_class_threats_det_design_effort_5000b_25bl_5th_3ch_fit.rds",
    stride = "100",
    prefix = "multiclass_lfMsNMix_class_threats"
  )
  if (length(x) == 0) {
    return(out)
  }
  for (arg in x) {
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(parts) == 2) {
      out[[parts[1]]] <- parts[2]
    }
  }
  out
}

cfg <- parse_args(args)

fmt_num <- function(x, digits = 2) formatC(x, digits = digits, format = "f")
fmt_ci <- function(mean, lwr, upr, digits = 2) {
  paste0(fmt_num(mean, digits), " [", fmt_num(lwr, digits), ", ", fmt_num(upr, digits), "]")
}

sum_draws <- function(x) {
  c(
    mean = mean(x),
    lwr = unname(stats::quantile(x, 0.025)),
    upr = unname(stats::quantile(x, 0.975))
  )
}

extract_species_params <- function(fit, samples, names_vec, rhat_mat, ess_mat) {
  samples <- as.matrix(samples)
  n_sp <- length(fit$sp.names)
  n_coef <- length(names_vec)
  out <- vector("list", n_sp * n_coef)
  idx <- 1
  for (k in seq_len(n_coef)) {
    cols <- ((k - 1) * n_sp + 1):(k * n_sp)
    for (s in seq_len(n_sp)) {
      sm <- sum_draws(samples[, cols[s]])
      out[[idx]] <- data.frame(
        class = fit$sp.names[s],
        term = names_vec[k],
        mean = sm["mean"],
        lwr = sm["lwr"],
        upr = sm["upr"],
        rhat = as.vector(rhat_mat)[cols[s]],
        ESS = as.vector(ess_mat)[cols[s]]
      )
      idx <- idx + 1
    }
  }
  do.call(rbind, out)
}

extract_kappa <- function(fit) {
  sm <- as.data.frame(t(apply(as.matrix(fit$kappa.samples), 2, sum_draws)))
  sm$class <- fit$sp.names
  sm$rhat <- fit$rhat$kappa
  sm$ESS <- fit$ESS$kappa
  sm[, c("class", "mean", "lwr", "upr", "rhat", "ESS")]
}

latent_mu_cor <- function(fit) {
  vals <- t(apply(fit$mu.samples, 1, function(v) {
    m <- matrix(v, nrow = length(fit$sp.names), byrow = TRUE)
    c(
      "fawn-male" = stats::cor(log(m[1, ]), log(m[2, ])),
      "fawn-female" = stats::cor(log(m[1, ]), log(m[3, ])),
      "male-female" = stats::cor(log(m[2, ]), log(m[3, ]))
    )
  }))
  out <- as.data.frame(t(apply(vals, 2, sum_draws)))
  out$pair <- rownames(out)
  out[, c("pair", "mean", "lwr", "upr")]
}

manual_ppc <- function(fit, stride = 100L) {
  set.seed(20260615)
  alpha_samples <- as.matrix(fit$alpha.samples)
  keep <- seq(1, nrow(alpha_samples), by = stride)
  sp_names <- fit$sp.names
  x_p_names <- fit$x.p.names
  y <- fit$y
  obs_mask <- !is.na(y[1, , ])
  y_grouped <- apply(y, c(1, 2), sum, na.rm = TRUE)
  n_sites <- dim(y)[2]
  n_reps <- dim(y)[3]
  n_coef <- length(x_p_names)
  fit_y <- matrix(NA_real_, length(keep), length(sp_names))
  fit_y_rep <- matrix(NA_real_, length(keep), length(sp_names))

  for (ii in seq_along(keep)) {
    samp <- keep[ii]
    for (s in seq_along(sp_names)) {
      cols <- ((s - 1) * n_coef + 1):(s * n_coef)
      alpha_vec <- alpha_samples[samp, cols]
      eta <- as.numeric(fit$X.p %*% alpha_vec)
      p_obs <- plogis(eta)
      p_mat <- matrix(NA_real_, n_sites, n_reps)
      p_mat[obs_mask] <- p_obs
      mu_vec <- fit$mu.samples[samp, s, ]
      E_grouped <- rowSums(sweep(p_mat, 1, mu_vec, `*`), na.rm = TRUE)
      N_rep <- stats::rnbinom(n_sites, mu = mu_vec, size = fit$kappa.samples[samp, s])
      y_rep_grouped <- numeric(n_sites)
      for (j in seq_len(n_sites)) {
        reps <- which(obs_mask[j, ])
        if (length(reps) > 0) {
          y_rep_grouped[j] <- sum(stats::rbinom(length(reps), size = N_rep[j], prob = p_mat[j, reps]))
        }
      }
      fit_y[ii, s] <- sum((sqrt(y_grouped[s, ]) - sqrt(E_grouped))^2)
      fit_y_rep[ii, s] <- sum((sqrt(y_rep_grouped) - sqrt(E_grouped))^2)
    }
  }

  data.frame(
    class = sp_names,
    fit_y = colMeans(fit_y),
    fit_y_rep = colMeans(fit_y_rep),
    bayes_p = colMeans(fit_y_rep > fit_y)
  )
}

extract_model <- function(file, model_name, stride) {
  fit <- readRDS(file)
  list(
    model = model_name,
    fit = fit,
    beta = extract_species_params(fit, fit$beta.samples, fit$x.names, fit$rhat$beta, fit$ESS$beta),
    alpha = extract_species_params(fit, fit$alpha.samples, fit$x.p.names, fit$rhat$alpha, fit$ESS$alpha),
    kappa = extract_kappa(fit),
    cor = latent_mu_cor(fit),
    ppc = manual_ppc(fit, stride = stride)
  )
}

models <- list(
  baseline = cfg$baseline,
  design = cfg$design,
  design_effort = cfg$design_effort,
  final = cfg$final
)
stride <- as.integer(cfg$stride)
res <- lapply(names(models), function(nm) extract_model(models[[nm]], nm, stride))
names(res) <- names(models)

beta_tbl <- do.call(rbind, lapply(res, function(x) transform(x$beta, model = x$model)))
alpha_tbl <- do.call(rbind, lapply(res, function(x) transform(x$alpha, model = x$model)))
kappa_tbl <- do.call(rbind, lapply(res, function(x) transform(x$kappa, model = x$model)))
cor_tbl <- do.call(rbind, lapply(res, function(x) transform(x$cor, model = x$model)))
ppc_tbl <- do.call(rbind, lapply(res, function(x) transform(x$ppc, model = x$model)))

utils::write.csv(beta_tbl, paste0(cfg$prefix, "_beta.csv"), row.names = FALSE)
utils::write.csv(alpha_tbl, paste0(cfg$prefix, "_alpha.csv"), row.names = FALSE)
utils::write.csv(kappa_tbl, paste0(cfg$prefix, "_kappa.csv"), row.names = FALSE)
utils::write.csv(cor_tbl, paste0(cfg$prefix, "_latent_mu_cor.csv"), row.names = FALSE)
utils::write.csv(ppc_tbl, paste0(cfg$prefix, "_ppc.csv"), row.names = FALSE)

final_beta <- subset(beta_tbl, model == "final" & term %in% c(
  "Nature.Reserves.z", "Open.Spaces.z", "null.boar.abundance.z",
  "null.jackal.abundance.z", "Null.Wolf.Abundance.z"
))
final_kappa <- subset(kappa_tbl, model == "final")
final_ppc <- subset(ppc_tbl, model == "final")
cor_compare <- subset(cor_tbl, model %in% c("baseline", "design_effort", "final"))
kappa_compare <- subset(kappa_tbl, model %in% c("baseline", "design_effort", "final"))

lines <- c(
  "# Multi-class lfMsNMix class-threat model",
  "",
  "## Model",
  "",
  "- Model type: `spAbundance::lfMsNMix`",
  "- Family: negative binomial",
  "- Latent factors: `2`",
  "- Shared abundance formula used in the fitted model:",
  "  `~ Nature.Reserves.z + Open.Spaces.z + null.boar.abundance.z + null.jackal.abundance.z + Null.Wolf.Abundance.z`",
  "- Shared detection formula in the selected final run:",
  "  `~ Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera + log_days_operational.z`",
  "",
  "Note: `lfMsNMix` does not support truly different abundance formulas by class. The requested class-specific threat model was approximated by fitting the union of the requested covariates and then interpreting the class-specific posterior coefficients.",
  "",
  "## Detection Comparison",
  "",
  "The 2000-batch comparison favored keeping `log_days_operational.z` in the detection model. Relative to the design-only detection model, the `design_effort` fit improved community-level abundance mixing and female `kappa` mixing, while leaving the main class-specific slope directions unchanged.",
  "",
  "## Final Class-specific Threat Coefficients",
  "",
  "| Class | Term | Posterior mean | 95% CrI | Rhat | ESS |",
  "|---|---|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(final_beta))) {
  lines <- c(
    lines,
    paste0(
      "| `", final_beta$class[i], "` | `", final_beta$term[i], "` | `",
      fmt_num(final_beta$mean[i]), "` | `", fmt_ci(final_beta$mean[i], final_beta$lwr[i], final_beta$upr[i]),
      "` | `", fmt_num(final_beta$rhat[i]), "` | `", round(final_beta$ESS[i]), "` |"
    )
  )
}

lines <- c(
  lines,
  "",
  "## Kappa",
  "",
  "| Model | Class | Mean | 95% CrI | Rhat | ESS |",
  "|---|---|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(kappa_compare))) {
  lines <- c(
    lines,
    paste0(
      "| `", kappa_compare$model[i], "` | `", kappa_compare$class[i], "` | `",
      fmt_num(kappa_compare$mean[i]), "` | `", fmt_ci(kappa_compare$mean[i], kappa_compare$lwr[i], kappa_compare$upr[i]),
      "` | `", fmt_num(kappa_compare$rhat[i]), "` | `", round(kappa_compare$ESS[i]), "` |"
    )
  )
}

lines <- c(
  lines,
  "",
  "## Posterior Predictive Check",
  "",
  "Manual site-grouped Freeman-Tukey PPC using every 100th posterior draw.",
  "",
  "| Class | Fit(y) | Fit(yrep) | Bayesian p |",
  "|---|---:|---:|---:|"
)

for (i in seq_len(nrow(final_ppc))) {
  lines <- c(
    lines,
    paste0(
      "| `", final_ppc$class[i], "` | `", fmt_num(final_ppc$fit_y[i]), "` | `",
      fmt_num(final_ppc$fit_y_rep[i]), "` | `", fmt_num(final_ppc$bayes_p[i], 3), "` |"
    )
  )
}

lines <- c(
  lines,
  "",
  "## Latent Mean-abundance Correlation",
  "",
  "These are posterior correlations among site-level log latent mean abundances (`log(mu)`) for each class pair. This is a directly comparable cross-class association metric across model variants.",
  "",
  "| Model | Pair | Mean | 95% CrI |",
  "|---|---|---:|---:|"
)

for (i in seq_len(nrow(cor_compare))) {
  lines <- c(
    lines,
    paste0(
      "| `", cor_compare$model[i], "` | `", cor_compare$pair[i], "` | `",
      fmt_num(cor_compare$mean[i]), "` | `", fmt_ci(cor_compare$mean[i], cor_compare$lwr[i], cor_compare$upr[i]), "` |"
    )
  )
}

writeLines(lines, paste0(cfg$prefix, "_results.md"))
