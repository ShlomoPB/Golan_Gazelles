args <- commandArgs(trailingOnly = TRUE)
arg_map <- list(
  fawn_data = "gazelle_fawn_spNMix_data.rds",
  male_data = "gazelle_male_spNMix_data.rds",
  female_data = "gazelle_female_spNMix_data.rds",
  output_file = "gazelle_multispecies_NMix_data.rds"
)

for (arg in args) {
  pieces <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(pieces) == 2) {
    key <- sub("^--", "", pieces[1])
    arg_map[[key]] <- pieces[2]
  }
}

fawn <- readRDS(arg_map$fawn_data)
male <- readRDS(arg_map$male_data)
female <- readRDS(arg_map$female_data)

same_matrix <- function(a, b) {
  identical(dim(a), dim(b)) && identical(dimnames(a), dimnames(b)) &&
    all((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b))
}

same_df <- function(a, b) {
  identical(names(a), names(b)) &&
    identical(rownames(a), rownames(b)) &&
    all(mapply(function(x, y) {
      if (is.numeric(x) && is.numeric(y)) {
        isTRUE(all.equal(x, y, tolerance = 1e-12, check.attributes = FALSE))
      } else {
        identical(x, y)
      }
    }, a, b))
}

if (!same_df(fawn$abund.covs, male$abund.covs) || !same_df(fawn$abund.covs, female$abund.covs)) {
  stop("Abundance covariates differ across class-specific data objects.")
}

for (nm in names(fawn$det.covs)) {
  if (!same_matrix(fawn$det.covs[[nm]], male$det.covs[[nm]]) ||
      !same_matrix(fawn$det.covs[[nm]], female$det.covs[[nm]])) {
    stop("Detection covariate ", nm, " differs across class-specific data objects.")
  }
}

if (!isTRUE(all.equal(fawn$coords, male$coords, tolerance = 1e-12, check.attributes = FALSE)) ||
    !isTRUE(all.equal(fawn$coords, female$coords, tolerance = 1e-12, check.attributes = FALSE))) {
  stop("Coordinates differ across class-specific data objects.")
}

y <- array(
  NA_integer_,
  dim = c(3, nrow(fawn$y), ncol(fawn$y)),
  dimnames = list(
    class = c("fawn", "male", "female"),
    site = rownames(fawn$y),
    replicate = colnames(fawn$y)
  )
)
y[1, , ] <- fawn$y
y[2, , ] <- male$y
y[3, , ] <- female$y

data_list <- list(
  y = y,
  abund.covs = fawn$abund.covs,
  det.covs = fawn$det.covs,
  coords = fawn$coords,
  metadata = list(
    classes = c("fawn", "male", "female"),
    source_files = c(
      fawn = arg_map$fawn_data,
      male = arg_map$male_data,
      female = arg_map$female_data
    )
  )
)

saveRDS(data_list, arg_map$output_file)

cat("Prepared multi-class N-mixture data:\n")
cat("  classes:", dim(y)[1], "\n")
cat("  sites:", dim(y)[2], "\n")
cat("  rounds:", dim(y)[3], "\n")
cat("  output:", arg_map$output_file, "\n")
