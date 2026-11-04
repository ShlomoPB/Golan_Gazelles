library(unmarked)
library(AICcmodavg)
setwd("C:/Users/preissbloom/Documents/Github Repositories/Golan_Gazelles/Gazelle VSCode")
# ============================================================================
# LOAD DATA
# ============================================================================

y_male <- read.csv("gazelle_male_abundance_60x7 - Copy.csv")
y_female <- read.csv("gazelle_female_abundance_60x7 - Copy.csv")
y_juvenile <- read.csv("gazelle_juvenile_abundance_60x7 - Copy.csv")

siteCovs <- read.csv("SiteCovs_060324.csv")
obsCovs <- read.csv("ObsCovs_050524.csv")

# Convert character columns to factors
siteCovs[sapply(siteCovs, is.character)] <- 
  lapply(siteCovs[sapply(siteCovs, is.character)], as.factor)
obsCovs[sapply(obsCovs, is.character)] <- 
  lapply(obsCovs[sapply(obsCovs, is.character)], as.factor)

# Create unmarked frames for each class
umf_male <- unmarkedFramePCount(y = y_male, siteCovs = siteCovs, obsCovs = obsCovs)
umf_female <- unmarkedFramePCount(y = y_female, siteCovs = siteCovs, obsCovs = obsCovs)
umf_juvenile <- unmarkedFramePCount(y = y_juvenile, siteCovs = siteCovs, obsCovs = obsCovs)

# ============================================================================
# DETECTION MODELS - MALE
# ============================================================================

male_det_models <- list(
  g20 = pcount(~1 ~ 1, umf_male),
  g21 = pcount(~Cattle.Sightings ~ 1, umf_male),
  g22 = pcount(~Camera.Model ~ 1, umf_male),
  g23 = pcount(~Camera.Sensitivity ~ 1, umf_male),
  g24 = pcount(~Livestock.Present. ~ 1, umf_male),
  g25 = pcount(~Capture.Delay ~ 1, umf_male),
  g26 = pcount(~Vegetation.Structure ~ 1, umf_male),
  g27 = pcount(~Path.Type ~ 1, umf_male),
  g28 = pcount(~Area.Behind.Camera ~ 1, umf_male),
  g29 = pcount(~Bottleneck. ~ 1, umf_male),
  g30 = pcount(~Agricultural.Presence ~ 1, umf_male),
  g3 = pcount(~Path.Type + Camera.Model ~ 1, umf_male),
  g5 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera ~ 1, umf_male),
  g8 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ 1, umf_male),
  g10 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 1, umf_male),
  g11 = pcount(~Path.Type + Camera.Model + Vegetation.Structure ~ 1, umf_male)
)

# ============================================================================
# DETECTION MODELS - FEMALE
# ============================================================================

female_det_models <- list(
  g20 = pcount(~1 ~ 1, umf_female),
  g21 = pcount(~Cattle.Sightings ~ 1, umf_female),
  g22 = pcount(~Camera.Model ~ 1, umf_female),
  g23 = pcount(~Camera.Sensitivity ~ 1, umf_female),
  g24 = pcount(~Livestock.Present. ~ 1, umf_female),
  g25 = pcount(~Capture.Delay ~ 1, umf_female),
  g26 = pcount(~Vegetation.Structure ~ 1, umf_female),
  g27 = pcount(~Path.Type ~ 1, umf_female),
  g28 = pcount(~Area.Behind.Camera ~ 1, umf_female),
  g29 = pcount(~Bottleneck. ~ 1, umf_female),
  g30 = pcount(~Agricultural.Presence ~ 1, umf_female),
  g3 = pcount(~Path.Type + Camera.Model ~ 1, umf_female),
  g5 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera ~ 1, umf_female),
  g8 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ 1, umf_female),
  g10 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 1, umf_female),
  g11 = pcount(~Path.Type + Camera.Model + Vegetation.Structure ~ 1, umf_female)
)

# ============================================================================
# DETECTION MODELS - JUVENILE
# ============================================================================

juvenile_det_models <- list(
  g20 = pcount(~1 ~ 1, umf_juvenile),
  g21 = pcount(~Cattle.Sightings ~ 1, umf_juvenile),
  g22 = pcount(~Camera.Model ~ 1, umf_juvenile),
  g23 = pcount(~Camera.Sensitivity ~ 1, umf_juvenile),
  g24 = pcount(~Livestock.Present. ~ 1, umf_juvenile),
  g25 = pcount(~Capture.Delay ~ 1, umf_juvenile),
  g26 = pcount(~Vegetation.Structure ~ 1, umf_juvenile),
  g27 = pcount(~Path.Type ~ 1, umf_juvenile),
  g28 = pcount(~Area.Behind.Camera ~ 1, umf_juvenile),
  g29 = pcount(~Bottleneck. ~ 1, umf_juvenile),
  g30 = pcount(~Agricultural.Presence ~ 1, umf_juvenile),
  g3 = pcount(~Path.Type + Camera.Model ~ 1, umf_juvenile),
  g5 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera ~ 1, umf_juvenile),
  g8 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ 1, umf_juvenile),
  g10 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 1, umf_juvenile),
  g11 = pcount(~Path.Type + Camera.Model + Vegetation.Structure ~ 1, umf_juvenile)
)

# ============================================================================
# MODEL SELECTION - DETECTION MODELS
# ============================================================================

male_det_sel <- modSel(fitList(fits = male_det_models), nullmod = "g20")
female_det_sel <- modSel(fitList(fits = female_det_models), nullmod = "g20")
juvenile_det_sel <- modSel(fitList(fits = juvenile_det_models), nullmod = "g20")

cat("\n===== MALE DETECTION MODELS =====\n")
print(male_det_sel)
cat("\n===== FEMALE DETECTION MODELS =====\n")
print(female_det_sel)
cat("\n===== JUVENILE DETECTION MODELS =====\n")
print(juvenile_det_sel)

write.csv(as(male_det_sel, "data.frame"), 
          "Gazelle_Detection_Models_MALE.csv", row.names = FALSE)
write.csv(as(female_det_sel, "data.frame"), 
          "Gazelle_Detection_Models_FEMALE.csv", row.names = FALSE)
write.csv(as(juvenile_det_sel, "data.frame"), 
          "Gazelle_Detection_Models_JUVENILE.csv", row.names = FALSE)

# ============================================================================
# ABUNDANCE MODELS - MALE
# ============================================================================

male_abund_models <- list(
  gm0 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ 1, umf_male),
  
  gm2 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ Null.Wolf.Abundance, umf_male),
  
  gm3 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ Built.and.Roads, umf_male),
  
  gm4 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ Open.Spaces, umf_male),
  
  gm5 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ null.jackal.abundance, umf_male),
  
  gm6 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ Nature.Reserves, umf_male),
  
  gm7 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ Firing.Zones, umf_male),
  
  gm8 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ Orchards, umf_male),
  
  gm9 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ Cropland, umf_male),
  
  gm10 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ Cattle.Abundance, umf_male),
  
  gm11 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ Null.Wolf.Abundance * Nature.Reserves + null.jackal.abundance, umf_male),
  
  gm12 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ Null.Wolf.Abundance + Nature.Reserves * null.jackal.abundance, umf_male),
  
  gm13 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ Null.Wolf.Abundance + null.boar.abundance + null.jackal.abundance, umf_male),
  
  gm14 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ null.boar.abundance * Nature.Reserves, umf_male),
  
  gm15 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ Null.Wolf.Abundance * Open.Spaces, umf_male),
  
  gm16 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ null.boar.abundance * Open.Spaces, umf_male),
  
  gm17 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ 
                  Null.Wolf.Abundance * null.boar.abundance, umf_male),
  
  gm18 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ 
                  Null.Wolf.Abundance * null.jackal.abundance, umf_male),
  
  gm19 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ 
                  Null.Wolf.Abundance * Cattle.Abundance, umf_male),
  
  gm20 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ 
                  Null.Wolf.Abundance * Built.and.Roads, umf_male),
  
  gm21 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ 
                  Cropland * Open.Spaces, umf_male),
  
  gm22 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ 
                  Null.Wolf.Abundance * Nature.Reserves * null.jackal.abundance, umf_male)
)

# ============================================================================
# ABUNDANCE MODELS - FEMALE
# ============================================================================

female_abund_models <- list(
  gm0 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 1, umf_female),
  
  gm2 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Null.Wolf.Abundance, umf_female),
  
  gm3 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Built.and.Roads, umf_female),
  
  gm4 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Open.Spaces, umf_female),
  
  gm5 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ null.jackal.abundance, umf_female),
  
  gm6 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Nature.Reserves, umf_female),
  
  gm7 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Firing.Zones, umf_female),
  
  gm8 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Orchards, umf_female),
  
  gm9 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Cropland, umf_female),
  
  gm10 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Cattle.Abundance, umf_female),
  
  gm11 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Null.Wolf.Abundance * Nature.Reserves + null.jackal.abundance, umf_female),
  
  gm12 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Null.Wolf.Abundance + Nature.Reserves * null.jackal.abundance, umf_female),
  
  gm13 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Null.Wolf.Abundance + null.boar.abundance + null.jackal.abundance, umf_female),
  
  gm14 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ null.boar.abundance * Nature.Reserves, umf_female),
  
  gm15 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Null.Wolf.Abundance * Open.Spaces, umf_female),
  
  gm16 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ null.boar.abundance * Open.Spaces, umf_female),
  
  gm17 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 
                  Null.Wolf.Abundance * null.boar.abundance, umf_female),
  
  gm18 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 
                  Null.Wolf.Abundance * null.jackal.abundance, umf_female),
  
  gm19 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 
                  Null.Wolf.Abundance * Cattle.Abundance, umf_female),
  
  gm20 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 
                  Null.Wolf.Abundance * Built.and.Roads, umf_female),
  
  gm21 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 
                  Cropland * Open.Spaces, umf_female),
  
  gm22 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 
                  Null.Wolf.Abundance * Nature.Reserves * null.jackal.abundance, umf_female)
)

# ============================================================================
# ABUNDANCE MODELS - JUVENILE
# ============================================================================

juvenile_abund_models <- list(
  gm0 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 1, umf_juvenile),
  
  gm2 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Null.Wolf.Abundance, umf_juvenile),
  
  gm3 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Built.and.Roads, umf_juvenile),
  
  gm4 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Open.Spaces, umf_juvenile),
  
  gm5 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ null.jackal.abundance, umf_juvenile),
  
  gm6 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Nature.Reserves, umf_juvenile),
  
  gm7 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Firing.Zones, umf_juvenile),
  
  gm8 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Orchards, umf_juvenile),
  
  gm9 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Cropland, umf_juvenile),
  
  gm10 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Cattle.Abundance, umf_juvenile),
  
  gm11 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Null.Wolf.Abundance * Nature.Reserves + null.jackal.abundance, umf_juvenile),
  
  gm12 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Null.Wolf.Abundance + Nature.Reserves * null.jackal.abundance, umf_juvenile),
  
  gm13 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Null.Wolf.Abundance + null.boar.abundance + null.jackal.abundance, umf_juvenile),
  
  gm14 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ null.boar.abundance * Nature.Reserves, umf_juvenile),

  gm15 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Null.Wolf.Abundance * Open.Spaces, umf_juvenile),
  
  gm16 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ null.boar.abundance * Open.Spaces, umf_juvenile),
  
  gm17 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 
                  Null.Wolf.Abundance * null.boar.abundance, umf_juvenile),
  
  gm18 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 
                  Null.Wolf.Abundance * null.jackal.abundance, umf_juvenile),
  
  gm19 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 
                  Null.Wolf.Abundance * Cattle.Abundance, umf_juvenile),
  
  gm20 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 
                  Null.Wolf.Abundance * Built.and.Roads, umf_juvenile),
  
  gm21 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 
                  Cropland * Open.Spaces, umf_juvenile),
  
  gm22 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 
                  Null.Wolf.Abundance * Nature.Reserves * null.jackal.abundance, umf_juvenile)
  
  )

# ============================================================================
# MODEL SELECTION - ABUNDANCE MODELS
# ============================================================================

male_abund_sel <- modSel(fitList(fits = male_abund_models), nullmod = "gm0")
female_abund_sel <- modSel(fitList(fits = female_abund_models), nullmod = "gm0")
juvenile_abund_sel <- modSel(fitList(fits = juvenile_abund_models), nullmod = "gm0")

cat("\n===== MALE ABUNDANCE MODELS =====\n")
print(male_abund_sel)
cat("\n===== FEMALE ABUNDANCE MODELS =====\n")
print(female_abund_sel)
cat("\n===== JUVENILE ABUNDANCE MODELS =====\n")
print(juvenile_abund_sel)

write.csv(as(male_abund_sel, "data.frame"), 
          "Gazelle_Abundance_Models_MALE.csv", row.names = FALSE)
write.csv(as(female_abund_sel, "data.frame"), 
          "Gazelle_Abundance_Models_FEMALE.csv", row.names = FALSE)
write.csv(as(juvenile_abund_sel, "data.frame"), 
          "Gazelle_Abundance_Models_JUVENILE.csv", row.names = FALSE)

# ============================================================================
# GOODNESS OF FIT (adjust model names based on your top models)
# ============================================================================

best_male <- male_abund_models$gm13
best_female <- female_abund_models$gm22
best_juvenile <- juvenile_abund_models$gm4

male_gof <- Nmix.gof.test(best_male, nsim = 100)
female_gof <- Nmix.gof.test(best_female, nsim = 100)
juvenile_gof <- Nmix.gof.test(best_juvenile, nsim = 100)

cat("\n===== MALE GOF =====\n")
print(male_gof, digits.vals = 4, digits.chisq = 4)
cat("\n===== FEMALE GOF =====\n")
print(female_gof, digits.vals = 4, digits.chisq = 4)
cat("\n===== JUVENILE GOF =====\n")
print(juvenile_gof, digits.vals = 4, digits.chisq = 4)

# ============================================================================
# COEFFICIENT COMPARISON (with different best models per class)
# ============================================================================

best_male <- male_abund_models$gm13
best_female <- female_abund_models$gm11
best_juvenile <- juvenile_abund_models$gm4

male_coefs <- coef(best_male)
female_coefs <- coef(best_female)
juvenile_coefs <- coef(best_juvenile)

# Print each separately with model name
cat("\n===== MALE COEFFICIENTS (model gm13) =====\n")
print(male_coefs)

cat("\n===== FEMALE COEFFICIENTS (model gm11) =====\n")
print(female_coefs)

cat("\n===== JUVENILE COEFFICIENTS (model gm4) =====\n")
print(juvenile_coefs)

# Write each to separate CSVs
write.csv(
  data.frame(Parameter = names(male_coefs), Coefficient = male_coefs, row.names = NULL),
  "Gazelle_TopModel_Coefficients_MALE_gm13.csv", row.names = FALSE
)
write.csv(
  data.frame(Parameter = names(female_coefs), Coefficient = female_coefs, row.names = NULL),
  "Gazelle_TopModel_Coefficients_FEMALE_gm11.csv", row.names = FALSE
)
write.csv(
  data.frame(Parameter = names(juvenile_coefs), Coefficient = juvenile_coefs, row.names = FALSE),
  "Gazelle_TopModel_Coefficients_JUVENILE_gm4.csv", row.names = FALSE
)

# If you want to compare ONLY the abundance-level coefficients that appear in all models:
# Extract just the lam() coefficients (abundance, not detection p())
male_abund_coefs <- male_coefs[grep("^lam", names(male_coefs))]
female_abund_coefs <- female_coefs[grep("^lam", names(female_coefs))]
juvenile_abund_coefs <- juvenile_coefs[grep("^lam", names(juvenile_coefs))]

cat("\n===== ABUNDANCE COEFFICIENTS COMPARISON =====\n")
cat("MALE (gm13):\n")
print(male_abund_coefs)
cat("\nFEMALE (gm11):\n")
print(female_abund_coefs)
cat("\nJUVENILE (gm4):\n")
print(juvenile_abund_coefs)


