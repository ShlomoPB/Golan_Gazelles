# Supplementary table quality control

Generated from saved final-analysis objects and summary files without refitting models.

## Numerical validation

| Table | Result | Manuscript value | Reproduced value | Within rounding | Source | Note |
| --- | --- | --- | --- | --- | --- | --- |
| S1 | Juvenile – adult female estimate | 0.863 | 0.863347425311037 | YES | gazelle_lf2_NMix_NB_open_wolf_jackal_interaction_order-fawn-male-female_2000b_25bl_10th_3ch_fit.rds |  |
| S1 | Juvenile – adult female lower CrI | 0.523 | 0.523338424639095 | YES | gazelle_lf2_NMix_NB_open_wolf_jackal_interaction_order-fawn-male-female_2000b_25bl_10th_3ch_fit.rds |  |
| S1 | Juvenile – adult female upper CrI | 0.998 | 0.998448827940644 | YES | gazelle_lf2_NMix_NB_open_wolf_jackal_interaction_order-fawn-male-female_2000b_25bl_10th_3ch_fit.rds |  |
| S1 | Juvenile – adult male estimate | 0.657 | 0.657368207277351 | YES | gazelle_lf2_NMix_NB_open_wolf_jackal_interaction_order-fawn-male-female_2000b_25bl_10th_3ch_fit.rds |  |
| S1 | Juvenile – adult male lower CrI | 0.435 | 0.435370942256979 | YES | gazelle_lf2_NMix_NB_open_wolf_jackal_interaction_order-fawn-male-female_2000b_25bl_10th_3ch_fit.rds |  |
| S1 | Juvenile – adult male upper CrI | 0.806 | 0.806226816547003 | YES | gazelle_lf2_NMix_NB_open_wolf_jackal_interaction_order-fawn-male-female_2000b_25bl_10th_3ch_fit.rds |  |
| S1 | Adult female – adult male estimate | 0.907 | 0.906536431528505 | YES | gazelle_lf2_NMix_NB_open_wolf_jackal_interaction_order-fawn-male-female_2000b_25bl_10th_3ch_fit.rds |  |
| S1 | Adult female – adult male lower CrI | 0.712 | 0.712480660992152 | YES | gazelle_lf2_NMix_NB_open_wolf_jackal_interaction_order-fawn-male-female_2000b_25bl_10th_3ch_fit.rds |  |
| S1 | Adult female – adult male upper CrI | 1 | 0.999565855843029 | YES | gazelle_lf2_NMix_NB_open_wolf_jackal_interaction_order-fawn-male-female_2000b_25bl_10th_3ch_fit.rds |  |
| S2 | Wolf–jackal site correlation | -0.183 | -0.182689289473995 | YES | SiteCovs_060324.csv |  |
| S3 | Juvenile low-wolf jackal ratio | 12.8 | 12.7932227205181 | YES | publication_wolf_jackal_interaction_nmix_spabundance_fits.rds |  |
| S3 | Adult female low-wolf jackal ratio | 12.2 | 12.1551174283531 | YES | publication_wolf_jackal_interaction_nmix_spabundance_fits.rds |  |
| S3 | Adult male low-wolf jackal ratio | 7.5 | 7.47361902050348 | YES | publication_wolf_jackal_interaction_nmix_spabundance_fits.rds |  |
| S3 | Adult female attenuation ratio-of-ratios | 7.97 | 7.970869863168 | YES | publication_wolf_jackal_interaction_nmix_spabundance_fits.rds |  |
| S3 | Adult male attenuation ratio-of-ratios | 4.9 | 4.89804401386477 | YES | publication_wolf_jackal_interaction_nmix_spabundance_fits.rds |  |
| S5 | Candidate AIC: Wolf × jackal; no survey-period term | 444.913 | 444.912522445422 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S5 | Candidate AIC: Additive wolf + jackal; no survey-period term | 443.222 | 443.222403879682 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S5 | Candidate AIC: Wolf × jackal × survey period | 444.194 | 444.194079240363 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S5 | Candidate AIC: Additive predator slopes varying by survey period | 444.105 | 444.104826538834 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S5 | Candidate AIC: Wolf × jackal + common survey-period shift | 442.722 | 442.721762661205 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S5 | Candidate AIC: Additive wolf + jackal + common survey-period shift | 441.075 | 441.075486303177 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S5 | Add common Early/Late shift to additive model LR | 4.147 | 4.14691757650428 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S5 | Add common Early/Late shift to additive model P | 0.0417 | 0.0417103873420457 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S5 | Add wolf × jackal to retained model LR | 0.35 | 0.353723641972238 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S5 | Add wolf × jackal to retained model P | 0.552 | 0.552012797700972 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S5 | Allow additive predator slopes to vary by Early/Late LR | 0.971 | 0.970659764343168 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S5 | Allow additive predator slopes to vary by Early/Late P | 0.616 | 0.615494122733071 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S6 | Late vs Early survey OR | 1.57 | 1.56924020015271 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S6 | Late vs Early survey lower CI | 1.01 | 1.01106099913428 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S6 | Late vs Early survey upper CI | 2.44 | 2.43557491376272 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S6 | Late vs Early survey P | 0.045 | 0.044536087827315 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S6 | Open-space cover OR | 1.55 | 1.55351958270117 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S6 | Open-space cover lower CI | 1.04 | 1.03998513938583 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S6 | Open-space cover upper CI | 2.32 | 2.32063228832414 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S6 | Open-space cover P | 0.031 | 0.0314413552923108 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S6 | Wolf activity OR | 0.93 | 0.928869284909848 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S6 | Wolf activity lower CI | 0.7 | 0.700654644598399 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S6 | Wolf activity upper CI | 1.23 | 1.23141715408662 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S6 | Wolf activity P | 0.608 | 0.608004624281759 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S6 | Jackal activity OR | 0.94 | 0.942659764417875 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S6 | Jackal activity lower CI | 0.63 | 0.63104519621082 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S6 | Jackal activity upper CI | 1.41 | 1.40815180400406 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S6 | Jackal activity P | 0.773 | 0.7730546248299 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S6 | Wolf × jackal sensitivity OR | 0.86 | 0.859651102065014 | YES | publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds |  |
| S7 | Juvenile selected model | Wolf × jackal | Wolf × jackal | YES | definitive_hierarchical_diel_ecological_model_comparison.csv |  |
| S7 | Adult female selected model | Wolf × jackal | Wolf × jackal | YES | definitive_hierarchical_diel_ecological_model_comparison.csv |  |
| S7 | Adult male selected model | Wolf × jackal | Wolf × jackal | YES | definitive_hierarchical_diel_ecological_model_comparison.csv |  |
| S7 | Pooled gazelles selected model | Wolf-only | Wolf-only | YES | definitive_hierarchical_diel_ecological_model_comparison.csv |  |
| S7 | Juvenile additive-model ΔAIC | 9.666 | 9.66555973452705 | YES | definitive_hierarchical_diel_ecological_model_comparison.csv |  |
| S7 | Adult female additive-model ΔAIC | 2.495 | 2.49458567954389 | YES | definitive_hierarchical_diel_ecological_model_comparison.csv |  |
| S7 | Adult male additive-model ΔAIC | 3.262 | 3.26156645745141 | YES | definitive_hierarchical_diel_ecological_model_comparison.csv |  |
| S7 | Pooled gazelles additive-model ΔAIC | 0.196 | 0.195946064604868 | YES | definitive_hierarchical_diel_ecological_model_comparison.csv |  |
| S8 | Juvenile interaction EDF | 2.067 | 2.06728414486985 | YES | definitive_hierarchical_diel_model_summary.csv; definitive_hierarchical_diel_effect_stability.csv |  |
| S8 | Juvenile interaction support | 7 | 7 | YES | definitive_hierarchical_diel_model_summary.csv; definitive_hierarchical_diel_effect_stability.csv |  |
| S8 | Adult-female wolf EDF | 3.916 | 3.91620164572939 | YES | definitive_hierarchical_diel_model_summary.csv; definitive_hierarchical_diel_effect_stability.csv |  |
| S8 | Adult-male interaction EDF | 2.991 | 2.99143793781118 | YES | definitive_hierarchical_diel_model_summary.csv; definitive_hierarchical_diel_effect_stability.csv |  |
| S8 | Pooled selected wolf-only EDF | 3.041 | 3.86964278459021 | NO | definitive_hierarchical_diel_model_summary.csv; definitive_hierarchical_diel_effect_stability.csv | Manuscript/prompt value 3.041 is from the full interaction robustness model; the final selected wolf-only model gives EDF 3.870 |
| S9 | Wolf × time omnibus P | 0.0778 | 0.0777866427862999 | YES | joint_demographic_penalty_calibrated_outputs/penalty_calibrated_omnibus_tests.csv |  |
| S9 | Jackal × time omnibus P | 0.7066 | 0.706564848318109 | YES | joint_demographic_penalty_calibrated_outputs/penalty_calibrated_omnibus_tests.csv |  |
| S9 | Wolf × jackal × time omnibus P | 0.006 | 0.00596393847454235 | YES | joint_demographic_penalty_calibrated_outputs/penalty_calibrated_omnibus_tests.csv |  |
| S9 | Juvenile – Adult female Holm P | 0.0121 | 0.0120739310273681 | YES | joint_demographic_penalty_calibrated_outputs/penalty_calibrated_pairwise_tests.csv |  |
| S9 | Juvenile – Adult male Holm P | 0.0111 | 0.0111119607447454 | YES | joint_demographic_penalty_calibrated_outputs/penalty_calibrated_pairwise_tests.csv |  |
| S9 | Adult female – Adult male Holm P | 0.2767 | 0.276733355215314 | YES | joint_demographic_penalty_calibrated_outputs/penalty_calibrated_pairwise_tests.csv |  |

## Structural checks

- Ten unique table numbers generated: **YES** (`S1`–`S10`).
- Tables generated in numerical order: **YES**.
- Appendix S1 contains Tables S1–S4: **YES**.
- Appendix S2 contains Tables S5–S6: **YES**.
- Appendix S3 contains Tables S7–S10: **YES**.
- Current manuscript appendix cross-references use the final S1/S2/S3 map and contain no Appendix S4/S5 references: **NO OR NOT CHECKED**.

## Unresolved or material issues

1. **Pooled diel EDF mismatch:** the selected pooled wolf-only GAM gives EDF = 3.869643 (Table S8). The value 3.041 in the current manuscript/prompt comes from the pooled full-interaction robustness fit and must not be reported as the selected-model EDF.
2. **Predator-index provenance:** the current manuscript records the final state/detection specifications and empirical-Bayes extraction, but the original wolf/jackal fit objects and candidate-fit statistics are absent. Table S2 therefore does not report unreproducible AIC, WAIC, or convergence statistics for those index models.
3. **Early/Late boundary wording:** code classifies round midpoints on or before 17 October 2020 as Early and later dates as Late; the manuscript says before 17 October versus 17 October onward. No midpoint occurred exactly on 17 October, so estimates are unaffected, but wording should be harmonized.
4. **Pairwise joint-GAM robustness:** the omnibus wolf × jackal demographic heterogeneity was audited in 7/7 specifications; pairwise juvenile–adult contrasts were tested only in the primary apparent-solar model and must not be described as 7/7 robust.

## Manuscript cross-reference occurrences

- Current manuscript DOCX was not found at the configured path; cross-references were not scanned.

The manuscript was inspected only; it was not modified.
