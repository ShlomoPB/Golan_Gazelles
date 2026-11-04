# Results

## Detection process

Detection varied among camera stations and survey rounds, supporting explicit modeling of the observation process in all analyses. In the primary class-specific abundance models, detection was modeled as a function of sampling effort (`log_days_operational`) and cattle activity. In the exploratory joint latent-factor models, the broader detection structure including path type, camera model, vegetation structure, area behind the camera, and sampling effort outperformed the design-only alternative and was retained for the final joint analysis.

## Class-specific responses to habitat, predator, and anthropogenic covariates

Demographic classes differed in their responses to ecological pressures, but predator covariates emerged as the most consistent predictors across the population (Table 1; Fig. 1). In the fawn model, open habitat was the strongest positive predictor of abundance (`Open.Spaces.z = 1.05`, 95% CrI `0.53` to `1.61`), indicating that juveniles were most abundant at more open sites. Jackal abundance showed a weaker negative association with fawn abundance (`null.jackal.abundance.z = -0.53`, 95% CrI `-1.15` to `0.07`), whereas boar abundance, human abundance, and wolf abundance were all uncertain.

Adult males showed the clearest response to infrastructure and mesopredator pressure (Table 1; Fig. 1). Male abundance declined with increasing built and road cover (`Built.and.Roads.z = -0.64`, 95% CrI `-1.04` to `-0.24`) and with increasing jackal abundance (`null.jackal.abundance.z = -0.73`, 95% CrI `-1.11` to `-0.36`). The wolf effect was positive but only borderline supported (`Null.Wolf.Abundance.z = 0.26`, 95% CrI `-0.01` to `0.52`), and neither boar abundance nor human abundance showed a clearly supported independent effect.

Adult females also responded most strongly to predator covariates (Table 1; Fig. 1). Female abundance declined with jackal abundance (`null.jackal.abundance.z = -0.81`, 95% CrI `-1.34` to `-0.30`) and increased with wolf abundance (`Null.Wolf.Abundance.z = 0.42`, 95% CrI `0.01` to `0.84`). In contrast, the reserve covariate was weak and not clearly supported in the primary female model (`Nature.Reserves.z = -0.10`, 95% CrI `-0.55` to `0.35`). As in the other classes, boar abundance and human abundance remained uncertain.

Taken together, the class-specific abundance models indicate that predator pressure was the dominant common axis of response, especially the negative association with jackals. Habitat context and human-associated landscape structure differed more strongly among classes: fawns were most strongly associated with open habitat, males were most negatively associated with infrastructure, and females were most strongly structured by predator covariates.

## Cross-class comparison

The strongest cross-class pattern was the consistency of the jackal effect (Table 1; Fig. 1). Jackal abundance was negative in all three class-specific models and was the clearest shared predictor across demographic classes. By contrast, wolves showed a more class-dependent pattern, with positive associations that were weak for fawns but stronger for adults, especially females. Human abundance and boar abundance showed consistent directional tendencies, with human coefficients tending to be negative and boar coefficients tending to be positive, but neither covariate was clearly supported in any class.

These results suggest that pooled population-level analyses would obscure biologically important demographic structure. Juveniles responded most strongly to habitat openness, whereas adults responded more strongly to predator gradients and, in males, to human-associated infrastructure.

## Sensitivity and joint-model results

Sensitivity analyses based on occupancy models recovered the same broad ecological story, particularly the importance of open habitat for fawns and the negative effect of jackals across classes. The exploratory latent-factor multispecies N-mixture model also broadly recovered negative jackal effects for all classes and positive wolf effects for adults, indicating that the core ecological inference was not an artifact of the single-class abundance framework.

However, the latent-factor joint model was less suitable as the main inferential model (Table 2). Although its species-level slope directions were broadly similar to those of the primary class-specific models, posterior predictive checks remained poor, especially for males and females, and the final residual cross-class correlations were weaker than in earlier exploratory fits. In the longer class-threat latent-factor model, posterior correlations among site-level log latent mean abundances were `0.19` (`-0.07` to `0.51`) for fawn-male, `0.26` (`-0.02` to `0.53`) for fawn-female, and `0.18` (`-0.09` to `0.51`) for male-female. These results suggest some shared structure among classes, but the strongest evidence for class-specific threat responses comes from the separate non-spatial negative-binomial abundance models.

## Summary

Overall, demographic classes did not respond uniformly to landscape and predator gradients. Jackals were the most consistent negative correlate of gazelle abundance across all classes. Wolves were positively associated with adult abundance, especially in females. Habitat openness was most important for fawns, whereas built and road cover most strongly constrained males. This demographic differentiation in threat response suggests that conservation assessments based on pooled population responses may miss ecologically important class-specific patterns.
