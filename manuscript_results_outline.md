# Results Outline

## 1. Primary Class-Specific Predator Responses

Lead with the finalized separate-class negative-binomial `spAbundance` N-mixture models. Report that jackal activity was negatively associated with every class, whereas positive wolf associations and positive wolf by jackal interactions were supported for adult females and males but not fawns. Explain conditional effects rather than interpreting interaction-model main effects as unconditional averages.

Key interpretation:
- The negative jackal association weakened as wolf activity increased for adult females and males.
- Equivalently, the positive wolf association strengthened where jackal activity was high.
- Fawns showed a negative jackal gradient but uncertain wolf and interaction coefficients.
- These observational patterns support a two-predator-context hypothesis but do not prove wolf-mediated jackal suppression.

Topic sentence:
Adult gazelle responses depended on the joint predator context, whereas fawn relative abundance was associated primarily with lower jackal activity.

## 2. Residual Demographic-Class Correlation

Present the focused two-factor negative-binomial `lfMsNMix` as a secondary joint analysis. State that the fawn, adult-male, and adult-female latent relative-abundance or site-use processes retained strong positive residual correlations after accounting for open-space cover, wolf activity, jackal activity, and their interaction.

Key results:
- Fawn-male: 0.66 (95% CrI = 0.43-0.81).
- Fawn-female: 0.86 (0.51-1.00).
- Male-female: 0.91 (0.73-1.00).
- Free latent loadings converged (`Rhat` = 1.006-1.009), but abundance and detection intercepts did not and were excluded from inference.
- Site-grouped posterior predictive checks were not extreme, although ungrouped checks indicated replicate-level lack of fit.
- The earlier broad additive model produced similar correlations, indicating stability to the focused covariate specification.

Do not use ecological slopes from the joint model as replacements for the finalized separate-class results. Describe the correlations as shared unexplained spatial structure rather than causal demographic interactions.

Topic sentence:
Despite their different predator responses, all three demographic classes retained strong residual spatial association after accounting for measured habitat and predator gradients.

## 3. Framework Sensitivity

Report that matched `unmarked` negative-binomial N-mixture models reproduced the class-specific predator-interaction pattern. Occupancy and independent-event analyses may be cited as additional sensitivity evidence where their estimands align, but they should not be presented as independent replication.

Topic sentence:
Agreement between Bayesian and maximum-likelihood N-mixture implementations indicated that the principal predator gradients were not specific to one fitting framework.

## 4. Class-Specific Synthesis

Fawns showed a negative jackal association but uncertain wolf and interaction responses. Adult males and females shared negative jackal, positive wolf, and positive interaction coefficients. The latent-correlation analysis indicates that these demographic differences occurred within a population whose classes nevertheless showed strongly coupled spatial patterns.

Topic sentence:
Shared demographic spatial structure coexisted with class-specific predator responses, particularly the weaker wolf association among fawns.

## 5. Closing Synthesis

End with the following linked conclusions:
- jackal activity was the most consistent negative predator correlate;
- wolf associations were positive for adults and conditional on jackal activity;
- fawn responses differed from those of adults;
- demographic classes nevertheless shared strong residual spatial structure;
- class-pooled analyses would obscure meaningful demographic differences.

Topic sentence:
Separate-class and joint analyses together revealed a spatially coupled demographic population whose classes differed in their responses to sympatric predators.
