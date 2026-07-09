# Demographic Predictors of Reform UK Vote Share (2024 General Election)

A multiple-regression analysis of which constituency-level demographic and
socioeconomic characteristics best predict Reform UK's share of the vote across the
650 UK parliamentary constituencies in the 2024 general election.

*Independent analysis, originally developed as MSc Data Science coursework at
Liverpool John Moores University and reworked here.*

## Question

Which features of a constituency (age structure, ethnic composition, education,
deprivation, income, region) independently predict how strongly it voted for
Reform UK, and how does the picture change when these correlated variables are
considered together rather than one at a time?

## Approach

- **Multi-source data assembly.** Election results are merged with 2021 census
  demographics across several files, with a deliberately focused set of
  non-redundant explanatory variables chosen to avoid collinearity and overfitting.
- **Univariate then multivariate modelling.** Each variable's marginal association
  with vote share is measured first, then all are combined in a multiple linear
  regression to isolate each variable's independent contribution, with
  multiple-comparison correction and full model diagnostics.
- **Confounding, suppression, and a sign reversal.** The analysis centres on how
  associations change between the univariate and multivariate views: variables that
  stay strong, a deprivation effect that only emerges once other variables are
  controlled for (suppression), and an age effect that reverses sign once
  confounders are removed. It is a concrete illustration of why simple pairwise
  relationships can mislead when predictors are correlated.

## Files

- `election_vote_share.Rmd` — the full analysis source.
- `election_vote_share.pdf` — the knitted report with all code, plots, and results.
- the constituency and demographic CSVs used by the analysis.

## Methods and skills

Multiple linear regression · correlation analysis · one-way ANOVA · multiple-
comparison correction · confounding and suppression effects · model diagnostics and
selection statistics (adjusted R², AIC, PRESS) · multi-source data merging ·
R (`dplyr`, `ggplot2`).

## Data

2024 general election results and 2021 census demographics, sourced from the House
of Commons Library (public data). The analysis is fully reproducible in R from the
original public files.