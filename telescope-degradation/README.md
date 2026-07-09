# Modelling Optical Degradation of the Liverpool Telescope

A time-series analysis of a decade of nightly reflectivity measurements from the
Liverpool Telescope, modelling how its mirrors degrade between cleanings,
separating genuine measurements from cloud-affected ones, testing whether the
degradation rate changes over time, and forecasting when cleaning will next be
required.

*Independent analysis, originally developed as MSc Data Science coursework at
Liverpool John Moores University and reworked here.*

## Problem

Every clear night, the telescope measures its optical transmission (the fraction
of incident light reaching the focus). Contamination causes this to decay
exponentially between cleanings, which restore it in sudden jumps. The data are
complicated by cloud-affected nights, which depress the measurement without any
real change in the mirrors, and by strong autocorrelation between nearby nights.

## Approach

- **Segmented exponential-decay model.** Transmission is log-transformed and
  modelled as a segmented linear decay, with cleaning events identified as abrupt
  upward steps and each segment given its own recovery level.
- **Iterative, one-sided outlier removal.** Cloud-affected points only ever reduce
  the measurement, so they are removed with a one-sided cut applied iteratively:
  the model is refitted after each pass so it tightens around the genuine data and
  progressively exposes cloud points a single pass would miss.
- **Autocorrelation-aware inference.** Because nightly measurements are strongly
  correlated over short timescales, testing whether the degradation rate changes
  over time is done on monthly-aggregated data, giving an honest effective sample
  size and avoiding the over-confident p-values a naive nightly regression
  produces.
- **Forecasting.** The final model predicts transmission on a target date and the
  date at which transmission will next fall below the operational threshold, both
  with confidence intervals.

## Files

- `telescope_degradation.Rmd` — the full analysis source.
- `telescope_degradation.pdf` — the knitted report with all code, plots, and results.
- `telescope.csv` — the nightly transmission data.

## Methods and skills

Time-series modelling · exponential-decay / log-linear models · segmented
regression · iterative robust outlier removal · handling autocorrelation via
temporal aggregation · model diagnostics · forecasting with confidence intervals ·
R (`ggplot2`, `dplyr`, `lubridate`).

## Running it

With `telescope.csv` in the same folder, knit `telescope_degradation.Rmd` in R to
reproduce the models, plots, and forecasts.

## Notes

Telescope data provided by Liverpool John Moores University and included with
permission. This project began as university coursework and has been reworked as an
independent analysis.