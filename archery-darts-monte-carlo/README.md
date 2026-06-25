# Target Practice — Monte Carlo Simulation of Archery & Darts

A simulation study of archery and darts that models player accuracy probabilistically,
computes expected scores and in-play win probabilities under formal competition rules, and
recovers a player's hidden skill parameters from their shots using Bayesian inference.

*Co-authored with Ken Kahuthu. MSc Data Science coursework, Liverpool John Moores
University. Awarded 85% (highest in cohort).*

## Core idea

A player's shot is modelled as a draw from a **bivariate normal distribution**: the mean
is where they aim, and the covariance encodes their spread and any horizontal/vertical
correlation. Converting each shot to polar coordinates (radius and angle) maps it onto the
scoring rings of an archery target or the sectors of a dartboard. From this single model,
the project builds outward in two directions — *forward* simulation (given a player's
skill, how do they score and how often do they win?) and *inverse* inference (given only a
player's shots, can we recover how they were aiming?).

## What's in it

- **Archery simulation** — score players against a ringed target, with Monte Carlo
  confidence intervals on mean scores and head-to-head win probabilities.
- **World Archery set-play matches** — simulate full matches under the official
  set-play scoring rules (including shoot-offs) to estimate match-win probabilities.
- **Parameter sensitivity** — quantify how aiming bias and spread each affect expected
  performance: what actually makes a better archer?
- **Darts simulation** — a full dartboard model (sectors, trebles, doubles, bull),
  simulated games of 301 (double-out), and an analysis of **optimal targeting** — when a
  player should aim for treble 20 versus a safer target, with the crossover computed as a
  function of accuracy.
- **Live in-play odds** — a ball-by-ball match simulator that produces a live win-
  probability "worm chart" as a match unfolds.
- **Bayesian inverse modelling (MCMC)** — the centrepiece: a hand-coded Metropolis-
  Hastings sampler that recovers a player's aim point and consistency from the spatial
  coordinates of their darts alone, characterising not just point estimates but the full
  posterior uncertainty. Includes burn-in handling and trace-plot convergence diagnostics.

## Methods and skills

Monte Carlo simulation · Bayesian inference · Metropolis-Hastings MCMC (implemented from
scratch) · bivariate normal modelling · convergence diagnostics · expected-value
optimisation · probabilistic match modelling · R, with `mvtnorm` and `ggplot2`.

## Files

- `target_practice.Rmd` — the full R Markdown report: simulation code, mathematical
  setup, visualisations, and written interpretation throughout.
- `target_practice.pdf` — the knitted report.

## Why it's interesting

The project connects a clean statistical model to applied sports questions on both ends.
The forward simulations answer competition and strategy questions (win probabilities,
optimal targeting, live odds), while the MCMC section solves the harder inverse problem —
inferring a player's hidden skill from observed outcomes — which has natural applications
in player profiling and performance analytics.