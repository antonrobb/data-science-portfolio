# Anton Robb — Data Science & Statistical Modelling Portfolio

A collection of projects from my MSc Data Science (Liverpool John Moores University), spanning Bayesian statistics, stochastic modelling, machine learning, and big data engineering. My focus is on probabilistic modelling, simulation, and turning messy real-world data into reliable insight, with a particular interest in applying statistical methods to sport.

Background: MSc Data Science (Distinction) · BSc Mathematics with Statistics (University of Nottingham) · over 2 years as a Data Analyst at BDO UK.

## Projects

### 1. Bayesian Inference for a Two-State Stochastic Velocity-Jump Model — MSc Dissertation
**R · Bayesian inference · adaptive MCMC · hidden Markov models · likelihood derivation**

[View project →](velocity-jump-bayesian-inference/)

A particle moves at one of two constant velocities, switching between them at random times, and is observed only through noisy position measurements taken at a fixed frame rate. The switching itself is never seen. This project derives an exact likelihood for that model, implements it alongside three truncated approximations from the published literature, and compares all four under a common Bayesian sampler.

* Derived an exact likelihood by extending a published family of truncated approximations to its limit, using the occupation-time density to collapse an infinite sum over switch counts into a closed-form expression
* Verified the derivation three ways, each computed by a route independent of the formula itself: against 400,000 simulated paths, against exact normalisation, and against a separately implemented approximation onto which it converges
* Built the full pipeline across seven base-R modules: Gillespie simulation, four likelihood implementations, and adaptive Metropolis MCMC recovering velocity, switching-rate and measurement-noise parameters
* Found that the zero-switch approximation is structurally misspecified, with biases that persist as the sample grows, while the truncated approximations' biases are sample-size artefacts that vanish by N = 1600. More data cannot fix a misspecified likelihood
* MSc dissertation, supervised by Dr Ivo Siekmann. Awarded 82%, the highest mark in the cohort. Currently being prepared for journal submission.

### 2. Target Practice — Monte Carlo Simulation of Archery & Darts
**R · Monte Carlo simulation · Bayesian inference · MCMC**

[View project →](archery-darts-monte-carlo/)

A Monte Carlo study modelling player accuracy in archery and darts using a bivariate normal distribution, computing expected scores and in-play win probabilities under formal competition rules, and recovering hidden player parameters from observed shots via Metropolis-Hastings MCMC.

* Models shot location as a bivariate normal and converts to polar coordinates to assign scores
* Simulates full matches to compute live, turn-by-turn win probabilities
* Implements Bayesian inverse modelling (hand-coded Metropolis-Hastings) to recover a player's aim and spread from their scores
* Co-authored with Ken Kahuthu. Awarded 85% (highest in cohort).

### 3. Human Activity Recognition — Machine Learning
**Python · XGBoost · feature engineering · scikit-learn**

[View project →](human-activity-recognition/)

A classification project identifying six human activities from smartphone accelerometer data. Diagnosed orientation bias in a 1D-CNN approach and pivoted to physics-informed feature engineering (FFT, jerk derivatives), training a gradient-boosted model that generalised strongly to unseen users.

* Engineered frequency and derivative-based features from raw sensor signals
* Trained and tuned an XGBoost classifier with near-identical public/private leaderboard scores (0.925 / 0.922), evidencing strong generalisation
* Co-authored with Ken Kahuthu. Awarded 77%.

### 4. Surface Effects on Playing Style in 2013 Grand Slam Tennis — Exploratory Data Analysis
**R · ggplot2 · dplyr · data wrangling · statistical visualisation**

[View project →](tennis-surface-analysis/)

An exploratory analysis of how court surface, gender, and tournament round shape professional tennis playing style, using 1,010 match-level statistics from the 2013 Grand Slams (UCI Machine Learning Repository).

* Recoded raw player labels (P1/P2) to winner/loser for interpretable margin statistics, sharpening downstream signal by 2–5× (R² improvements on key relationships)
* Built a data quality assessment visualising missingness by tournament, identifying and handling sparse US Open round information
* Designed a multi-stage analysis: simple exploratory plots, complex multi-panel figures with faceting, regression slopes, and annotations
* Found that winning margins in winners and unforced errors are positively related everywhere, but the relationship is weakest on grass (R² = 0.33), supporting the hypothesis that faster surfaces reward serve-decided points

### 5. Modelling Optical Degradation of the Liverpool Telescope — Time Series
**R · time-series modelling · segmented regression · forecasting**

[View project →](telescope-degradation/)

A time-series analysis of a decade of nightly reflectivity measurements from the Liverpool Telescope, modelling how the mirrors degrade between cleanings, separating genuine measurements from cloud-affected ones, and forecasting when cleaning will next be required.

* Modelled exponential decay between cleaning events using a segmented log-linear model, with cleanings identified as abrupt upward steps in transmission
* Removed cloud-affected measurements with an iterative, one-sided outlier procedure, refitting after each pass to progressively expose points a single pass would miss
* Tested for a changing degradation rate on monthly-aggregated data to account for the strong autocorrelation of nightly measurements, giving honest standard errors rather than the over-confident results a naive nightly regression produces
* Forecast future transmission and the date it will next fall below the operational threshold, with confidence intervals

### 6. Demographic Predictors of Reform UK Vote Share — Multiple Regression
**R · multiple regression · ANOVA · statistical inference**

[View project →](election-vote-share/)

A multiple-regression analysis of which constituency-level demographic and socioeconomic characteristics best predict Reform UK's vote share across the 650 UK parliamentary constituencies in the 2024 general election.

* Merged 2024 election results with 2021 census demographics across multiple sources, selecting a focused, non-redundant set of explanatory variables
* Built a multiple linear regression with model diagnostics, multiple-comparison correction, and selection statistics (adjusted R², AIC, PRESS)
* Uncovered a suppression effect (deprivation becoming significant only after controlling for other variables) and a sign reversal (the older-population effect flipping between the univariate and multivariate views), illustrating why pairwise relationships mislead when predictors are correlated

### 7. Sentiment in Online Music Communities — Big Data
**PySpark · Spark · MongoDB · NLP (VADER, TextBlob) · real-time dashboards**

[View project →](sentiment-analysis-big-data/)

An end-to-end big data system analysing sentiment in online music communities across two platforms: a batch pipeline processing historical YouTube comments, and a real-time streaming pipeline monitoring Reddit music subreddits. Sentiment is compared across genres (Hip Hop, Indie, Pop) and artists.

* Built a PySpark batch pipeline for cleaning, language filtering, deduplication, and sentiment scoring of YouTube comments collected via the YouTube Data API
* Built a live Reddit streaming pipeline (PRAW) with comments buffered in MongoDB and scored on retrieval
* Scored sentiment with two models (VADER and TextBlob) and surfaced results through an interactive dashboard with rolling sentiment curves and word clouds
* Documented methodology, limitations, and ethical considerations (privacy, representativeness, algorithmic and linguistic bias)

## Skills demonstrated

**Languages:** R, Python

**Statistics & modelling:** Bayesian inference, MCMC and adaptive MCMC, likelihood derivation, Monte Carlo simulation, hidden Markov models, stochastic processes, regression, time-series analysis

**Machine learning:** XGBoost / gradient boosting, feature engineering, neural networks

**Big data:** PySpark, Spark Streaming, MongoDB, NLP, large-scale data pipelines

**Data analysis:** exploratory data analysis (EDA), data wrangling, statistical visualisation

## Contact

* **LinkedIn:** https://www.linkedin.com/in/antonrobb
* **Email:** antonrobb00@gmail.com
