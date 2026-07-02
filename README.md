Anton Robb — Data Science & Statistical Modelling Portfolio

A collection of projects from my MSc Data Science (Liverpool John Moores University), spanning Bayesian statistics, machine learning, and big data engineering. My focus is on probabilistic modelling, simulation, and turning messy real-world data into reliable insight, with a particular interest in applying statistical methods to sport.

Background: MSc Data Science (predicted Distinction) · BSc Mathematics with Statistics (University of Nottingham) · 2.5 years as a Data Analyst at BDO UK.

## Projects

### 1. Target Practice — Monte Carlo Simulation of Archery & Darts
**R · Monte Carlo simulation · Bayesian inference · MCMC**

A Monte Carlo study modelling player accuracy in archery and darts using a bivariate normal distribution, computing expected scores and in-play win probabilities under formal competition rules, and recovering hidden player parameters from observed shots via Metropolis-Hastings MCMC.

* Models shot location as a bivariate normal and converts to polar coordinates to assign scores
* Simulates full matches to compute live, turn-by-turn win probabilities
* Implements Bayesian inverse modelling (hand-coded Metropolis-Hastings) to recover a player's aim and spread from their scores
* Co-authored with Ken Kahuthu. Awarded 85% (highest in cohort).

### 2. Human Activity Recognition — Machine Learning
**Python · XGBoost · feature engineering · scikit-learn**

A classification project identifying six human activities from smartphone accelerometer data. Diagnosed orientation bias in a 1D-CNN approach and pivoted to physics-informed feature engineering (FFT, jerk derivatives), training a gradient-boosted model that generalised strongly to unseen users.

* Engineered frequency and derivative-based features from raw sensor signals
* Trained and tuned an XGBoost classifier with near-identical public/private leaderboard scores (0.925 / 0.922), evidencing strong generalisation
* Co-authored with Ken Kahuthu. Awarded 77%.

### 3. Surface Effects on Playing Style in 2013 Grand Slam Tennis — Exploratory Data Analysis
**R · ggplot2 · dplyr · data wrangling · statistical visualisation**

An exploratory analysis of how court surface, gender, and tournament round shape professional tennis playing style, using 1,010 match-level statistics from the 2013 Grand Slams (UCI Machine Learning Repository).

* Recoded raw player labels (P1/P2) to winner/loser for interpretable margin statistics, sharpening downstream signal by 2–5× (R² improvements on key relationships)
* Built a data quality assessment visualising missingness by tournament, identifying and handling sparse US Open round information
* Designed a multi-stage analysis: simple exploratory plots, complex multi-panel figures with faceting, regression slopes, and annotations
* Found that winning margins in winners and unforced errors are positively related everywhere, but the relationship is weakest on grass (R² = 0.33), supporting the hypothesis that faster surfaces reward serve-decided points

### 4. Sentiment in Online Music Communities — Big Data
**PySpark · Spark · MongoDB · NLP (VADER, TextBlob) · real-time dashboards**

An end-to-end big data system analysing sentiment in online music communities across two platforms: a batch pipeline processing historical YouTube comments, and a real-time streaming pipeline monitoring Reddit music subreddits. Sentiment is compared across genres (Hip Hop, Indie, Pop) and artists.

* Built a PySpark batch pipeline for cleaning, language filtering, deduplication, and sentiment scoring of YouTube comments collected via the YouTube Data API
* Built a live Reddit streaming pipeline (PRAW) with comments buffered in MongoDB and scored on retrieval
* Scored sentiment with two models (VADER and TextBlob) and surfaced results through an interactive dashboard with rolling sentiment curves and word clouds
* Documented methodology, limitations, and ethical considerations (privacy, representativeness, algorithmic and linguistic bias)

## Skills demonstrated

**Languages:** R, Python

**Statistics & modelling:** Bayesian inference, MCMC, Monte Carlo simulation, hidden Markov models, regression

**Machine learning:** XGBoost / gradient boosting, feature engineering, neural networks

**Big data:** PySpark, Spark Streaming, MongoDB, NLP, large-scale data pipelines

**Data analysis:** exploratory data analysis (EDA), data wrangling, statistical visualisation

## Contact

* **LinkedIn:** https://www.linkedin.com/in/antonrobb
* **Email:** ntnrbb@gmail.com
