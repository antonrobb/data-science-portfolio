# Music Sentiment Analysis — Batch & Real-Time Big Data Pipeline

Sentiment analysis of online music communities across two platforms, combining a
large-scale **batch** pipeline over historical YouTube comments with a **real-time
streaming** pipeline over live Reddit comments. Sentiment is compared across genres
(Hip Hop, Indie, Pop) and two artists (Drake, Kendrick Lamar).

*MSc Data Science coursework, Liverpool John Moores University.*

## What it does

The project demonstrates an end-to-end big data system with two complementary halves:

**Batch (YouTube → PySpark).** Collects historical comments via the YouTube Data API,
then runs a distributed PySpark pipeline: deduplication, language filtering, text
normalisation, dual sentiment scoring, and genre/artist classification, ending with
summary charts and an interactive dashboard.

**Streaming (Reddit → MongoDB → live dashboard).** Ingests an unbounded stream of new
Reddit comments via PRAW, buffers them in MongoDB, and visualises rolling sentiment in
real time. MongoDB decouples ingestion from analysis, so the stream keeps running even if
the analysis kernel restarts.

Both halves score sentiment with **two models** — VADER (tuned for social-media text) and
TextBlob — so the analysis can show where the models agree and where they diverge, rather
than relying on a single tool.

## Files

- `batch_analysis.ipynb` — YouTube collection and the full PySpark batch pipeline, through
  to genre/artist results and an interactive dashboard.
- `streaming_ingestion.ipynb` — the live Reddit → MongoDB streamer and a self-updating
  real-time sentiment plot.
- `live_dashboard.ipynb` — interactive `ipywidgets` views over the live MongoDB buffer
  (per-subreddit, cross-genre, word clouds, and artist comparison).
- `report.pdf` — full write-up: methodology, results, and a discussion of ethical and
  methodological limitations.

## Selected findings

Indie and Pop comments skewed most consistently positive, while Hip Hop was more
polarised. Drake commentary was more emotionally volatile than Kendrick Lamar's, which was
steadier and more positive — a pattern that held across both the batch and streaming data.
The report discusses these alongside the known limitations of lexicon-based sentiment
tools (sarcasm, slang, and AAVE in particular).

## Tech stack

PySpark · Spark SQL · MongoDB · PRAW · YouTube Data API · NLTK (VADER) · TextBlob ·
Matplotlib · ipywidgets · wordcloud

## Running it

The notebooks expect API credentials and services that are not bundled with this repo:

- **Batch:** a YouTube Data API key (placeholder `API_KEY` in `batch_analysis.ipynb`). Once
  comments have been collected to `comments.jsonl`, the rest of the pipeline runs without
  a key.
- **Streaming/dashboard:** a running local MongoDB instance and Reddit API credentials
  (placeholders in `streaming_ingestion.ipynb`). Run the streaming notebook first and
  leave it running, then open the dashboard.

The raw collected comment data is intentionally not included in this repository.