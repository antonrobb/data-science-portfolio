# Human Activity Recognition from Smartphone Accelerometer Data

Classifying physical activities from raw tri-axial accelerometer signals using engineered
time- and frequency-domain features and a gradient-boosted tree model (XGBoost).

*Co-authored with Ken Kahuthu. MSc Data Science coursework, Liverpool John Moores University.*

## Problem

Given short snippets of smartphone accelerometer data (`x`, `y`, `z` over time), predict
which activity the user was performing. The central challenge is **generalising to people
the model has never seen**: raw signals vary with how each person moves and how the phone
is oriented, so a model can easily memorise individuals rather than learn the activities.

## Approach

- **Feature engineering (140 features per snippet).** Each snippet is summarised by
  statistical and spectral features over eight derived signals — the three raw axes, an
  orientation-invariant magnitude, and the jerk (first derivative) of each. Frequency-domain
  features from the FFT capture the periodic rhythm that distinguishes activities such as
  walking and jogging from static ones.
- **User-grouped validation.** Cross-validation uses `GroupKFold` grouped by user, so every
  fold is validated on users held entirely out of training. Reported accuracy is therefore an
  honest estimate of unseen-user performance, not inflated by leakage across the split.
- **Regularised XGBoost.** A conservative configuration (low learning rate with many trees and
  early stopping, shallow trees, subsampling, L1/L2 penalties) guards against overfitting on a
  limited number of users.

## Files

- `human_activity_recognition.ipynb` — the full pipeline: feature engineering, user-grouped
  cross-validation, and prediction, with explanatory notes throughout.

## Notes

This was the final competition entry. An earlier variant explored SHAP-based feature selection
down to ~69 features; it achieved slightly higher raw leaderboard accuracy but generalised less
consistently between the public and private splits, so the full 140-feature model shown here was
preferred. Earlier exploration also tested a 1D-CNN on the raw signals, which proved sensitive to
device orientation and motivated the move to the orientation-invariant engineered features used here.

## Running it

The notebook expects the dataset CSVs (`signals*.csv`, `metadata*.csv`, and the Kaggle test files)
in the same directory. With those present, run the cells top to bottom to reproduce the feature
matrices, cross-validation scores, and the `predictions_AK-67.csv` submission file.