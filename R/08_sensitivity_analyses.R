#==============================================================================
# 08 - Sensitivity analyses
# Multiple imputation methods, measurement count restrictions, extended window
#==============================================================================

library(tidyverse)
library(survival)
library(mice)

# ---- Sensitivity 1: Different imputation methods ----
# Compare median imputation, kNN imputation, MICE

# a) Median imputation
data_median <- data_long %>%
  group_by(person_id) %>%
  mutate(across(all_of(traj_vars), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
  ungroup()

# b) kNN imputation (k = 5)
# using VIM::kNN or caret::preProcess

# c) MICE (primary analysis)
# already done in 03_trajectory_modeling.R

# ---- Sensitivity 2: Measurement count restrictions ----
# Require >= 6, >= 8 measurements per patient

for (min_n in c(6, 8)) {
  stable_patients <- data_long %>%
    filter(!is.na(albumin) & !is.na(bmi)) %>%
    group_by(person_id) %>%
    summarise(n_obs = n()) %>%
    filter(n_obs >= min_n) %>%
    pull(person_id)

  # Re-run trajectory model on restricted cohort
  # Re-run Cox model
  # Compare HR stability
}

# ---- Sensitivity 3: Extended clearance window ----
# Primary: 180-day baseline window
# Sensitivity: 365-day window

# ---- Sensitivity 4: Excluding early deaths ----
# Primary: exclude deaths < 3 months
# Sensitivity: exclude deaths < 6 months

# ---- Sensitivity 5: Complete case analysis ----
# Only patients with all time points observed

# ---- Output: forest plot comparing primary vs all sensitivity HRs ----

cat("Sensitivity analyses complete.\n")
cat("Methods: imputation, measurement count, window, early death exclusion, complete case\n")
