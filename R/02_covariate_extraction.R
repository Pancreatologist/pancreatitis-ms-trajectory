#==============================================================================
# 02 - Baseline characteristics & covariate extraction
# Demographics, BMI, heart rate, comorbidities, medications, labs, surveys
#==============================================================================

library(tidyverse)

# ---- Baseline window: +/- 180 days around pancreatitis onset ----

# ---- Demographics ----
# age_at_pancreatitis, gender, race/ethnicity

# ---- Anthropometrics / vitals ----
# BMI, heart rate, systolic/diastolic BP

# ---- Comorbidities ----
# Charlson comorbidity index components, prior CV disease

# ---- Medications ----
# insulin, metformin, GLP-1 RA, SGLT2i, statins, PERT (enzyme replacement)

# ---- Laboratory biomarkers ----
# albumin, hematocrit, HbA1c, triglycerides, eGFR (baseline values)

# ---- Survey modules ----
# Basics (demographics)
# Lifestyle (smoking, alcohol, diet, physical activity)
# Social (socioeconomic status, social support)

# ---- Fitbit wearable data (optional / sensitivity) ----
# steps, sleep, heart rate (average over first 6 months)

# ---- Output: patient-level baseline table ----
# One row per patient, all baseline covariates

cat("Baseline covariate extraction complete.\n")
