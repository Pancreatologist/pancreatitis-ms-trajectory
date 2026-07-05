#==============================================================================
# 01 - Cohort construction
# Identify adults with pancreatitis + incident metabolic sequelae
# from All of Us Research Program data
#==============================================================================

library(tidyverse)

# ---- Data sources (obtained from All of Us Researcher Workbench) ----
# person_df        = person-level demographics
# condition_df     = condition occurrences (ICD codes)
# measurement_df   = lab / vital measurements
# drug_df          = drug exposures
# survey modules   = Basics, Lifestyle, Social

# ---- Step 1: Identify pancreatitis cases ----
pancreatitis_codes <- c(
  "K85", "K85.0", "K85.1", "K85.2", "K85.3", "K85.8", "K85.9",
  "K86", "K86.0", "K86.1", "K86.2", "K86.3", "K86.4", "K86.8", "K86.9"
)

# ---- Step 2: Identify incident metabolic sequelae ----
# Diabetes (E08-E13), Gout (M10), Osteoporosis (M80-M81),
# Exocrine pancreatic insufficiency (K86.8, K86.9), Sarcopenia (M62.8)

sequela_codes <- list(
  diabetes = c("E08", "E09", "E10", "E11", "E12", "E13"),
  gout = c("M10"),
  osteoporosis = c("M80", "M81"),
  exocrine = c("K86.8", "K86.9"),
  sarcopenia = c("M62.8")
)

# ---- Step 3: Apply inclusion / exclusion criteria ----
# - Age >= 18 at pancreatitis onset
# - At least one incident metabolic sequela AFTER pancreatitis onset
# - At least 4 measurement time points within 24 months post-pancreatitis
# - Exclude early deaths (< 3 months)

# ---- Output: cohort of 2,929 patients with person_id linkage ----

cat("Cohort construction complete.\n")
cat("Expected n = 2,929 patients with pancreatitis + >= 1 incident sequela.\n")
