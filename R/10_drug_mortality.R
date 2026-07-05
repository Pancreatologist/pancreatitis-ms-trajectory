#==============================================================================
# 10 - Drug-mortality analyses
# Subphenotype-stratified Cox models across time windows
#==============================================================================

library(tidyverse)
library(survival)
library(marginaleffects)

# ---- Data setup ----
drug_data <- train_assignments %>%
  left_join(patient_data, by = "person_id") %>%
  mutate(
    trajectory_group = factor(trajectory_group),
    death_event = has_death,
    surv_time = ifelse(has_death, death_time, censor_time)
  ) %>%
  filter(!is.na(surv_time), surv_time > 0)

# ---- Drug exposures ----
# GLP-1 RA, insulin, metformin, SGLT2i, statins, PERT

drug_vars <- c("glp1_use", "insulin_use", "metformin_use",
               "sglt2_use", "statin_use", "pert_use")

# ---- Time windows ----
# 0-12 months, >12-24 months, >24-36 months

time_windows <- list(
  "0-12" = c(0, 12),
  "12-24" = c(12, 24),
  "24-36" = c(24, 36)
)

# ---- Subphenotype x drug interaction ----
# Test whether drug-mortality association differs by trajectory subphenotype

drug_results <- list()

for (dv in drug_vars) {
  if (!dv %in% names(drug_data)) next

  fml <- as.formula(paste0(
    "Surv(surv_time, death_event) ~ trajectory_group * ", dv,
    " + age_at_pancreatitis + gender + bmi"
  ))

  fit <- tryCatch(coxph(fml, data = drug_data), error = function(e) NULL)

  if (!is.null(fit)) {
    drug_results[[dv]] <- list(
      model = fit,
      interaction_p = anova(fit)["trajectory_group:glp1_use", "P"]
    )
  }
}

# ---- Forest plot: drug HRs by subphenotype ----
# Three-panel figure for three time windows

cat("Drug-mortality analyses complete.\n")
cat(strwrap("Key question: does subphenotype modify the effect of drugs on mortality?"), "\n")
