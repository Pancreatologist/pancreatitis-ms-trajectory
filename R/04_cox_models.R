#==============================================================================
# 04 - Cox proportional hazards models for all-cause mortality
# Five sequentially adjusted models with trajectory subphenotype as exposure
#==============================================================================

library(tidyverse)
library(survival)
library(survminer)

# ---- Data setup ----
# patient_data: one row per patient with death_date, censor_time, covariates
# train_assignments: person_id + trajectory_group

cox_data <- train_assignments %>%
  left_join(patient_data, by = "person_id") %>%
  mutate(
    trajectory_group = factor(trajectory_group),
    death_event = has_death,
    surv_time = ifelse(has_death, death_time, censor_time)
  ) %>%
  filter(!is.na(surv_time), surv_time > 0)

# Reference group = lowest risk group (T3: Metabolic Rebound)
risk_by_group <- cox_data %>%
  group_by(trajectory_group) %>%
  summarise(risk = mean(death_event)) %>%
  arrange(risk)
ref_group <- as.character(risk_by_group$trajectory_group[1])

# ---- Run sequential Cox models ----
run_cox_model <- function(data, covars, model_label, ref_grp) {
  cox_vars <- intersect(covars, names(data))
  cox_ready <- data %>%
    select(person_id, surv_time, death_event, all_of(cox_vars)) %>%
    filter(complete.cases(.)) %>%
    droplevels()

  cox_ready$trajectory_group <- relevel(factor(cox_ready$trajectory_group), ref = ref_grp)
  cox_formula <- as.formula(paste0("Surv(surv_time, death_event) ~ ",
                                    paste(cox_vars, collapse = " + ")))
  cox_model <- coxph(cox_formula, data = cox_ready)
  cox_summary <- summary(cox_model)

  hr_df <- data.frame(
    model = model_label,
    variable = rownames(cox_summary$conf.int),
    HR = cox_summary$conf.int[, "exp(coef)"],
    HR_low = cox_summary$conf.int[, "lower .95"],
    HR_high = cox_summary$conf.int[, "upper .95"],
    p_value = cox_summary$coefficients[, "Pr(>|z|)"],
    n = nrow(cox_ready),
    n_event = sum(cox_ready$death_event),
    stringsAsFactors = FALSE
  )

  ph_test <- cox.zph(cox_model)
  list(hr = hr_df, ph_p = ph_test$table["GLOBAL", "p"], model = cox_model)
}

# Model 1: Unadjusted
m1 <- run_cox_model(cox_data, "trajectory_group", "Model1", ref_group)

# Model 2: + age, sex
m2 <- run_cox_model(cox_data, c("trajectory_group", "age_at_pancreatitis", "gender"),
                    "Model2", ref_group)

# Model 3: + BMI
m3 <- run_cox_model(cox_data, c("trajectory_group", "age_at_pancreatitis", "gender", "bmi"),
                    "Model3", ref_group)

# Model 4: + heart rate
m4 <- run_cox_model(cox_data,
                    c("trajectory_group", "age_at_pancreatitis", "gender", "bmi", "heart_rate"),
                    "Model4", ref_group)

# Model 5: + GLP-1 RA use
m5 <- run_cox_model(cox_data,
                    c("trajectory_group", "age_at_pancreatitis", "gender", "bmi", "heart_rate", "glp1_use"),
                    "Model5", ref_group)

# ---- Combine results ----
all_hr_results <- bind_rows(m1$hr, m2$hr, m3$hr, m4$hr, m5$hr)

# ---- Output ----
# all_hr_results: HR table for all 5 models
# Forest plot saved to figures/

cat("Cox models complete. 5 sequential models fitted.\n")
cat("Primary comparison: T1 (Late Metabolic Decline) vs T3 (Metabolic Rebound)\n")
