#==============================================================================
# 07 - Subgroup analyses
# Stratified by sex, age, BMI; interaction tests
#==============================================================================

library(tidyverse)
library(survival)

# ---- Data setup ----
sub_data <- train_assignments %>%
  left_join(patient_data, by = "person_id") %>%
  mutate(
    trajectory_group = factor(trajectory_group),
    death_event = has_death,
    surv_time = ifelse(has_death, death_time, censor_time)
  ) %>%
  filter(!is.na(surv_time), surv_time > 0)

# ---- Subgroup 1: Sex ----
by_sex <- list()
for (sx in c("male", "female")) {
  sub <- sub_data %>% filter(gender == sx)
  if (nrow(sub) < 50 || sum(sub$death_event) < 5) next
  fit <- coxph(Surv(surv_time, death_event) ~ trajectory_group + age_at_pancreatitis + bmi,
               data = sub)
  by_sex[[sx]] <- summary(fit)
}

# ---- Subgroup 2: Age (< 60 vs >= 60) ----
sub_data <- sub_data %>%
  mutate(age_group = ifelse(age_at_pancreatitis < 60, "<60", ">=60"))

by_age <- list()
for (ag in unique(sub_data$age_group)) {
  sub <- sub_data %>% filter(age_group == ag)
  if (nrow(sub) < 50 || sum(sub$death_event) < 5) next
  fit <- coxph(Surv(surv_time, death_event) ~ trajectory_group + gender + bmi,
               data = sub)
  by_age[[ag]] <- summary(fit)
}

# ---- Subgroup 3: BMI (< 24 vs >= 24 kg/m^2) ----
sub_data <- sub_data %>%
  mutate(bmi_group = ifelse(bmi < 24, "<24", ">=24"))

by_bmi <- list()
for (bg in unique(sub_data$bmi_group)) {
  sub <- sub_data %>% filter(bmi_group == bg)
  if (nrow(sub) < 50 || sum(sub$death_event) < 5) next
  fit <- coxph(Surv(surv_time, death_event) ~ trajectory_group + age_at_pancreatitis + gender,
               data = sub)
  by_bmi[[bg]] <- summary(fit)
}

# ---- Interaction tests ----
# Trajectory group x sex
fit_int_sex <- coxph(Surv(surv_time, death_event) ~ trajectory_group * gender +
                       age_at_pancreatitis + bmi, data = sub_data)
anova_int_sex <- anova(fit_int_sex)

# Trajectory group x age (continuous)
fit_int_age <- coxph(Surv(surv_time, death_event) ~ trajectory_group * age_at_pancreatitis +
                       gender + bmi, data = sub_data)
anova_int_age <- anova(fit_int_age)

# Trajectory group x BMI
fit_int_bmi <- coxph(Surv(surv_time, death_event) ~ trajectory_group * bmi +
                       age_at_pancreatitis + gender, data = sub_data)
anova_int_bmi <- anova(fit_int_bmi)

# ---- Output: forest plot of subgroup HRs ----

cat("Subgroup analyses complete.\n")
cat("Interaction tests: sex, age, BMI\n")
