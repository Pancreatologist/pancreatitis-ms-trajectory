#==============================================================================
# 05 - Kaplan-Meier survival curves & landmark analyses
# Survival at 6, 12, 36 months; landmark at 6, 12, 24 months
#==============================================================================

library(tidyverse)
library(survival)
library(survminer)

# ---- Data setup ----
km_data <- train_assignments %>%
  left_join(patient_data, by = "person_id") %>%
  mutate(
    trajectory_group = factor(trajectory_group),
    death_event = has_death,
    surv_time = ifelse(has_death, death_time, censor_time)
  ) %>%
  filter(!is.na(surv_time), surv_time > 0)

# ---- Kaplan-Meier at 6, 12, and 36 months ----
time_points <- c(6, 12, 36)

for (t_max in time_points) {
  km_sub <- km_data %>%
    mutate(
      surv_time_t = pmin(surv_time, t_max),
      death_event_t = ifelse(surv_time > t_max, 0, death_event)
    )

  fit_km <- survfit(Surv(surv_time_t, death_event_t) ~ trajectory_group, data = km_sub)

  # Log-rank test
  surv_diff <- survdiff(Surv(surv_time_t, death_event_t) ~ trajectory_group, data = km_sub)
  logrank_p <- 1 - pchisq(surv_diff$chisq, df = length(unique(km_sub$trajectory_group)) - 1)

  cat(sprintf("KM at %d months: log-rank P = %.4f, events = %d\n",
              t_max, logrank_p, sum(km_sub$death_event_t)))

  # Plot with ggsurvplot
  g <- ggsurvplot(
    fit_km,
    data = km_sub,
    xlim = c(0, t_max),
    pval = TRUE,
    pval.method = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    palette = c("#E41A1C", "#377EB8", "#4DAF4A"),
    legend = "bottom",
    xlab = "Time (months)",
    ylab = "Survival Probability",
    ggtheme = theme_classic()
  )

  # Save figure
}

# ---- Landmark analyses (to address immortal-time bias) ----
# Landmarks at 6, 12, 24 months
# Only include patients alive and in follow-up at the landmark time

landmark_times <- c(6, 12, 24)
landmark_results <- list()

for (lm_time in landmark_times) {
  lm_data <- km_data %>%
    filter(surv_time > lm_time | death_event) %>%
    mutate(
      surv_time_lm = pmax(0, surv_time - lm_time),
      death_event_lm = ifelse(surv_time <= lm_time, 0, death_event)
    ) %>%
    filter(!is.na(surv_time_lm), surv_time_lm >= 0)

  fit_lm <- coxph(Surv(surv_time_lm, death_event_lm) ~ trajectory_group +
                    age_at_pancreatitis + gender, data = lm_data)

  landmark_results[[as.character(lm_time)]] <- summary(fit_lm)
  cat(sprintf("Landmark %d months: n = %d, events = %d\n",
              lm_time, nrow(lm_data), sum(lm_data$death_event_lm)))
}

cat("KM and landmark analyses complete.\n")
