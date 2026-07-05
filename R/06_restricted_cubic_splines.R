#==============================================================================
# 06 - Restricted cubic splines (RCS) for dose-response analysis
# Composite metabolic score (Cox-weighted and equal-weighted) vs mortality
#==============================================================================

library(tidyverse)
library(survival)
library(rms)
library(patchwork)

# ---- Data setup ----
rcs_data <- train_assignments %>%
  left_join(patient_data, by = "person_id") %>%
  mutate(
    death_event = has_death,
    surv_time = ifelse(has_death, death_time, censor_time)
  ) %>%
  filter(!is.na(surv_time), surv_time > 0) %>%
  filter(!is.na(albumin), !is.na(bmi), !is.na(hematocrit))

# ---- Build composite metabolic scores ----

# Score 1: Cox-weighted (linear predictor from Cox model of 3 variables)
cox_base <- coxph(Surv(surv_time, death_event) ~ albumin + bmi + hematocrit,
                  data = rcs_data)
rcs_data$composite_cox <- predict(cox_base, type = "lp")

# Score 2: Z-score equal-weighted (negated so higher = worse = higher risk)
rcs_data <- rcs_data %>%
  mutate(
    composite_z_raw = (-scale(albumin)[, 1]) + (-scale(bmi)[, 1]) + (-scale(hematocrit)[, 1]),
    composite_z = scale(composite_z_raw)[, 1]
  )

# ---- datadist setup for rms ----
dd_dat <- datadist(rcs_data)
options(datadist = "dd_dat")

n_knots <- 4 # 4 knots = 3 spline segments

# ---- RCS: Cox-weighted composite score ----
fit_cox <- cph(Surv(surv_time, death_event) ~ rcs(composite_cox, n_knots) +
                 age_at_pancreatitis + gender,
               data = rcs_data, x = TRUE, y = TRUE)

anova_cox <- anova(fit_cox)
nonlinear_p_cox <- anova_cox[grep("Nonlinear", rownames(anova_cox))[1], "P"]
overall_p_cox <- anova_cox[grep("^composite_cox", rownames(anova_cox))[1], "P"]

pred_cox <- Predict(fit_cox, composite_cox, fun = exp, ref.zero = TRUE)
pred_df_cox <- data.frame(
  x = pred_cox$composite_cox,
  yhat = as.numeric(pred_cox$yhat),
  lower = as.numeric(pred_cox$lower),
  upper = as.numeric(pred_cox$upper)
)

# ---- RCS: Z-score equal-weighted composite score ----
fit_z <- cph(Surv(surv_time, death_event) ~ rcs(composite_z, n_knots) +
               age_at_pancreatitis + gender,
             data = rcs_data, x = TRUE, y = TRUE)

anova_z <- anova(fit_z)
nonlinear_p_z <- anova_z[grep("Nonlinear", rownames(anova_z))[1], "P"]
overall_p_z <- anova_z[grep("^composite_z", rownames(anova_z))[1], "P"]

pred_z <- Predict(fit_z, composite_z, fun = exp, ref.zero = TRUE)
pred_df_z <- data.frame(
  x = pred_z$composite_z,
  yhat = as.numeric(pred_z$yhat),
  lower = as.numeric(pred_z$lower),
  upper = as.numeric(pred_z$upper)
)

# ---- Plotting function ----
make_rcs_plot <- function(df, ref_val, overall_p, nonlinear_p, xlab, title, color) {
  p_label_overall <- ifelse(overall_p < 0.001, "P_overall < 0.001",
                            sprintf("P_overall = %.3f", overall_p))
  p_label_nonlinear <- ifelse(nonlinear_p < 0.001, "P_nonlinear < 0.001",
                              sprintf("P_nonlinear = %.3f", nonlinear_p))

  ggplot(df, aes(x = x)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = color) +
    geom_line(aes(y = yhat), color = color, linewidth = 1.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
    geom_vline(xintercept = ref_val, linetype = "dotted", color = "grey40", linewidth = 0.8) +
    annotate("text", x = min(df$x) + 0.05 * diff(range(df$x)),
             y = max(df$upper) * 0.95, label = p_label_overall,
             hjust = 0, size = 6, fontface = "italic") +
    annotate("text", x = min(df$x) + 0.05 * diff(range(df$x)),
             y = max(df$upper) * 0.80, label = p_label_nonlinear,
             hjust = 0, size = 6, fontface = "italic", color = "#D55E00") +
    scale_y_log10() +
    labs(x = xlab, y = "HR (95% CI)", title = title) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14)
    )
}

# ---- Generate plots ----
p_cox <- make_rcs_plot(pred_df_cox, median(rcs_data$composite_cox),
                       overall_p_cox, nonlinear_p_cox,
                       "Composite Score (Cox-weighted)",
                       "RCS: Composite Score & Mortality", "#E41A1C")

p_z <- make_rcs_plot(pred_df_z, median(rcs_data$composite_z),
                     overall_p_z, nonlinear_p_z,
                     "Composite Score (Z-score equal-weighted)",
                     "RCS: Composite Score & Mortality", "#377EB8")

p_combined <- p_cox / p_z

# ---- Save results ----
# ggsave to figures/ directory

rcs_summary <- data.frame(
  method = c("Cox-weighted", "Z-score equal-weighted"),
  overall_P = c(signif(overall_p_cox, 3), signif(overall_p_z, 3)),
  nonlinear_P = c(signif(nonlinear_p_cox, 3), signif(nonlinear_p_z, 3)),
  knots = n_knots
)
# write_csv to tables/ directory

cat("RCS analysis complete.\n")
cat(sprintf("Cox-weighted: overall P = %.3g, nonlinear P = %.3g\n",
            overall_p_cox, nonlinear_p_cox))
cat(sprintf("Z-score equal-weighted: overall P = %.3g, nonlinear P = %.3g\n",
            overall_p_z, nonlinear_p_z))
