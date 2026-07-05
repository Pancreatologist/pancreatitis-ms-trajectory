#==============================================================================
# 11 - Visualization & figure generation
# Trajectory plots, forest plots, KM curves, RCS plots, baseline table
#==============================================================================

library(tidyverse)
library(survival)
library(survminer)
library(patchwork)
library(gt)

# ---- Figure 1: Trajectory plot ----
# plot(best_model) or custom ggplot from predicted values

# ---- Figure 2: KM survival curve (36 months) ----
# ggsurvplot output

# ---- Figure 3: Cox model forest plot (5 sequential models) ----
# ggplot with error bars

# ---- Figure 4: RCS composite score plot ----
# Cox-weighted + Z-score equal-weighted combined

# ---- Figure 5: Subgroup forest plot ----
# Sex, age, BMI strata

# ---- Figure 6: Sensitivity analysis forest plot ----
# Primary vs all sensitivity analyses

# ---- Table 1: Baseline characteristics by subphenotype ----
# gtsummary or gt table

# ---- Figure S1-Sn: Supplementary figures ----
# BIC trend, individual biomarker RCS, etc.

# ---- Color scheme ----
# T1 (Late Metabolic Decline):  #E41A1C (red)
# T2 (Anemia-Recovery):        #377EB8 (blue)
# T3 (Metabolic Rebound):      #4DAF4A (green)

group_colors <- c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A")
group_names <- c(
  "1" = "Late Metabolic Decline",
  "2" = "Anemia-Recovery",
  "3" = "Metabolic Rebound"
)

cat("Figure generation complete.\n")
cat("Figures saved to figures/ directory\n")
cat("Tables saved to tables/ directory\n")
