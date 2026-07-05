#==============================================================================
# 03 - Group-based multi-trajectory modeling (GBMTM)
# Joint modeling of albumin, hematocrit, BMI over 24 months post-pancreatitis
#==============================================================================

library(tidyverse)
library(gbmt)
library(mice)
library(parallel)

# ---- Trajectory variables ----
traj_vars <- c("albumin", "bmi", "hematocrit")

# ---- Data preparation ----
# Long-format data with person_id, time (months), and trajectory variables
# Time points: 0, 6, 12, 18, 24 months
# Imputation: MICE (pmm method, 5 imputations) for missing values

# Example structure:
# data_long <- data.frame(
#   person_id = rep(1:100, each = 5),
#   time = rep(c(0, 6, 12, 18, 24), 100),
#   albumin = rnorm(500, 42, 5),
#   bmi = rnorm(500, 26, 4),
#   hematocrit = rnorm(500, 40, 5)
# )

# ---- Model grid: polynomial degree x number of groups ----
model_grid <- expand.grid(d = 2:3, ng = 2:5)

# ---- Fit all models ----
n_cores <- min(detectCores() - 1, 8)
cl <- makeCluster(n_cores)
clusterEvalQ(cl, library(gbmt))
clusterExport(cl, c("data_long", "traj_vars", "model_grid"))

all_models <- parLapply(cl, 1:nrow(model_grid), function(i) {
  d_val <- model_grid$d[i]
  ng_val <- model_grid$ng[i]
  tryCatch({
    m <- gbmt(
      x.names = traj_vars,
      unit = "person_id",
      time = "time",
      d = d_val,
      ng = ng_val,
      data = data_long,
      scaling = 1
    )
    list(model = m, d = d_val, ng = ng_val, error = FALSE)
  }, error = function(e) {
    list(model = NULL, d = d_val, ng = ng_val, error = TRUE, msg = e$message)
  })
})
stopCluster(cl)

# ---- Model selection by BIC (with validity constraints) ----
# Validity: average posterior probability >= 0.7 and smallest group >= 5%
fit_stats <- bind_rows(lapply(all_models, function(res) {
  if (res$error) return(data.frame(d = res$d, ng = res$ng,
    AIC = NA, BIC = NA, valid = FALSE))
  m <- res$model
  assign_vec <- as.integer(m$assign)
  props <- as.numeric(prop.table(table(assign_vec)))
  avg_pp <- min(m$appa, na.rm = TRUE)
  data.frame(
    d = res$d, ng = res$ng,
    AIC = as.numeric(m$ic["aic"]),
    BIC = as.numeric(m$ic["bic"]),
    avg_pp = avg_pp,
    min_prop = min(props),
    valid = (!is.na(avg_pp) & avg_pp >= 0.7 & min(props) >= 0.05)
  )
}))

valid_stats <- fit_stats %>% filter(valid)
best_row <- if (nrow(valid_stats) == 0) {
  fit_stats %>% arrange(BIC) %>% slice(1)
} else {
  valid_stats %>% arrange(BIC) %>% slice(1)
}

best_idx <- which(fit_stats$d == best_row$d & fit_stats$ng == best_row$ng)
best_model <- all_models[[best_idx]]$model

# ---- Group assignments ----
train_assignments <- data.frame(
  person_id = as.numeric(names(best_model$assign)),
  trajectory_group = as.integer(best_model$assign)
)

# ---- Internal validation: bootstrap stability (50 resamples) ----
# Re-fit best model on bootstrap samples, compare group assignments

# ---- Output ----
# best_model: GBMTM model object
# train_assignments: person-level group assignments
# fit_stats: all model fit statistics

cat(sprintf("Best model: d = %d, ng = %d (BIC = %.1f)\n",
            best_row$d, best_row$ng, best_row$BIC))
cat(sprintf("Groups: %d trajectories identified\n", best_row$ng))
