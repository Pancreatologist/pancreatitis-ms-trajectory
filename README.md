# Metabolic Trajectories and All-Cause Mortality in Post-Pancreatitis Metabolic Sequelae

Analytic code for the study **"Association of metabolic trajectories with all-cause mortality among patients with post-pancreatitis metabolic sequelae: a longitudinal cohort study."**

This repository contains the R code used to identify longitudinal metabolic subphenotypes among patients with post-pancreatitis metabolic sequelae (PPMS) and to quantify their association with all-cause mortality, using data from the [All of Us Research Program](https://www.researchallofus.org).

---

## Overview

Pancreatitis frequently gives rise to a cluster of **post-pancreatitis metabolic sequelae (PPMS)** — diabetes, gout, osteoporosis, exocrine pancreatic insufficiency, and sarcopenia — that extend well beyond diabetes alone. Prior work has largely treated these as fixed, isolated endpoints. This study instead asks whether their *coordinated longitudinal trajectory* defines distinct subphenotypes with different prognoses.

Using **group-based multi-trajectory modeling** of three routinely collected biomarkers — serum albumin, hematocrit, and body mass index (BMI) — measured over the 24 months after pancreatitis onset (n = 2,929), we identified three metabolic subphenotypes and assessed their link to all-cause mortality with sequentially adjusted Cox models, plus extensive sensitivity, landmark, and restricted-cubic-spline analyses.

**Key finding:** the *Late Metabolic Decline* subphenotype — low BMI with persistently low hematocrit — carried a roughly three-fold higher risk of all-cause mortality than the *Metabolic Rebound* subphenotype, independent of age, sex, BMI, heart rate, and GLP-1 receptor agonist use. The prognostic signal resided in the overall recovery trajectory rather than in any single incident diagnosis.

---

## Study summary

| Item | Detail |
|---|---|
| **Design** | Retrospective longitudinal cohort study |
| **Data source** | *All of Us* Research Program (Controlled Tier, v8) |
| **Cohort** | 2,929 adults (≥18 y) with pancreatitis + ≥1 incident metabolic sequela |
| **Exposure** | Trajectory-defined metabolic subphenotype (T1–T3) |
| **Primary outcome** | All-cause mortality |
| **Secondary outcomes** | Incident diabetes, gout, osteoporosis, exocrine pancreatic insufficiency |
| **Software** | R 4.5.3 |

### The three subphenotypes

| Subphenotype | Label | n (%) | Trajectory pattern | Mortality |
|---|---|---|---|---|
| **T1** | Late Metabolic Decline | 963 (32.9%) | Lowest BMI; hematocrit stays low without recovery | **5.2%** |
| **T2** | Anemia-Recovery | 1,637 (55.9%) | Preserved hematocrit; highest BMI | 2.8% |
| **T3** | Metabolic Rebound *(reference)* | 329 (11.2%) | Hematocrit rises the most; broad normalization | 1.5% |

### Headline results

- **T1 vs T3 mortality** — adjusted HR **3.18** (95% CI 1.13–8.97; *P* = 0.029) in the fully adjusted model; stable across all five sequential models (HR 3.18–3.60).
- **Timing** — survival curves diverged by **36 months** (log-rank *P* = 0.001) but not at 6 or 12 months, indicating separation emerges over the medium term.
- **Dose-response** — a composite metabolic score showed a **linear** dose-response with mortality (*P* < 0.001; non-linearity *P* > 0.8) under both Cox-weighted and equal-weighted derivations.
- **Robustness** — the association held across median / kNN / MICE imputation, cohort restrictions, and landmark analyses at 6 and 12 months.
- **Not diagnosis-driven** — subphenotype was **not** significantly associated with any individual incident sequela, so the prognostic signal lies in the coordinated trajectory, not a single downstream diagnosis.

---

## Analysis pipeline

The code reproduces the following analytical steps described in the manuscript:

1. **Cohort construction** — identify adults with acute/chronic pancreatitis and ≥1 incident metabolic sequela (ICD codes: diabetes E08–E13, gout M10, osteoporosis M80–M81, exocrine pancreatic insufficiency K86.8/K86.9, sarcopenia M62.8); apply sequential exclusions to isolate *new-onset* sequelae and remove early (<3-month) deaths.
2. **Baseline / covariate extraction** — demographics, BMI, heart rate, comorbidities, medication exposure (insulin, metformin, GLP-1 RA, SGLT2i, statins), and labs within a ±180-day window; linkage of Social / Lifestyle / Basics survey modules and Fitbit wearable metrics.
3. **Group-based multi-trajectory modeling** — joint modeling of albumin, hematocrit, and BMI; model selection over 2–5 groups and polynomial degrees 1–3 by BIC (final: quadratic, 3 groups); posterior-probability class assignment.
4. **Baseline characteristics** — ANOVA / Kruskal-Wallis and χ² / Fisher tests across subphenotypes.
5. **Cox proportional hazards** — five sequentially adjusted models (unadjusted → + age/sex → + BMI → + heart rate → + GLP-1 RA), with T3 as reference.
6. **Kaplan-Meier** — survival compared at 6, 12, and 36 months with window-specific right-censoring and log-rank tests.
7. **Subgroup / interaction analyses** — stratified by sex, age (<60 vs ≥60), and BMI (<24 vs ≥24 kg/m²).
8. **Sensitivity analyses** — median, kNN, and MICE imputation; measurement-count restrictions; extended clearance window.
9. **Restricted cubic splines** — composite metabolic score (Cox-weighted and equal-weighted) modeled with 4-knot RCS; nonlinearity tested by likelihood-ratio test.
10. **Landmark analyses** — landmarks at 6, 12, and 24 months to address immortal-time bias.
11. **Component outcomes** — separate Cox models for each incident sequela.
12. **Drug–mortality analyses** — subphenotype-stratified Cox models across time windows.

---

## Requirements

- **R** ≥ 4.5.3

Key packages used in the manuscript:

```r
install.packages(c(
  "gbmt",      # group-based multi-trajectory modeling
  "rms",       # restricted cubic splines
  "mice",      # multiple imputation by chained equations
  "survival",  # Cox models, Kaplan-Meier
  "survminer"  # survival curve visualization
))
```

> Additional utility packages (e.g. for data wrangling, tables, and figures) may be required depending on the scripts; see the top of each script for its `library()` calls.

---

## Repository structure

```
.
├── README.md
├── LICENSE
├── R/
│   ├── 01_cohort_construction.R         # Cohort identification & exclusions
│   ├── 02_covariate_extraction.R        # Baseline covariates, labs, meds, surveys
│   ├── 03_trajectory_modeling.R         # GBMTM trajectory model (albumin + BMI + Hct)
│   ├── 04_cox_models.R                  # 5 sequential Cox models for mortality
│   ├── 05_km_landmark.R                 # Kaplan-Meier + landmark analyses
│   ├── 06_restricted_cubic_splines.R    # RCS dose-response (Cox-weighted + equal-weighted)
│   ├── 07_subgroup_analyses.R           # Sex / age / BMI subgroups + interaction tests
│   ├── 08_sensitivity_analyses.R        # Imputation, measurement count, window
│   ├── 09_component_outcomes.R          # Individual sequela outcomes
│   ├── 10_drug_mortality.R              # Drug × subphenotype interaction
│   └── 11_visualization.R               # Figure generation
├── figures/                              # Generated figures
└── tables/                               # Generated tables
```

---

## Data availability

The data analyzed in this study were obtained from the *All of Us* Research Program and are available to authorized researchers through the [Researcher Workbench](https://www.researchallofus.org). Under the program's Data Use and Registration Agreement and Data User Code of Conduct, **participant-level data cannot be redistributed**; access is granted by the program to investigators who complete the required registration and training.

This repository therefore contains **analytic code only**, not participant data.

---

## Citation

If you use this code, please cite the associated publication:

> Wu D, Huang Y, Chen C, Chen M, Evans A, Lv Y, Mukherjee R, Xiao J, Cai W, Huang W, Sutton R, Peng J. *Association of metabolic trajectories with all-cause mortality among patients with post-pancreatitis metabolic sequelae: a longitudinal cohort study.*

*(Update with the full journal citation, year, volume, and DOI once available.)*

---

## Funding

Supported by the China Postdoctoral Foundation (2025M772297), the Jiangsu Funding Program for Excellent Postdoctoral Talent, and the Young Scholars Fostering Fund of the First Affiliated Hospital of Nanjing Medical University (PY2025026). The *All of Us* Research Program is funded by the National Institutes of Health. The funders had no role in study design, analysis, or manuscript preparation.

---

## Acknowledgements

We gratefully acknowledge the *All of Us* Research Program participants and the National Institutes of Health's *All of Us* Research Program for making the participant data available.

---

## License

MIT License — see [LICENSE](LICENSE) for details.
