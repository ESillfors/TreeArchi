# TreeArchi

TreeArchi is an R package for analysing tree architectural trait relationships using size-controlled linear mixed-effects models.

The package provides a reproducible workflow for model fitting, automatic backward BIC model selection, nonlinear and interaction modelling, local marginal effect estimation, diagnostic evaluation, and publication-ready visualisations.

TreeArchi was developed for terrestrial laser scanning (TLS)-derived tree architectural traits, but the workflow can be applied to other datasets containing continuous structural variables.

---


## Features

- Size-controlled linear mixed-effects models
- Automatic backward BIC model selection
- Quadratic effects
- Two-way interaction effects
- Local marginal effect estimation
- Coefficient plots
- Extended model diagnostics
- Local effect lineplots
- Group-specific analyses
- Height–DBH shift grouping
- Sensitivity analyses
- Organized output folders

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("ESillfors/TreeArchi")
```

---

## Basic workflow

```r
library(TreeArchi)

data_path <- "path/to/your/data.csv"
out_base  <- "path/to/output/folder"

res <- run_treearchi_mlr(
  data_path = data_path,
  out_base = out_base,
  response_var = "projected_area_m2",
  forced_var = "tree_vol_m3",
  predictor_vars = c(
    "tree_height_m",
    "alpha_volume_m3",
    "csh_raw",
    "branch_len",
    "tree_vol_m3"
  ),
  random_effect_var = "genus",
  id_col = "tls_id",
  analysis_tag = "BASIC_MLR"
)
```

The function creates a timestamped output folder containing model summaries, coefficient plots, local marginal effect tables, diagnostic plots, and local effect visualisations.

---

## Group-specific analysis

TreeArchi can also run group-specific mixed-effects models, for example across plots, treatment groups, or height–DBH shift groups.

```r
run_treearchi_group_mlr(
  data_path = data_path,
  out_base = out_base,
  response_var = "projected_area_m2",
  group_var = "plot",
  predictor_vars = c(
    "tree_height_m",
    "alpha_volume_m3",
    "csh_raw",
    "branch_len",
    "tree_vol_m3"
  ),
  group_reference = NULL,
  random_effect_var = "genus",
  id_col = "tls_id",
  quadratic_vars = c(
    "tree_height_m",
    "alpha_volume_m3",
    "tree_vol_m3"
  ),
  interaction_pairs = list(
    c("tree_height_m", "alpha_volume_m3"),
    c("tree_height_m", "tree_vol_m3"),
    c("alpha_volume_m3", "tree_vol_m3")
  ),
  analysis_tag = "GROUP_MLR",
  run_diagnostics = TRUE
)
```

---

## Height–DBH shift groups

TreeArchi can classify trees according to whether they are shorter or taller than expected for their DBH.

```r
dat <- read.csv("path/to/your/data.csv")

dat <- make_treearchi_shift_groups(
  dat,
  height_var = "tree_height_m",
  dbh_var = "dbh_m",
  mode = "binary"
)

table(dat$height_dbh_shift_group, useNA = "ifany")
```

Available modes:

- `"binary"`: creates two groups, `shorter` and `taller`
- `"three_class"`: creates three groups, `shorter`, `expected`, and `taller`
- `"extreme"`: creates two groups, `shorter` and `taller`, while trees close to the expected height–DBH relationship are set to `NA`

The resulting group variable can be used in `run_treearchi_group_mlr()`.

---

## Sensitivity analysis

TreeArchi includes a sensitivity framework for evaluating how trimming selected observations influences local marginal effects and interaction patterns.

```r
run_treearchi_sensitivity(
  data_path = data_path,
  out_base = out_base,
  response_var = "projected_area_m2",
  group_var = "plot",
  predictor_vars = c(
    "tree_height_m",
    "alpha_volume_m3",
    "csh_raw",
    "branch_len",
    "tree_vol_m3"
  ),
  trim_vars = c("branch_len", "csh_raw"),
  trim_props = c(0, 0.10, 0.20, 0.30),
  trim_tails = c("lower", "upper"),
  random_effect_var = "genus",
  id_col = "tls_id",
  analysis_tag = "SENSITIVITY"
)
```

---

## Output structure

Each `run_treearchi_mlr()` analysis creates a timestamped output folder with the following structure:

```text
01_model/
02_coefficients/
03_local_effects/
04_diagnostics/
05_local_effect_lineplots/
```

| Folder | Contents |
|---|---|
| `01_model/` | Selected model, fixed-effect estimates, model summary, BIC selection history |
| `02_coefficients/` | Coefficient plot and coefficient table |
| `03_local_effects/` | Local marginal effect tables |
| `04_diagnostics/` | Residual diagnostics, Q–Q plot, leverage plot, diagnostic values |
| `05_local_effect_lineplots/` | Local marginal effect visualisations |

---

## Main functions

```r
run_treearchi_mlr()
run_treearchi_group_mlr()
run_treearchi_sensitivity()
make_treearchi_shift_groups()
collect_treearchi_sensitivity_effects()
select_treearchi_sensitivity_top_cases()
select_treearchi_convergence_cases()
```

---

## Status

TreeArchi is under active development and is currently intended for research workflows involving TLS-derived tree architectural traits, structural ecology, and size-controlled trait relationship modelling.

---

## Citation

If you use TreeArchi in scientific work, please cite the accompanying manuscript once available.
