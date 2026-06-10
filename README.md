# TreeArchi

TreeArchi is an R package for exploring tree architectural trait networks using size-controlled multiple linear regression and mixed-effects models. The package is designed to identify nonlinear, interaction-dependent, and group-specific relationships among architectural traits.

## Current features

- Size-controlled multiple linear regression
- Mixed-effects models with genus-level random effects
- Group-specific MLR analyses
- Height-DBH shift grouping
- Local effect tables
- Local effect lineplots
- Diagnostic plots
- Organized output folders

## Installation

```r
# install.packages("devtools")
devtools::install_github("YOUR_USERNAME/TreeArchi")
```

## Basic workflow

```r
library(TreeArchi)

data_path <- "path/to/your/data.csv"
out_base <- "path/to/output/folder"

run_treearchi_mlr(
  data_path = data_path,
  out_base = out_base,
  response_var = "projected_area_m2",
  predictor_vars = c(
    "tree_height_m",
    "alpha_volume_m3",
    "csh_raw",
    "branch_len",
    "tree_vol_m3"
  ),
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
  analysis_tag = "BASIC_MLR"
)
```

## Group-specific analysis

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

## Height-DBH shift groups

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
- `"extreme"`: creates two groups, `shorter` and `taller`, while trees close to the expected height-DBH relationship are set to `NA`

Example:

```r
dat <- make_treearchi_shift_groups(
  dat,
  height_var = "tree_height_m",
  dbh_var = "dbh_m",
  mode = "binary"
)

shift_data_path <- "path/to/data_with_shift_groups.csv"

write.csv(
  dat,
  shift_data_path,
  row.names = FALSE
)

run_treearchi_group_mlr(
  data_path = shift_data_path,
  out_base = "path/to/output/folder",
  response_var = "projected_area_m2",
  group_var = "height_dbh_shift_group",
  predictor_vars = c(
    "tree_height_m",
    "alpha_volume_m3",
    "csh_raw",
    "branch_len",
    "tree_vol_m3"
  ),
  group_reference = "shorter",
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
  analysis_tag = "HEIGHT_DBH_SHIFT_BINARY",
  run_diagnostics = TRUE
)
```

## Output structure

TreeArchi creates organized output folders:

```text
01_model
02_coefficients
03_local_effects
04_diagnostics
05_local_effect_lineplots
```

These folders contain model summaries, coefficient tables, local effect tables, diagnostic plots, and local effect visualizations.

## Notes

TreeArchi is currently under active development. The package is intended for research workflows involving tree architectural traits, TLS-derived structural variables, and size-controlled trait network analyses.
