`%||%` <- function(a, b) if (!is.null(a)) a else b

pretty_names <- c(
  tls_id = "TLS ID",
  sba_degrees = "Stem Branch angle (°)",
  ba2 = "Second Branch angle (°)",
  bar = "Branch angle ratio",
  cdhr = "Crown spread",
  rvr = "Relative volume ratio",
  dhr = "DBH / height",
  sbd = "Stem–branch distance (m)",
  csh_raw = "Crown start height (m)",
  ch_raw = "Crown depth (m)",
  sbr = "Stem–branch radius (m)",
  sbl = "Stem–branch length (m)",
  clvr_m.2 = "Branch slenderness (m²)",
  branch_len = "Total branch length (m)",
  branch_vol_m3 = "Total branch volume (m³)",
  cd_raw = "Crown diameter (m)",
  projected_area_m2 = "Projected crown area (m²)",
  alpha_volume_m3 = "Alpha volume (m³)",
  dbh_m = "DBH (m)",
  tree_height_m = "Tree height (m)",
  tree_vol_m3 = "Tree volume (m³)",
  trunk_vol_m3 = "Trunk volume (m³)",
  base_vol_0_10 = "Base volume (0–10%)",
  AGB_TLS = "AGB (TLS)",
  genus = "Genus",
  plot = "Plot"
)

core_tree_size_vars <- c(
  "dbh_m",
  "tree_vol_m3",
  "trunk_vol_m3",
  "base_vol_0_10",
  "AGB_TLS",
  "sbr"
)

height_vars <- c("tree_height_m")

crown_size_vars <- c(
  "projected_area_m2",
  "alpha_volume_m3",
  "cd_raw",
  "branch_len",
  "branch_vol_m3",
  "clvr_m.2"
)

mediator_vars <- c(
  "sbl",
  "sbd",
  "csh_raw",
  "ch_raw"
)

shape_vars <- c(
  "sba_degrees",
  "ba2",
  "bar",
  "cdhr",
  "rvr",
  "dhr"
)

angles_vars <- c("sba_degrees", "ba2")
unitless_vars <- c("bar", "rvr", "cdhr")

metric_vars_allow <- c(
  "AGB_TLS",
  "dbh_m", "tree_height_m", "csh_raw", "ch_raw", "sbl", "sbd", "sbr",
  "branch_len", "clvr_m.2", "branch_vol_m3", "cd_raw",
  "projected_area_m2", "alpha_volume_m3",
  "tree_vol_m3", "trunk_vol_m3", "base_vol_0_10"
)
