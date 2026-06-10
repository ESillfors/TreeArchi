# -------------------------
# Selective log transform
# -------------------------
log_transform_vars <- function(df, vars, allow = metric_vars_allow) {

  for (v in vars) {
    if (v %in% allow && v %in% names(df)) {

      # skip non-numeric columns
      if (!is.numeric(df[[v]])) next

      # avoid log(0)
      min_val <- min(df[[v]], na.rm = TRUE)

      if (!is.finite(min_val)) next

      if (min_val <= 0) {
        shift <- abs(min_val) + 1e-6
        df[[v]] <- log(df[[v]] + shift)
      } else {
        df[[v]] <- log(df[[v]])
      }
    }
  }

  return(df)
}

# -------------------------
# Center variables
# -------------------------
center_vars <- function(df, vars) {

  for (v in vars) {
    if (v %in% names(df)) {

      # skip non-numeric columns
      if (!is.numeric(df[[v]])) next

      df[[v]] <- as.numeric(scale(df[[v]], center = TRUE, scale = FALSE))
    }
  }

  return(df)
}

# -------------------------
# Prepare dataset for model
# -------------------------
prepare_model_data <- function(df, response, predictors) {

  vars <- unique(c(response, predictors))
  df_sub <- df[, vars, drop = FALSE]

  # convert character columns that should be numeric
  for (v in vars) {
    if (v %in% names(df_sub) && is.character(df_sub[[v]])) {
      suppressWarnings({
        x_num <- as.numeric(df_sub[[v]])
      })

      if (sum(!is.na(x_num)) > 0.8 * length(x_num)) {
        df_sub[[v]] <- x_num
      }
    }
  }

  # log transform metric variables only
  df_sub <- log_transform_vars(df_sub, vars)

  # center numeric variables
  df_sub <- center_vars(df_sub, vars)

  return(df_sub)
}
