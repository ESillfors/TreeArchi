#' Make height-DBH shift groups
#'
#' Creates height-DBH deviation groups based on residuals from a height-DBH model.
#'
#' @param dat Data frame.
#' @param height_var Height variable name.
#' @param dbh_var DBH variable name.
#' @param group_var Name of the new grouping column.
#' @param mode Grouping mode: "binary", "three_class", or "extreme".
#' @param expected_sd Width of the expected group in residual SD units.
#' @param min_n Minimum number of complete observations.
#'
#' @return Data frame with a new grouping column.
#' @export
make_treearchi_shift_groups <- function(dat,
                                        height_var = "tree_height_m",
                                        dbh_var = "dbh_m",
                                        group_var = "height_dbh_shift_group",
                                        mode = "binary",
                                        expected_sd = 0.5,
                                        min_n = 10) {

  mode <- match.arg(mode, c("binary", "three_class", "extreme"))

  if (!is.data.frame(dat)) {
    stop("dat must be a data frame.", call. = FALSE)
  }

  required_vars <- c(height_var, dbh_var)
  missing_vars <- setdiff(required_vars, names(dat))

  if (length(missing_vars) > 0) {
    stop(
      "Missing required variable(s): ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  dat[[height_var]] <- as.numeric(gsub(",", ".", dat[[height_var]]))
  dat[[dbh_var]] <- as.numeric(gsub(",", ".", dat[[dbh_var]]))

  complete_idx <- stats::complete.cases(dat[, required_vars, drop = FALSE])

  if (sum(complete_idx) < min_n) {
    stop("Too few complete observations.", call. = FALSE)
  }

  model_dat <- dat[complete_idx, required_vars, drop = FALSE]

  fit <- stats::lm(
    stats::as.formula(paste(height_var, "~", dbh_var)),
    data = model_dat
  )

  res <- stats::residuals(fit)
  res_sd <- stats::sd(res, na.rm = TRUE)

  group_values <- rep(NA_character_, nrow(dat))

  if (mode == "binary") {
    group_values[complete_idx] <- ifelse(res < 0, "shorter", "taller")
  }

  if (mode == "three_class") {
    group_values[complete_idx] <- ifelse(
      res < -expected_sd * res_sd,
      "shorter",
      ifelse(res > expected_sd * res_sd, "taller", "expected")
    )
  }

  if (mode == "extreme") {
    group_values[complete_idx] <- ifelse(
      res < -expected_sd * res_sd,
      "shorter",
      ifelse(res > expected_sd * res_sd, "taller", NA_character_)
    )
  }

  if (mode == "three_class") {
    dat[[group_var]] <- factor(group_values, levels = c("shorter", "expected", "taller"))
  } else {
    dat[[group_var]] <- factor(group_values, levels = c("shorter", "taller"))
  }

  dat
}
