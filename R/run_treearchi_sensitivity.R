#' Run TreeArchi sensitivity analysis
#'
#' Re-runs group MLR analyses after trimming lower or upper tails of selected
#' variables within selected target groups.
#'
#' @param data_path Path to input CSV data.
#' @param out_base Output directory.
#' @param response_var Response variable.
#' @param group_var Grouping variable.
#' @param predictor_vars Predictor variables.
#' @param trim_vars Variables used for trimming.
#' @param target_groups Groups where trimming is applied.
#' @param trim_props Proportions to trim.
#' @param trim_tails Tails to trim: "lower" and/or "upper".
#' @param group_reference Reference group.
#' @param random_effect_var Random effect variable.
#' @param id_col ID column.
#' @param quadratic_vars Quadratic variables.
#' @param interaction_pairs Interaction pairs.
#' @param analysis_tag Analysis tag.
#' @param run_diagnostics Logical.
#'
#' @return Invisibly returns a data frame with sensitivity run status.
#' @export
run_treearchi_sensitivity <- function(
    data_path,
    out_base,
    response_var,
    group_var,
    predictor_vars,
    trim_vars,
    target_groups = NULL,
    trim_props = c(0, 0.10, 0.20, 0.30),
    trim_tails = c("lower", "upper"),
    group_reference = NULL,
    random_effect_var = "genus",
    id_col = "tls_id",
    quadratic_vars = NULL,
    interaction_pairs = NULL,
    analysis_tag = "SENSITIVITY",
    run_diagnostics = TRUE
) {

  if (!file.exists(data_path)) {
    stop("File not found: ", data_path, call. = FALSE)
  }
  trim_props <- sort(unique(c(0, trim_props)))

  dat <- read_messy_csv(data_path)

  needed_cols <- unique(c(
    response_var,
    group_var,
    predictor_vars,
    trim_vars,
    random_effect_var,
    id_col
  ))

  missing_cols <- setdiff(needed_cols, names(dat))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  dat[[group_var]] <- as.character(dat[[group_var]])

  if (is.null(target_groups)) {
    target_groups <- sort(unique(dat[[group_var]][!is.na(dat[[group_var]])]))
  }

  run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  out_root <- file.path(
    out_base,
    paste0(analysis_tag, "_", run_timestamp)
  )

  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

  status_rows <- list()
  run_i <- 0

  for (target_group_i in target_groups) {
    for (trim_var_i in trim_vars) {
      for (trim_tail_i in trim_tails) {
        for (trim_prop_i in trim_props) {

          run_i <- run_i + 1

          run_id <- paste0(
            "TG", target_group_i,
            "_", trim_var_i,
            "_", substr(trim_tail_i, 1, 1),
            round(trim_prop_i * 100)
          )

          message("[", run_i, "] Running ", run_id)

          dat_i <- dat

          dat_i[[trim_var_i]] <- suppressWarnings(
            as.numeric(gsub(",", ".", dat_i[[trim_var_i]]))
          )

          target_idx <- dat_i[[group_var]] == target_group_i &
            is.finite(dat_i[[trim_var_i]])

          n_before_target <- sum(target_idx, na.rm = TRUE)

          if (n_before_target < 3) {
            status_rows[[length(status_rows) + 1]] <- data.frame(
              run_id = run_id,
              target_group = target_group_i,
              trim_var = trim_var_i,
              trim_tail = trim_tail_i,
              trim_prop = trim_prop_i,
              n_before_total = nrow(dat_i),
              n_after_total = NA_integer_,
              n_before_target = n_before_target,
              n_removed = NA_integer_,
              fit_ok = FALSE,
              message = "Too few finite observations in target group",
              stringsAsFactors = FALSE
            )
            next
          }

          x <- dat_i[[trim_var_i]][target_idx]

          if (trim_tail_i == "lower") {
            cutoff <- stats::quantile(
              x,
              probs = trim_prop_i,
              na.rm = TRUE,
              type = 7
            )

            remove_idx <- target_idx & dat_i[[trim_var_i]] < cutoff
          }

          if (trim_tail_i == "upper") {
            cutoff <- stats::quantile(
              x,
              probs = 1 - trim_prop_i,
              na.rm = TRUE,
              type = 7
            )

            remove_idx <- target_idx & dat_i[[trim_var_i]] > cutoff
          }

          n_removed <- sum(remove_idx, na.rm = TRUE)

          dat_trimmed <- dat_i[!remove_idx, , drop = FALSE]

          run_dir <- file.path(
            out_root,
            safe_filename(run_id, max_chars = 120)
          )

          dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

          trimmed_path <- file.path(run_dir, "DATA_TRIMMED.csv")

          utils::write.csv(
            dat_trimmed,
            trimmed_path,
            row.names = FALSE
          )

          fit_result <- tryCatch(
            {
              run_treearchi_group_mlr(
                data_path = trimmed_path,
                out_base = run_dir,
                response_var = response_var,
                group_var = group_var,
                predictor_vars = predictor_vars,
                group_reference = group_reference,
                random_effect_var = random_effect_var,
                id_col = id_col,
                quadratic_vars = quadratic_vars,
                interaction_pairs = interaction_pairs,
                analysis_tag = "sens",
                run_diagnostics = run_diagnostics
              )

              TRUE
            },
            error = function(e) e
          )

          fit_ok <- isTRUE(fit_result)
          msg <- if (fit_ok) "ok" else fit_result$message

          status_rows[[length(status_rows) + 1]] <- data.frame(
            run_id = run_id,
            target_group = target_group_i,
            trim_var = trim_var_i,
            trim_tail = trim_tail_i,
            trim_prop = trim_prop_i,
            cutoff_value = as.numeric(cutoff),
            n_before_total = nrow(dat_i),
            n_after_total = nrow(dat_trimmed),
            n_before_target = n_before_target,
            n_removed = n_removed,
            fit_ok = fit_ok,
            message = msg,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  status_tbl <- do.call(rbind, status_rows)

  utils::write.csv(
    status_tbl,
    file.path(out_root, "SENSITIVITY_STATUS.csv"),
    row.names = FALSE
  )

  message("Sensitivity analysis finished.")
  message("Output folder: ", out_root)

  invisible(status_tbl)
}
