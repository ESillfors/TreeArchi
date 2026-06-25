#' Run TreeArchi sensitivity analysis
#'
#' Re-runs group MLR analyses after trimming lower or upper tails of selected
#' variables within selected target groups. After the model runs, the function
#' can automatically collect local effects, select top sensitivity cases, redraw
#' sensitivity lineplots, and optionally rank convergence cases against a
#' reference group profile.
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
#' @param group_reference Reference group used in group MLR models.
#' @param reference_group One or more reference groups used for convergence ranking.
#' @param make_convergence_ranking Logical. If TRUE, run convergence ranking after collecting local effects.
#' @param convergence_top_n Number of top convergence cases per interaction.
#' @param random_effect_var Random effect variable.
#' @param id_col ID column.
#' @param quadratic_vars Quadratic variables.
#' @param interaction_pairs Interaction pairs.
#' @param analysis_tag Analysis tag.
#' @param run_diagnostics Logical.
#' @param output_level Either "focused" or "full". Focused draws only top-case figures.
#' @param top_n Number of top sensitivity cases per target group x interaction.
#' @param create_summary_outputs Logical. If TRUE, collect effects and create top-case tables.
#' @param redraw_lineplots Logical. If TRUE, redraw sensitivity lineplots after collection.
#'
#' @return Invisibly returns a list with run status, top cases, convergence cases, and output paths.
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
    reference_group = NULL,
    make_convergence_ranking = FALSE,
    convergence_top_n = 5,
    random_effect_var = "genus",
    id_col = "tls_id",
    quadratic_vars = NULL,
    interaction_pairs = NULL,
    analysis_tag = "SENSITIVITY",
    run_diagnostics = TRUE,
    output_level = c("focused", "full"),
    top_n = 5,
    create_summary_outputs = TRUE,
    redraw_lineplots = TRUE
) {

  output_level <- match.arg(output_level)

  if (!file.exists(data_path)) {
    stop("File not found: ", data_path, call. = FALSE)
  }

  trim_tails <- match.arg(
    trim_tails,
    choices = c("lower", "upper"),
    several.ok = TRUE
  )

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

  target_groups <- as.character(target_groups)

  if (isTRUE(make_convergence_ranking) && is.null(reference_group)) {
    warning(
      "make_convergence_ranking = TRUE, but reference_group is NULL. ",
      "Convergence ranking will be skipped.",
      call. = FALSE
    )
  }

  run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  out_root <- file.path(
    out_base,
    paste0(analysis_tag, "_", run_timestamp)
  )

  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

  dir_raw_runs <- file.path(out_root, "raw_model_runs")
  dir_lineplots <- file.path(out_root, "lineplots")
  dir_convergence <- file.path(out_root, "convergence_ranking")

  dir.create(dir_raw_runs, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_lineplots, recursive = TRUE, showWarnings = FALSE)

  run_info_file <- file.path(out_root, "RUN_INFO.txt")

  base::writeLines(
    c(
      "TreeArchi sensitivity analysis",
      paste0("Created: ", Sys.time()),
      paste0("Data path: ", data_path),
      paste0("Output root: ", out_root),
      paste0("Response: ", response_var),
      paste0("Group variable: ", group_var),
      paste0("Target groups: ", paste(target_groups, collapse = ", ")),
      paste0("Predictors: ", paste(predictor_vars, collapse = ", ")),
      paste0("Trim variables: ", paste(trim_vars, collapse = ", ")),
      paste0("Trim proportions: ", paste(trim_props, collapse = ", ")),
      paste0("Trim tails: ", paste(trim_tails, collapse = ", ")),
      paste0("Group reference: ", ifelse(is.null(group_reference), "none", group_reference)),
      paste0("Convergence ranking: ", make_convergence_ranking),
      paste0("Convergence reference group(s): ", ifelse(is.null(reference_group), "none", paste(reference_group, collapse = ", "))),
      paste0("Convergence top n: ", convergence_top_n),
      paste0("Random effect variable: ", random_effect_var),
      paste0("ID column: ", id_col),
      paste0("Quadratic vars: ", ifelse(is.null(quadratic_vars), "none", paste(quadratic_vars, collapse = ", "))),
      paste0(
        "Interaction pairs: ",
        ifelse(
          is.null(interaction_pairs) || length(interaction_pairs) == 0,
          "none",
          paste(vapply(interaction_pairs, paste, collapse = ":", FUN.VALUE = character(1)), collapse = ", ")
        )
      ),
      paste0("Output level: ", output_level),
      paste0("Top n: ", top_n),
      paste0("Create summary outputs: ", create_summary_outputs),
      paste0("Redraw lineplots: ", redraw_lineplots),
      "",
      "Output folders:",
      "raw_model_runs",
      "lineplots",
      "convergence_ranking, if make_convergence_ranking = TRUE",
      "Main result files are written directly to the sensitivity output root.",
      "",
      "Output level interpretation:",
      "- focused: all sensitivity models are fitted, but lineplots are drawn only for top-ranked cases.",
      "- full: all sensitivity models are fitted, and lineplots are drawn for all available sensitivity combinations."
    ),
    con = run_info_file
  )

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
              cutoff_value = NA_real_,
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

          cutoff <- NA_real_
          remove_idx <- rep(FALSE, nrow(dat_i))

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
            dir_raw_runs,
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
                run_diagnostics = run_diagnostics,
                make_coef_plot = FALSE,
                make_group_lineplots = FALSE
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

  status_file <- file.path(out_root, "SENSITIVITY_STATUS.csv")

  utils::write.csv(
    status_tbl,
    status_file,
    row.names = FALSE
  )

  effects_tbl <- NULL
  top_cases <- NULL
  convergence_cases <- NULL
  effects_file <- file.path(out_root, "ALL_SENSITIVITY_LOCAL_EFFECTS.csv")
  top_cases_dir <- out_root
  lineplot_dir <- NULL
  convergence_dir <- NULL

  if (isTRUE(create_summary_outputs)) {

    collect_result <- tryCatch(
      {
        collect_treearchi_sensitivity_effects(
          sensitivity_dir = out_root,
          out_file = effects_file
        )
      },
      error = function(e) e
    )

    if (inherits(collect_result, "error")) {
      warning(
        "Could not collect sensitivity effects: ",
        collect_result$message,
        call. = FALSE
      )
    } else {
      effects_tbl <- collect_result

      top_result <- tryCatch(
        {
          select_treearchi_sensitivity_top_cases(
            effects_tbl = effects_tbl,
            out_dir = top_cases_dir,
            top_n = top_n
          )
        },
        error = function(e) e
      )

      if (inherits(top_result, "error")) {
        warning(
          "Could not select top sensitivity cases: ",
          top_result$message,
          call. = FALSE
        )
      } else {
        top_cases <- top_result
      }

      if (isTRUE(make_convergence_ranking) && !is.null(reference_group)) {

        convergence_dir <- dir_convergence

        convergence_result <- tryCatch(
          {
            select_treearchi_convergence_cases(
              effects_tbl = effects_tbl,
              out_dir = convergence_dir,
              target_group = target_groups[1],
              reference_group = reference_group,
              top_n = convergence_top_n
            )
          },
          error = function(e) e
        )

        if (inherits(convergence_result, "error")) {
          warning(
            "Could not create convergence ranking: ",
            convergence_result$message,
            call. = FALSE
          )
          convergence_dir <- NULL
        } else {
          convergence_cases <- convergence_result
        }
      }

      if (isTRUE(redraw_lineplots)) {

        lineplot_dir <- dir_lineplots

        redraw_result <- tryCatch(
          {
            redraw_treearchi_sensitivity_lineplots(
              sensitivity_dir = out_root,
              effects_file = effects_file,
              out_dir = lineplot_dir,
              top_cases = top_cases,
              output_level = output_level
            )
          },
          error = function(e) e
        )

        if (inherits(redraw_result, "error")) {
          warning(
            "Could not redraw sensitivity lineplots: ",
            redraw_result$message,
            call. = FALSE
          )
        } else {
          lineplot_dir <- redraw_result
        }
      }
    }
  }

  base::writeLines(
    c(
      "",
      "Run completed:",
      paste0("Finished: ", Sys.time()),
      paste0("Status file: ", status_file),
      paste0("Effects file: ", effects_file),
      paste0("Top cases directory: ", top_cases_dir),
      paste0("Lineplot directory: ", ifelse(is.null(lineplot_dir), "not created", lineplot_dir)),
      paste0("Convergence directory: ", ifelse(is.null(convergence_dir), "not created", convergence_dir))
    ),
    con = run_info_file,
    sep = "\n",
    useBytes = TRUE
  )

  run_paths <- data.frame(
    output_root = out_root,
    raw_model_runs_dir = dir_raw_runs,
    lineplot_dir = ifelse(is.null(lineplot_dir), NA_character_, lineplot_dir),
    convergence_dir = ifelse(is.null(convergence_dir), NA_character_, convergence_dir),
    stringsAsFactors = FALSE
  )

  utils::write.csv(
    run_paths,
    file.path(out_root, "RUN_PATHS.csv"),
    row.names = FALSE
  )

  message("Sensitivity analysis finished.")
  message("Output folder: ", out_root)

  invisible(list(
    status = status_tbl,
    effects = effects_tbl,
    top_cases = top_cases,
    convergence_cases = convergence_cases,
    output_root = out_root,
    status_file = status_file,
    effects_file = effects_file,
    top_cases_dir = top_cases_dir,
    lineplot_dir = lineplot_dir,
    convergence_dir = convergence_dir,
    run_paths = run_paths
  ))
}
