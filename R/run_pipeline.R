run_treearchi_mlr <- function(data_path,
                              out_base,
                              response_var,
                              forced_var,
                              predictor_vars,
                              random_effect_var = "genus",
                              id_col = "tls_id",
                              analysis_tag = "TREEARCHI_MLR",
                              make_coef_plot = TRUE,
                              make_diagnostics = TRUE,
                              make_local_effect_lineplots = TRUE,
                              label_n = 5) {

  paths <- set_analysis_paths(data_path, out_base)

  run_tag <- paste0(
    safe_filename(analysis_tag, 25), "_",
    safe_filename(response_var, 25), "_FORCED_",
    safe_filename(forced_var, 25), "_",
    format(Sys.time(), "%Y%m%d_%H%M%S")
  )

  out_root <- file.path(paths$out_base, run_tag)
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

  dir_model <- file.path(out_root, "01_model")
  dir_coef  <- file.path(out_root, "02_coefficients")
  dir_local <- file.path(out_root, "03_local_effects")
  dir_diag  <- file.path(out_root, "04_diagnostics")
  dir_lines <- file.path(out_root, "05_local_effect_lineplots")

  dir.create(dir_model, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_coef, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_local, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_diag, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_lines, recursive = TRUE, showWarnings = FALSE)

  dat_raw <- read_messy_csv(paths$data_path)

  needed <- unique(c(response_var, forced_var, predictor_vars, random_effect_var, id_col))
  missing_vars <- setdiff(needed, names(dat_raw))

  if (length(missing_vars) > 0) {
    stop("Missing variables in data: ", paste(missing_vars, collapse = ", "))
  }

  dat0 <- dat_raw[, needed, drop = FALSE]
  dat0 <- tidyr::drop_na(dat0)

  rand_terms <- paste0("(1|", random_effect_var, ")")

  dat_model <- prepare_model_data(
    df = dat0,
    response = response_var,
    predictors = predictor_vars
  )

  dat_model[[random_effect_var]] <- dat0[[random_effect_var]]
  dat_model[[id_col]] <- dat0[[id_col]]

  model_res <- run_interaction_quadratic_bic(
    datQ = dat_model,
    response = response_var,
    forced_var = forced_var,
    x_vars = predictor_vars,
    rand_terms = rand_terms
  )

  selected_model <- model_res$selected$model

  saveRDS(selected_model, file.path(dir_model, "model_selected.rds"))
  saveRDS(model_res, file.path(dir_model, "model_result_full.rds"))

  safe_write_csv(
    tibble::tibble(
      selected_terms = model_res$selected$terms
    ),
    file.path(dir_model, "selected_terms.csv")
  )

  safe_write_csv(
    model_res$selected$history,
    file.path(dir_model, "selection_history.csv")
  )

  fixed_effects <- broom.mixed::tidy(
    selected_model,
    effects = "fixed",
    conf.int = TRUE
  )

  safe_write_csv(
    fixed_effects,
    file.path(dir_model, "fixed_effects.csv")
  )

  model_info <- tibble::tibble(
    analysis_tag = analysis_tag,
    response_var = response_var,
    forced_var = forced_var,
    random_effect_var = random_effect_var,
    id_col = id_col,
    n_used = nrow(stats::model.frame(selected_model)),
    n_after_drop_na = nrow(dat_model),
    bic = model_res$selected$BIC,
    predictors_start = paste(predictor_vars, collapse = ", "),
    selected_terms = paste(model_res$selected$terms, collapse = " + "),
    model_type = "LMM with all 2-way interactions and quadratic terms; backward BIC selection",
    transform_note = "Selective log transform for metric variables; centered after transform"
  )

  safe_write_csv(
    model_info,
    file.path(dir_model, "model_info.csv")
  )

  safe_write_lines(
    capture.output(summary(selected_model)),
    file.path(dir_model, "model_summary.txt")
  )

  local_effects <- compute_local_effects_table(
    model = selected_model,
    dat = dat_model
  )

  local_summary <- summarise_local_effects(local_effects)

  safe_write_csv(
    local_effects,
    file.path(dir_local, "local_effects_all.csv")
  )

  safe_write_csv(
    local_summary,
    file.path(dir_local, "local_effects_summary.csv")
  )

  local_effect_lineplots <- NULL

  if (isTRUE(make_local_effect_lineplots)) {
    local_effect_lineplots <- save_local_effect_lineplots(
      local_effects = local_effects,
      out_dir = dir_lines,
      response_label = label_term(response_var)
    )
  }

  diagnostics <- NULL

  if (isTRUE(make_diagnostics)) {
    diagnostics <- make_diag_plots_lmm(
      model = selected_model,
      dat = dat_model,
      response_label = label_term(response_var),
      out_dir = dir_diag,
      id_col = id_col,
      label_n = label_n
    )
  }

  coef_plot <- NULL

  if (isTRUE(make_coef_plot)) {
    if (exists("COL", mode = "list") && exists("theme_modern", mode = "function")) {
      coef_plot <- coef_plot_lmm(
        model = selected_model,
        model_name = analysis_tag,
        out_dir = dir_coef,
        labels_map = pretty_names,
        dat = dat_model,
        response = response_var,
        local_effect_summary = local_summary
      )
    } else {
      message("Skipping coef plot: COL/theme_modern not found. Add R/theme_helpers.R first.")
    }
  }

  safe_write_csv(
    tibble::tibble(
      output_root = out_root,
      model_dir = dir_model,
      coef_dir = dir_coef,
      local_effect_dir = dir_local,
      diagnostics_dir = dir_diag,
      local_effect_lineplots_dir = dir_lines
    ),
    file.path(out_root, "run_paths.csv")
  )

  list(
    output_root = out_root,
    data_used = dat_model,
    model_result = model_res,
    selected_model = selected_model,
    fixed_effects = fixed_effects,
    local_effects = local_effects,
    local_summary = local_summary,
    local_effect_lineplots = local_effect_lineplots,
    diagnostics = diagnostics,
    coef_plot = coef_plot,
    model_info = model_info
  )
}
