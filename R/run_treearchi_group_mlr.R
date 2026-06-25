#' Run group-specific TreeArchi MLR workflow
#'
#' @export
run_treearchi_group_mlr <- function(data_path,
                                    out_base,
                                    response_var,
                                    group_var,
                                    predictor_vars,
                                    group_reference = NULL,
                                    random_effect_var = "genus",
                                    id_col = "tls_id",
                                    quadratic_vars = NULL,
                                    interaction_pairs = NULL,
                                    analysis_tag = "GROUP_MLR",
                                    run_diagnostics = TRUE,
                                    make_coef_plot = TRUE,
                                    make_group_lineplots = TRUE) {

  requireNamespace("dplyr")
  requireNamespace("tidyr")
  requireNamespace("purrr")
  requireNamespace("lme4")
  requireNamespace("lmerTest")
  requireNamespace("broom.mixed")
  requireNamespace("performance")
  requireNamespace("readr")
  requireNamespace("ggplot2")
  requireNamespace("ggrepel")
  requireNamespace("patchwork")

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  safe_filename <- function(x, max_chars = 80) {
    x <- as.character(x)
    x <- gsub("[^A-Za-z0-9_\\.-]", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    ifelse(nchar(x) > max_chars, substr(x, 1, max_chars), x)
  }

  safe_write_csv <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(x, path)
  }

  safe_write_lines <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    readr::write_lines(as.character(x), path)
  }

  safe_ggsave <- function(filename, plot, width, height, dpi = 420) {
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(
      filename = filename,
      plot = plot,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white",
      limitsize = FALSE
    )
  }

  read_messy_csv <- function(path) {
    if (!file.exists(path)) stop("File not found: ", path)

    first_line <- readLines(path, n = 1, warn = FALSE)
    sep <- ifelse(
      length(strsplit(first_line, ";", fixed = TRUE)[[1]]) >
        length(strsplit(first_line, ",", fixed = TRUE)[[1]]),
      ";", ","
    )

    read.table(
      path,
      header = TRUE,
      sep = sep,
      stringsAsFactors = FALSE,
      fileEncoding = "UTF-8-BOM",
      check.names = FALSE
    ) |>
      tibble::as_tibble()
  }

  make_group_factor <- function(x, ref = NULL) {
    x <- as.character(x)
    levs <- sort(unique(x[!is.na(x)]))

    if (!is.null(ref)) {
      ref <- as.character(ref)
      if (!(ref %in% levs)) {
        stop("group_reference not found in group_var: ", ref)
      }
      levs <- c(ref, setdiff(levs, ref))
    }

    factor(x, levels = levs)
  }

  strip_units_label <- function(x) {
    x <- gsub("\\s*\\[[23]D\\]", "", x)
    x <- gsub("\\s*\\([^\\)]*\\)", "", x)
    x <- gsub("\\s+", " ", x)
    trimws(x)
  }

  label_term <- function(x) {
    pretty_names <- c(
      tls_id = "TLS ID",
      plot = "Plot",
      plot_factor = "Plot",
      genus = "Genus",
      tree_vol_m3 = "Tree volume (m³)",
      tree_height_m = "Tree height (m)",
      alpha_volume_m3 = "Alpha volume (m³)",
      csh_raw = "Crown start height (m)",
      ch_raw = "Crown depth (m)",
      branch_len = "Total branch length (m)",
      projected_area_m2 = "Projected crown area (m²)",
      dbh_m = "DBH (m)",
      sbd = "Stem–branch distance (m)",
      sba_degrees = "Stem branch angle (°)",
      ba2 = "Second branch angle (°)",
      bar = "Branch angle ratio",
      cdhr = "Crown spread"
    )

    vapply(x, function(xx) {
      if (xx == "alpha_volume_m3") return("Alpha volume (m³) [3D]")
      if (xx == "projected_area_m2") return("Projected crown area (m²) [2D]")
      if (xx %in% names(pretty_names)) return(unname(pretty_names[xx]))
      xx
    }, character(1))
  }

  metric_vars_allow <- c(
    "AGB_TLS",
    "dbh_m", "tree_height_m", "csh_raw", "ch_raw", "sbl", "sbd", "sbr",
    "branch_len", "clvr_m.2", "branch_vol_m3", "cd_raw",
    "projected_area_m2", "alpha_volume_m3",
    "tree_vol_m3", "trunk_vol_m3", "base_vol_0_10"
  )

  angles_vars   <- c("sba_degrees", "ba2")
  unitless_vars <- c("bar", "rvr", "cdhr")

  log_transform_cols <- function(dat, cols, eps = 1e-6) {
    out <- dat
    info <- list()

    for (v in cols) {
      if (!v %in% names(out)) next
      if (!is.numeric(out[[v]])) next

      x <- out[[v]]
      x_ok <- x[is.finite(x)]
      if (length(x_ok) == 0) next

      minx <- min(x_ok, na.rm = TRUE)
      shift <- ifelse(minx <= 0, abs(minx) + eps, 0)
      out[[v]] <- log(x + shift)

      info[[v]] <- tibble::tibble(
        variable = v,
        min_before = minx,
        shift_added = shift
      )
    }

    list(data = out, info = dplyr::bind_rows(info))
  }

  log_transform_selective <- function(dat, cols) {
    cols <- unique(cols)
    cols <- cols[cols %in% names(dat)]
    cols <- cols[cols %in% metric_vars_allow]
    cols <- setdiff(cols, c(angles_vars, unitless_vars))
    log_transform_cols(dat, cols)
  }

  center_predictors <- function(dat, vars) {
    out <- dat

    for (v in vars) {
      if (!v %in% names(out)) next
      if (!is.numeric(out[[v]])) next

      ok <- is.finite(out[[v]])
      out[[v]][ok] <- out[[v]][ok] - mean(out[[v]][ok], na.rm = TRUE)
    }

    out
  }

  make_3class_means <- function(x) {
    out_class <- rep(NA_character_, length(x))
    ok <- is.finite(x)

    if (sum(ok) < 3) {
      return(c(low = NA_real_, medium = NA_real_, high = NA_real_))
    }

    grp <- dplyr::ntile(x[ok], 3)
    out_class[ok] <- c("low", "medium", "high")[grp]
    out_class <- factor(out_class, levels = c("low", "medium", "high"))

    means <- tapply(x, out_class, mean, na.rm = TRUE)

    c(
      low = unname(means["low"]),
      medium = unname(means["medium"]),
      high = unname(means["high"])
    )
  }

  make_3class_means_for_group <- function(dat, var, group_level) {
    x <- dat |>
      dplyr::filter(.data[[group_var]] == group_level) |>
      dplyr::pull(.data[[var]])

    make_3class_means(x)
  }

  group_mean <- function(dat, var, group_level) {
    dat |>
      dplyr::filter(.data[[group_var]] == group_level) |>
      dplyr::summarise(value = mean(.data[[var]], na.rm = TRUE)) |>
      dplyr::pull(value)
  }

  fit_lmer <- function(formula_obj, dat) {
    lme4::lmer(
      formula_obj,
      data = dat,
      REML = FALSE,
      na.action = na.omit,
      control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5),
        check.scaleX = "warning"
      )
    )
  }

  quad_term <- function(var) {
    paste0("I(", var, "^2)")
  }

  build_group_formula <- function(response_var,
                                  group_var,
                                  predictor_vars,
                                  quadratic_vars,
                                  interaction_pairs,
                                  random_effect_var) {

    main_terms <- paste0(group_var, ":", predictor_vars)

    quadratic_terms <- character(0)
    if (!is.null(quadratic_vars) && length(quadratic_vars) > 0) {
      quadratic_terms <- paste0(group_var, ":I(", quadratic_vars, "^2)")
    }

    interaction_terms <- character(0)
    if (!is.null(interaction_pairs) && length(interaction_pairs) > 0) {
      interaction_terms <- vapply(interaction_pairs, function(x) {
        paste0(group_var, ":", x[1], ":", x[2])
      }, character(1))
    }

    fixed_terms <- c(
      paste0("0 + ", group_var),
      main_terms,
      quadratic_terms,
      interaction_terms
    )

    formula_text <- paste0(
      response_var,
      " ~ ",
      paste(fixed_terms, collapse = " + "),
      " + (1|",
      random_effect_var,
      ")"
    )

    as.formula(formula_text)
  }

  split_term <- function(term) {
    strsplit(term, ":", fixed = TRUE)[[1]]
  }

  get_coef_by_components <- function(model, components) {
    cf <- lme4::fixef(model)
    cn <- names(cf)

    components <- sort(components)

    for (nm in cn) {
      parts <- sort(split_term(nm))
      if (length(parts) == length(components) && all(parts == components)) {
        return(unname(cf[[nm]]))
      }
    }

    0
  }

  group_component <- function(group_level) {
    paste0(group_var, group_level)
  }

  interaction_partners <- function(focal_var, interaction_pairs) {
    if (is.null(interaction_pairs) || length(interaction_pairs) == 0) return(character(0))

    out <- character(0)

    for (pp in interaction_pairs) {
      if (focal_var %in% pp) {
        out <- c(out, setdiff(pp, focal_var))
      }
    }

    unique(out)
  }

  calc_group_local_effect <- function(model,
                                      dat,
                                      focal_var,
                                      group_level,
                                      focal_value,
                                      values = list()) {

    gcomp <- group_component(group_level)

    beta_main <- get_coef_by_components(
      model,
      components = c(gcomp, focal_var)
    )

    beta_quad <- get_coef_by_components(
      model,
      components = c(gcomp, quad_term(focal_var))
    )

    effect <- beta_main + 2 * beta_quad * focal_value

    partners <- interaction_partners(focal_var, interaction_pairs)

    for (partner in partners) {

      if (partner %in% names(values)) {
        partner_value <- values[[partner]]
      } else {
        partner_value <- group_mean(dat, partner, group_level)
      }

      beta_int <- get_coef_by_components(
        model,
        components = c(gcomp, focal_var, partner)
      )

      effect <- effect + beta_int * partner_value
    }

    effect
  }

  plot_cols_default <- c(
    "#faa916", "#05668d", "#7D8FE6", "#6e1423",
    "#034078", "#db7c26", "#5c8001", "#8B5CF6",
    "#0F766E", "#B45309"
  )

  theme_effect <- function(base_size = 19) {
    ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(
        text = ggplot2::element_text(family = "sans", face = "bold", color = "black"),
        plot.background  = ggplot2::element_rect(fill = "white", color = NA),
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(color = "grey84", linewidth = 0.65),
        plot.title = ggplot2::element_text(face = "bold", size = 27, color = "black"),
        plot.subtitle = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold", size = 18, color = "black"),
        axis.title = ggplot2::element_text(size = 21, face = "bold", color = "black"),
        axis.text  = ggplot2::element_text(size = 17, face = "bold", color = "black"),
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 19, face = "bold", color = "black"),
        legend.text  = ggplot2::element_text(size = 17, face = "bold", color = "black"),
        panel.spacing = grid::unit(1.35, "lines")
      )
  }

  make_diag_plots_lmm <- function(model, dat, out_dir, id_col = "tls_id", label_n = 5) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    fitted <- stats::fitted(model)
    resid <- stats::resid(model)

    stdres <- tryCatch(
      as.numeric(stats::residuals(model)),
      error = function(e) {
        sdr <- stats::sd(resid, na.rm = TRUE)
        if (is.finite(sdr) && sdr > 0) resid / sdr else resid
      }
    )

    mf <- model.frame(model)
    used_rows <- rownames(mf)

    if (!is.null(used_rows) && length(used_rows) == length(fitted)) {
      suppressWarnings(used_idx <- as.integer(used_rows))
      if (all(is.finite(used_idx)) && length(used_idx) == length(fitted)) {
        dd <- dat[used_idx, , drop = FALSE]
      } else {
        dd <- dat[seq_along(fitted), , drop = FALSE]
      }
    } else {
      dd <- dat[seq_along(fitted), , drop = FALSE]
    }

    X <- model.matrix(model)
    y <- model.response(model.frame(model))
    dfX <- as.data.frame(X)
    dfX$.y <- y

    lm_approx <- tryCatch(stats::lm(.y ~ . - 1, data = dfX), error = function(e) NULL)

    leverage <- if (!is.null(lm_approx)) {
      tryCatch(stats::hatvalues(lm_approx), error = function(e) rep(NA_real_, length(fitted)))
    } else {
      rep(NA_real_, length(fitted))
    }

    cooks <- if (!is.null(lm_approx)) {
      tryCatch(stats::cooks.distance(lm_approx), error = function(e) rep(NA_real_, length(fitted)))
    } else {
      rep(NA_real_, length(fitted))
    }

    dd <- tibble::as_tibble(dd) |>
      dplyr::mutate(
        .id = .data[[id_col]],
        .fitted = fitted,
        .resid = resid,
        .stdres = stdres,
        .absstd = abs(stdres),
        .sqrtabsstd = sqrt(abs(stdres)),
        .leverage = leverage,
        .cooks = cooks
      )

    safe_write_csv(dd, file.path(out_dir, "DIAGNOSTIC_VALUES.csv"))

    label_df <- dd |>
      dplyr::arrange(dplyr::desc(.absstd)) |>
      dplyr::slice_head(n = label_n)

    lev_label_df <- dd |>
      dplyr::arrange(dplyr::desc(.cooks)) |>
      dplyr::slice_head(n = label_n)

    p1 <- ggplot2::ggplot(dd, ggplot2::aes(x = .fitted, y = .resid)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
      ggplot2::geom_point(color = "#5E8A6A", alpha = 0.50, size = 2.3) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 1.1, color = "#2F5D46") +
      ggplot2::geom_smooth(method = "loess", se = FALSE, linewidth = 1.1, linetype = "dashed", color = "grey35") +
      ggrepel::geom_text_repel(data = label_df, ggplot2::aes(label = .id), size = 3.4, fontface = "bold") +
      ggplot2::labs(title = "Residuals vs fitted", x = "Fitted value", y = "Residual") +
      theme_effect(base_size = 13) +
      ggplot2::theme(legend.position = "none")

    p2 <- ggplot2::ggplot(dd, ggplot2::aes(sample = .stdres)) +
      ggplot2::stat_qq(color = "#5E8A6A", alpha = 0.60, size = 2.1) +
      ggplot2::stat_qq_line(color = "#2F5D46", linewidth = 1.0) +
      ggplot2::labs(title = "Q-Q plot", x = "Theoretical quantiles", y = "Scaled residual") +
      theme_effect(base_size = 13) +
      ggplot2::theme(legend.position = "none")

    p3 <- ggplot2::ggplot(dd, ggplot2::aes(x = .fitted, y = .sqrtabsstd)) +
      ggplot2::geom_point(color = "#5E8A6A", alpha = 0.50, size = 2.3) +
      ggplot2::geom_smooth(method = "loess", se = FALSE, linewidth = 1.1, linetype = "dashed", color = "#2F5D46") +
      ggrepel::geom_text_repel(data = label_df, ggplot2::aes(label = .id), size = 3.4, fontface = "bold") +
      ggplot2::labs(title = "Scale-location", x = "Fitted value", y = expression(sqrt("|scaled residual|"))) +
      theme_effect(base_size = 13) +
      ggplot2::theme(legend.position = "none")

    p4 <- ggplot2::ggplot(dd, ggplot2::aes(x = .leverage, y = .stdres)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
      ggplot2::geom_point(ggplot2::aes(size = .cooks), color = "#5E8A6A", alpha = 0.60) +
      ggplot2::scale_size_continuous(range = c(2.0, 5.4), guide = "none") +
      ggrepel::geom_text_repel(data = lev_label_df, ggplot2::aes(label = .id), size = 3.4, fontface = "bold") +
      ggplot2::labs(title = "Residuals vs leverage", x = "Leverage", y = "Scaled residual") +
      theme_effect(base_size = 13) +
      ggplot2::theme(legend.position = "none")

    p_combined <- (p1 + p2) / (p3 + p4) +
      patchwork::plot_annotation(title = "Extended diagnostics")

    safe_ggsave(file.path(out_dir, "diag_resid_fit.png"), p1, 8.6, 5.6, 340)
    safe_ggsave(file.path(out_dir, "diag_qq.png"), p2, 8.6, 5.6, 340)
    safe_ggsave(file.path(out_dir, "diag_scale_location.png"), p3, 8.6, 5.6, 340)
    safe_ggsave(file.path(out_dir, "diag_resid_leverage.png"), p4, 8.6, 5.6, 340)
    safe_ggsave(file.path(out_dir, "diag_extended_4panel.png"), p_combined, 15.5, 11.5, 340)

    invisible(dd)
  }

  run_tag <- paste0(
    safe_filename(analysis_tag, 30),
    "_",
    safe_filename(response_var, 25),
    "_BY_",
    safe_filename(group_var, 25),
    "_",
    format(Sys.time(), "%Y%m%d_%H%M%S")
  )

  out_root <- file.path(out_base, run_tag)

  dir_model <- file.path(out_root, "01_model")
  dir_coef  <- file.path(out_root, "02_coefficients")
  dir_eff   <- file.path(out_root, "03_group_local_effects")
  dir_diag  <- file.path(out_root, "04_diagnostics")
  dir_lines <- file.path(out_root, "05_group_lineplots")

  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_model, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_coef, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_eff, recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(run_diagnostics)) {
    dir.create(dir_diag, recursive = TRUE, showWarnings = FALSE)
  }

  if (isTRUE(make_group_lineplots)) {
    dir.create(dir_lines, recursive = TRUE, showWarnings = FALSE)
  }

  raw <- read_messy_csv(data_path)

  needed_cols <- unique(c(
    id_col,
    random_effect_var,
    group_var,
    response_var,
    predictor_vars
  ))

  missing_cols <- setdiff(needed_cols, names(raw))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }

  dat0 <- raw |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c(response_var, predictor_vars)),
        ~ suppressWarnings(as.numeric(.x))
      ),
      !!id_col := as.character(.data[[id_col]]),
      !!random_effect_var := as.factor(.data[[random_effect_var]]),
      !!group_var := make_group_factor(.data[[group_var]], ref = group_reference)
    ) |>
    dplyr::filter(
      is.finite(.data[[response_var]]),
      dplyr::if_all(dplyr::all_of(predictor_vars), is.finite),
      !is.na(.data[[random_effect_var]]),
      !is.na(.data[[group_var]])
    )

  safe_write_csv(dat0, file.path(dir_model, "DATA_USED_BEFORE_TRANSFORM.csv"))

  log_vars <- c(response_var, predictor_vars)
  log_out <- log_transform_selective(dat0, log_vars)
  dat <- log_out$data

  if (!is.null(log_out$info) && nrow(log_out$info) > 0) {
    safe_write_csv(log_out$info, file.path(dir_model, "SHIFT_LOGSELECTIVE.csv"))
  }

  dat <- dat |>
    dplyr::select(dplyr::all_of(c(id_col, random_effect_var, group_var, response_var, predictor_vars))) |>
    tidyr::drop_na()

  dat <- center_predictors(dat, predictor_vars)

  safe_write_csv(dat, file.path(dir_model, "DATA_USED_AFTER_LOG_AND_CENTERING.csv"))

  if (is.null(quadratic_vars)) {
    quadratic_vars <- character(0)
  }

  if (is.null(interaction_pairs)) {
    interaction_pairs <- list()
  }

  model_formula <- build_group_formula(
    response_var = response_var,
    group_var = group_var,
    predictor_vars = predictor_vars,
    quadratic_vars = quadratic_vars,
    interaction_pairs = interaction_pairs,
    random_effect_var = random_effect_var
  )

  model <- fit_lmer(model_formula, dat)

  safe_write_lines(capture.output(summary(model)), file.path(dir_model, "MODEL_SUMMARY.txt"))
  safe_write_lines(deparse(model_formula), file.path(dir_model, "MODEL_FORMULA.txt"))

  fixed_tbl <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  random_tbl <- broom.mixed::tidy(model, effects = "ran_pars", conf.int = TRUE)

  safe_write_csv(fixed_tbl, file.path(dir_coef, "TIDY_FIXED.csv"))
  safe_write_csv(random_tbl, file.path(dir_model, "TIDY_RANDOM.csv"))

  fit_stats <- tibble::tibble(
    AIC = stats::AIC(model),
    BIC = stats::BIC(model),
    logLik = as.numeric(stats::logLik(model)),
    N = stats::nobs(model),
    groups = paste(levels(dat[[group_var]]), collapse = ", ")
  )

  safe_write_csv(fit_stats, file.path(dir_model, "FIT_STATS.csv"))

  r2_tbl <- tryCatch(
    as.data.frame(performance::r2_nakagawa(model)),
    error = function(e) data.frame(error = e$message)
  )
  safe_write_csv(r2_tbl, file.path(dir_model, "R2.csv"))

  coef_tbl <- fixed_tbl |>
    dplyr::filter(term != "(Intercept)") |>
    dplyr::mutate(
      term_type = dplyr::case_when(
        grepl(":I\\(", term) | grepl("I\\(", term) ~ "Quadratic",
        stringr::str_count(term, ":") >= 2 ~ "Interaction",
        grepl(":", term, fixed = TRUE) ~ "Main effect",
        TRUE ~ "Group intercept"
      ),
      term_label = term
    ) |>
    dplyr::arrange(estimate) |>
    dplyr::mutate(term_label = factor(term_label, levels = term_label))

  if (isTRUE(make_coef_plot)) {
    p_coef <- ggplot2::ggplot(coef_tbl, ggplot2::aes(x = estimate, y = term_label, color = term_type)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 1.1) +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = conf.low, xmax = conf.high),
        height = 0,
        linewidth = 1.5,
        lineend = "round"
      ) +
      ggplot2::geom_point(size = 4.2) +
      ggplot2::scale_color_manual(
        values = c(
          "Group intercept" = "grey35",
          "Main effect" = "#006D77",
          "Quadratic" = "#BC6C25",
          "Interaction" = "#B22234"
        )
      ) +
      ggplot2::labs(
        title = "Group-specific fixed-effect coefficients",
        x = "Estimate",
        y = NULL,
        color = NULL
      ) +
      theme_effect(base_size = 13) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 9.5, face = "bold", color = "black"),
        legend.position = "bottom"
      )

    safe_ggsave(file.path(dir_coef, "coef_fixed.png"), p_coef, 15.5, 13.5, 420)
  }

  safe_write_csv(coef_tbl, file.path(dir_coef, "coef_fixed_table.csv"))

  if (isTRUE(run_diagnostics)) {
    make_diag_plots_lmm(
      model = model,
      dat = dat,
      out_dir = dir_diag,
      id_col = id_col,
      label_n = 5
    )
  }

  group_levels <- levels(dat[[group_var]])
  all_effects <- list()

  if (length(interaction_pairs) > 0) {

    for (pair in interaction_pairs) {

      focal_var <- pair[1]
      partner_var <- pair[2]

      grid_df <- expand.grid(
        focal_class = c("low", "medium", "high"),
        partner_class = c("low", "medium", "high"),
        group_level = group_levels,
        stringsAsFactors = FALSE
      ) |>
        tibble::as_tibble()

      grid_df <- grid_df |>
        dplyr::rowwise() |>
        dplyr::mutate(
          focal_value = unname(make_3class_means_for_group(dat, focal_var, group_level)[focal_class]),
          partner_value = unname(make_3class_means_for_group(dat, partner_var, group_level)[partner_class])
        ) |>
        dplyr::ungroup()

      grid_df$local_beta <- purrr::pmap_dbl(
        list(grid_df$focal_value, grid_df$partner_value, grid_df$group_level),
        function(focal_value, partner_value, group_level) {
          values_tmp <- list()
          values_tmp[[partner_var]] <- partner_value

          calc_group_local_effect(
            model = model,
            dat = dat,
            focal_var = focal_var,
            group_level = group_level,
            focal_value = focal_value,
            values = values_tmp
          )
        }
      )

      grid_df <- grid_df |>
        dplyr::mutate(
          effect_type = "interaction",
          focal_var = focal_var,
          partner_var = partner_var,
          focal_pretty = label_term(focal_var),
          partner_pretty = label_term(partner_var),
          focal_clean = strip_units_label(focal_pretty),
          partner_clean = strip_units_label(partner_pretty),
          approx_pct = 100 * local_beta,
          exact_pct = 100 * (exp(local_beta) - 1),
          focal_class = factor(focal_class, levels = c("low", "medium", "high")),
          partner_class = factor(partner_class, levels = c("low", "medium", "high")),
          group_level = factor(group_level, levels = group_levels),
          group_label = paste0(group_var, " ", group_level),
          facet_label = factor(
            paste0(partner_clean, ": ", as.character(partner_class)),
            levels = paste0(unique(partner_clean)[1], ": ", c("low", "medium", "high"))
          )
        )

      out_name <- paste0(
        "EFFECT_",
        safe_filename(focal_var),
        "_BY_",
        safe_filename(partner_var)
      )

      safe_write_csv(grid_df, file.path(dir_eff, paste0(out_name, ".csv")))
      all_effects[[out_name]] <- grid_df

      pal <- plot_cols_default[seq_along(group_levels)]
      names(pal) <- group_levels

      p <- ggplot2::ggplot(
        grid_df,
        ggplot2::aes(
          x = as.numeric(focal_class),
          y = local_beta,
          color = group_level,
          group = group_level
        )
      ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 1.0) +
        ggplot2::geom_line(linewidth = 2.0, alpha = 0.85, lineend = "butt", linejoin = "mitre") +
        ggplot2::facet_wrap(~ facet_label, nrow = 1) +
        ggplot2::scale_x_continuous(
          breaks = c(1, 2, 3),
          labels = c("Low", "Medium", "High"),
          expand = ggplot2::expansion(mult = c(0.08, 0.08))
        ) +
        ggplot2::scale_color_manual(
          values = pal,
          name = group_var
        ) +
        ggplot2::labs(
          title = paste0(unique(grid_df$focal_clean)[1], " effect to ", strip_units_label(label_term(response_var))),
          x = paste0(unique(grid_df$focal_clean)[1], " class"),
          y = paste0("Local marginal effect to ", strip_units_label(label_term(response_var)))
        ) +
        theme_effect(base_size = 19)

      if (isTRUE(make_group_lineplots)) {
        safe_ggsave(file.path(dir_lines, paste0(out_name, ".png")), p, 16.5, 7.2, 600)
      }
    }
  }

  for (focal_var in predictor_vars) {

    grid_df <- expand.grid(
      focal_class = c("low", "medium", "high"),
      group_level = group_levels,
      stringsAsFactors = FALSE
    ) |>
      tibble::as_tibble()

    grid_df <- grid_df |>
      dplyr::rowwise() |>
      dplyr::mutate(
        focal_value = unname(make_3class_means_for_group(dat, focal_var, group_level)[focal_class])
      ) |>
      dplyr::ungroup()

    grid_df$local_beta <- purrr::pmap_dbl(
      list(grid_df$focal_value, grid_df$group_level),
      function(focal_value, group_level) {
        calc_group_local_effect(
          model = model,
          dat = dat,
          focal_var = focal_var,
          group_level = group_level,
          focal_value = focal_value,
          values = list()
        )
      }
    )

    grid_df <- grid_df |>
      dplyr::mutate(
        effect_type = "single",
        focal_var = focal_var,
        partner_var = NA_character_,
        focal_pretty = label_term(focal_var),
        focal_clean = strip_units_label(focal_pretty),
        approx_pct = 100 * local_beta,
        exact_pct = 100 * (exp(local_beta) - 1),
        focal_class = factor(focal_class, levels = c("low", "medium", "high")),
        group_level = factor(group_level, levels = group_levels),
        group_label = paste0(group_var, " ", group_level)
      )

    out_name <- paste0("SINGLE_EFFECT_", safe_filename(focal_var))

    safe_write_csv(grid_df, file.path(dir_eff, paste0(out_name, ".csv")))
    all_effects[[out_name]] <- grid_df

    pal <- plot_cols_default[seq_along(group_levels)]
    names(pal) <- group_levels

    p <- ggplot2::ggplot(
      grid_df,
      ggplot2::aes(
        x = as.numeric(focal_class),
        y = local_beta,
        color = group_level,
        group = group_level
      )
    ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 1.0) +
      ggplot2::geom_line(linewidth = 2.0, alpha = 0.85, lineend = "butt", linejoin = "mitre") +
      ggplot2::scale_x_continuous(
        breaks = c(1, 2, 3),
        labels = c("Low", "Medium", "High"),
        expand = ggplot2::expansion(mult = c(0.08, 0.08))
      ) +
      ggplot2::scale_color_manual(
        values = pal,
        name = group_var
      ) +
      ggplot2::labs(
        title = paste0(unique(grid_df$focal_clean)[1], " effect to ", strip_units_label(label_term(response_var))),
        x = paste0(unique(grid_df$focal_clean)[1], " class"),
        y = paste0("Local marginal effect to ", strip_units_label(label_term(response_var)))
      ) +
      theme_effect(base_size = 19)

    if (isTRUE(make_group_lineplots)) {
      safe_ggsave(file.path(dir_lines, paste0(out_name, ".png")), p, 11.5, 7.2, 600)
    }
  }

  all_effects_tbl <- dplyr::bind_rows(all_effects, .id = "effect_id")
  safe_write_csv(all_effects_tbl, file.path(dir_eff, "ALL_GROUP_LOCAL_EFFECTS.csv"))

  safe_write_lines(
    c(
      "TreeArchi group-specific MLR",
      paste0("Created: ", Sys.time()),
      paste0("Data path: ", data_path),
      paste0("Output root: ", out_root),
      paste0("Response: ", response_var),
      paste0("Group variable: ", group_var),
      paste0("Group reference: ", group_reference %||% "none"),
      paste0("Predictors: ", paste(predictor_vars, collapse = ", ")),
      paste0("Quadratic vars: ", paste(quadratic_vars, collapse = ", ")),
      paste0(
        "Interaction pairs: ",
        paste(vapply(interaction_pairs, paste, collapse = ":", FUN.VALUE = character(1)), collapse = ", ")
      ),
      "",
      "Model formula:",
      paste(deparse(model_formula), collapse = " "),
      "",
      "Metric note:",
      "alpha_volume_m3 = 3D alpha-volume proxy.",
      "projected_area_m2 = 2D projected crown area.",
      "",
      "Interpretation:",
      "- Response and metric predictors are log-transformed selectively.",
      "- Predictors are centered after log transformation.",
      "- Low / medium / high values are group-specific class means from transformed and centered data.",
      "- local_beta is on the log-response scale.",
      "- approx_pct = 100 * local_beta.",
      "- exact_pct = 100 * (exp(local_beta) - 1).",
      "- This function fits one shared model with group-specific terms, not separate models per group."
    ),
    file.path(out_root, "RUN_INFO.txt")
  )

  message("\nDONE: TreeArchi group MLR")
  message("Output root: ", out_root)

  invisible(list(
    model = model,
    formula = model_formula,
    data = dat,
    fixed_effects = fixed_tbl,
    local_effects = all_effects_tbl,
    out_root = out_root
  ))
}
