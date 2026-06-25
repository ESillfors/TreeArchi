# -------------------------
# Diagnostics helpers
# -------------------------

get_fixedpart_leverage <- function(model) {
  y <- model.response(model.frame(model))
  X <- model.matrix(model)

  dfX <- as.data.frame(X)
  dfX$.y <- y

  lm_approx <- tryCatch(
    stats::lm(.y ~ . - 1, data = dfX),
    error = function(e) NULL
  )

  if (is.null(lm_approx)) {
    n <- NROW(y)
    return(list(
      leverage = rep(NA_real_, n),
      cooks = rep(NA_real_, n)
    ))
  }

  list(
    leverage = tryCatch(
      stats::hatvalues(lm_approx),
      error = function(e) rep(NA_real_, NROW(y))
    ),
    cooks = tryCatch(
      stats::cooks.distance(lm_approx),
      error = function(e) rep(NA_real_, NROW(y))
    )
  )
}


make_qq_envelope_normal <- function(stdres, conf = 0.99, B = 2000, seed = 123) {

  stdres <- stdres[is.finite(stdres)]
  n <- length(stdres)

  if (n < 10) return(NULL)

  set.seed(seed)

  theo <- stats::qnorm(stats::ppoints(n))
  samp <- sort(stdres)

  q_theo <- stats::qnorm(c(0.25, 0.75))
  q_samp <- stats::quantile(
    samp,
    probs = c(0.25, 0.75),
    na.rm = TRUE,
    names = FALSE
  )

  slope <- unname((q_samp[2] - q_samp[1]) / (q_theo[2] - q_theo[1]))
  inter <- unname(q_samp[1] - slope * q_theo[1])

  sim_mat <- replicate(B, sort(stats::rnorm(n)))
  alpha <- (1 - conf) / 2

  lower_std <- apply(sim_mat, 1, stats::quantile, probs = alpha, na.rm = TRUE)
  upper_std <- apply(sim_mat, 1, stats::quantile, probs = 1 - alpha, na.rm = TRUE)
  med_std   <- apply(sim_mat, 1, stats::median, na.rm = TRUE)

  tibble::tibble(
    theoretical = theo,
    sample = samp,
    lower = inter + slope * lower_std,
    upper = inter + slope * upper_std,
    median_sim = inter + slope * med_std,
    qq_intercept = inter,
    qq_slope = slope
  )
}


# -------------------------
# Diagnostics
# -------------------------

make_diag_plots_lmm <- function(model,
                                dat,
                                response_label,
                                out_dir,
                                id_col = "tls_id",
                                label_n = 5) {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  fitted <- stats::fitted(model)
  resid  <- stats::resid(model)
  mf <- stats::model.frame(model)
  row_ids <- as.integer(rownames(mf))

  if (all(is.finite(row_ids)) && max(row_ids, na.rm = TRUE) <= nrow(dat)) {
    dat <- dat[row_ids, , drop = FALSE]
  } else {
    dat <- dat[seq_len(length(fitted)), , drop = FALSE]
  }

  stdres <- tryCatch(
    as.numeric(stats::residuals(model)),
    error = function(e) {
      sdr <- stats::sd(resid, na.rm = TRUE)
      if (is.finite(sdr) && sdr > 0) resid / sdr else resid
    }
  )

  lev_obj <- get_fixedpart_leverage(model)

  if (!id_col %in% names(dat)) {
    dat[[id_col]] <- seq_len(nrow(dat))
  }

  dd <- dat |>
    dplyr::mutate(
      .id = as.character(.data[[id_col]]),
      .fitted = fitted,
      .resid = resid,
      .stdres = stdres,
      .absstd = abs(stdres),
      .sqrt_abs_stdres = sqrt(abs(stdres)),
      .leverage = lev_obj$leverage,
      .cooks = lev_obj$cooks
    )

  safe_write_csv(
    dd,
    file.path(out_dir, "diagnostic_values.csv")
  )

  label_df <- dd |>
    dplyr::filter(is.finite(.absstd)) |>
    dplyr::arrange(dplyr::desc(.absstd)) |>
    dplyr::slice_head(n = label_n)

  # 1. Residuals vs fitted
  p1 <- ggplot2::ggplot(dd, ggplot2::aes(x = .fitted, y = .resid)) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = COL$zero,
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      color = COL$point,
      alpha = 0.55,
      size = 2
    ) +
    ggplot2::geom_smooth(
      method = "loess",
      se = FALSE,
      color = COL$loess,
      linewidth = 0.8
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = .id),
      color = COL$ink,
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0
    ) +
    ggplot2::labs(
      title = "Residuals vs fitted",
      subtitle = response_label,
      x = "Fitted value",
      y = "Residual"
    ) +
    theme_modern() +
    ggplot2::theme(legend.position = "none")

  safe_ggsave(
    file.path(out_dir, "diag_resid_fit.png"),
    p1,
    width = 8.6,
    height = 5.6,
    dpi = 340
  )

  # 2. Q-Q plot with simulated 99% envelope
  qq_env <- make_qq_envelope_normal(
    dd$.stdres,
    conf = 0.99,
    B = 2000,
    seed = 123
  )

  if (is.null(qq_env)) {
    qq_tmp <- dd |>
      dplyr::filter(is.finite(.stdres)) |>
      dplyr::mutate(.theoretical = stats::qqnorm(.stdres, plot.it = FALSE)$x)

    q_theo <- stats::qnorm(c(0.25, 0.75))
    q_samp <- stats::quantile(
      sort(qq_tmp$.stdres),
      probs = c(0.25, 0.75),
      na.rm = TRUE,
      names = FALSE
    )

    slope <- unname((q_samp[2] - q_samp[1]) / (q_theo[2] - q_theo[1]))
    inter <- unname(q_samp[1] - slope * q_theo[1])

    qq_env <- tibble::tibble(
      theoretical = sort(qq_tmp$.theoretical),
      sample = sort(qq_tmp$.stdres),
      lower = NA_real_,
      upper = NA_real_,
      median_sim = NA_real_,
      qq_intercept = inter,
      qq_slope = slope
    )
  }

  qq_ids <- dd |>
    dplyr::filter(is.finite(.stdres)) |>
    dplyr::arrange(.stdres) |>
    dplyr::mutate(.theoretical = qq_env$theoretical) |>
    dplyr::arrange(dplyr::desc(abs(.stdres))) |>
    dplyr::slice_head(n = label_n)

  p2 <- ggplot2::ggplot(qq_env, ggplot2::aes(x = theoretical, y = sample)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = COL$qq_band,
      alpha = 0.28
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = median_sim),
      color = COL$qq_line,
      linewidth = 0.50,
      alpha = 0.50
    ) +
    ggplot2::geom_abline(
      intercept = qq_env$qq_intercept[1],
      slope = qq_env$qq_slope[1],
      color = COL$qq_line,
      linewidth = 0.55,
      alpha = 0.75
    ) +
    ggplot2::geom_point(
      color = COL$point,
      alpha = 0.48,
      size = 1.8
    ) +
    ggrepel::geom_text_repel(
      data = qq_ids,
      ggplot2::aes(x = .theoretical, y = .stdres, label = .id),
      inherit.aes = FALSE,
      color = COL$ink,
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0
    ) +
    ggplot2::labs(
      title = "Q–Q plot (scaled residuals)",
      subtitle = paste0(response_label, " • fitted qqline with simulated 99% envelope"),
      x = "Theoretical quantiles",
      y = "Scaled residual"
    ) +
    theme_modern() +
    ggplot2::theme(legend.position = "none")

  safe_ggsave(
    file.path(out_dir, "diag_qq_boot99.png"),
    p2,
    width = 8.6,
    height = 5.6,
    dpi = 340
  )

  # 3. Scale-location
  p3 <- ggplot2::ggplot(dd, ggplot2::aes(x = .fitted, y = .sqrt_abs_stdres)) +
    ggplot2::geom_point(
      color = COL$point,
      alpha = 0.55,
      size = 2
    ) +
    ggplot2::geom_smooth(
      method = "loess",
      se = FALSE,
      color = COL$loess,
      linewidth = 0.8
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = .id),
      color = COL$ink,
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0
    ) +
    ggplot2::labs(
      title = "Scale-location",
      subtitle = response_label,
      x = "Fitted value",
      y = expression(sqrt("|scaled residual|"))
    ) +
    theme_modern() +
    ggplot2::theme(legend.position = "none")

  safe_ggsave(
    file.path(out_dir, "diag_scale_location.png"),
    p3,
    width = 8.6,
    height = 5.6,
    dpi = 340
  )

  # 4. Residuals vs leverage
  lev_label_df <- dd |>
    dplyr::filter(is.finite(.leverage), is.finite(.absstd)) |>
    dplyr::arrange(dplyr::desc(.cooks), dplyr::desc(.absstd)) |>
    dplyr::slice_head(n = label_n)

  p4 <- ggplot2::ggplot(dd, ggplot2::aes(x = .leverage, y = .stdres)) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = COL$zero,
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      ggplot2::aes(size = .cooks),
      color = COL$point,
      alpha = 0.55
    ) +
    ggplot2::scale_size_continuous(
      range = c(1.8, 5.0),
      guide = "none"
    ) +
    ggrepel::geom_text_repel(
      data = lev_label_df,
      ggplot2::aes(label = .id),
      color = COL$ink,
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0
    ) +
    ggplot2::labs(
      title = "Residuals vs leverage",
      subtitle = paste0(response_label, " • leverage/Cook's D from fixed-part approximation"),
      x = "Leverage",
      y = "Scaled residual"
    ) +
    theme_modern() +
    ggplot2::theme(legend.position = "none")

  safe_ggsave(
    file.path(out_dir, "diag_resid_leverage.png"),
    p4,
    width = 8.6,
    height = 5.6,
    dpi = 340
  )

  # Combined 4-panel diagnostic figure
  p_combined <- patchwork::wrap_plots(
    p1, p2, p3, p4,
    ncol = 2
  ) +
    patchwork::plot_annotation(
      title = "Extended diagnostics",
      subtitle = response_label,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14, color = COL$ink),
        plot.subtitle = ggplot2::element_text(size = 10, color = COL$muted),
        plot.background = ggplot2::element_rect(fill = COL$bg, color = NA)
      )
    )

  safe_ggsave(
    file.path(out_dir, "diag_extended_4panel.png"),
    p_combined,
    width = 15.5,
    height = 11.5,
    dpi = 340
  )

  safe_ggsave(
    file.path(out_dir, "diag_extended_4panel.pdf"),
    p_combined,
    width = 15.5,
    height = 11.5,
    dpi = 340
  )

  invisible(list(
    p_resid = p1,
    p_qq = p2,
    p_scale_location = p3,
    p_resid_leverage = p4,
    p_combined = p_combined
  ))
}
