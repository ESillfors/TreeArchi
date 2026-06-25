coef_plot_lmm <- function(model, model_name, out_dir, labels_map, dat = NULL, response = NULL,
                          local_effect_summary = NULL) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  td0 <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(
      term_pretty = vapply(term, pretty_term, character(1), labels_map = labels_map),
      term_class = dplyr::case_when(
        grepl(":", term, fixed = TRUE) ~ "Interaction",
        grepl("^I\\(.+\\^2\\)$", term) ~ "Quadratic",
        TRUE ~ "Main"
      )
    )

  if (!is.null(dat) && !is.null(response)) {
    ft <- flip_test(model, dat, response = response, x_vars = names(dat))
  } else {
    ft <- tibble::tibble(
      predictor = character(),
      raw_slope = numeric(),
      adj_slope = numeric(),
      flip = logical()
    )
  }

  td_plot <- td0 %>%
    dplyr::left_join(ft, by = c("term" = "predictor"))

  if (!is.null(local_effect_summary) && nrow(local_effect_summary) > 0) {
    td_plot <- td_plot %>%
      dplyr::left_join(local_effect_summary, by = c("term" = "term"))
  } else {
    td_plot <- td_plot %>%
      dplyr::mutate(
        overall_local_beta = NA_real_,
        low_local_beta = NA_real_,
        mid_local_beta = NA_real_,
        high_local_beta = NA_real_
      )
  }

  if (!"p.value" %in% names(td_plot)) {
    td_plot$p.value <- NA_real_
  }

  td_plot <- td_plot %>%
    dplyr::mutate(
      p_label = ifelse(is.na(p.value), "", paste0("p ", safe_p(p.value))),
      estimate_for_order = dplyr::coalesce(estimate, overall_local_beta, raw_slope),
      x_adjusted = estimate,
      x_raw = raw_slope - 0.10,
      x_overall = overall_local_beta + 0.10,
      x_low = low_local_beta - 0.22,
      x_mid = mid_local_beta,
      x_high = high_local_beta + 0.22
    )

  if (nrow(td_plot) == 0) return(invisible(NULL))

  td_plot <- td_plot %>%
    dplyr::arrange(dplyr::if_else(is.na(estimate_for_order), Inf, estimate_for_order)) %>%
    dplyr::mutate(term_pretty = factor(term_pretty, levels = unique(term_pretty)))

  x_all <- c(
    td_plot$conf.low, td_plot$conf.high, td_plot$x_adjusted, td_plot$x_raw,
    td_plot$x_overall, td_plot$x_low, td_plot$x_mid, td_plot$x_high
  )
  x_all <- x_all[is.finite(x_all)]

  if (length(x_all) == 0) {
    xr <- c(-1, 1)
  } else {
    xr <- range(x_all, na.rm = TRUE)
    if (!is.finite(diff(xr)) || diff(xr) == 0) xr <- xr + c(-0.5, 0.5)
  }

  xpad <- max(0.06 * diff(xr), 0.08)

  td_plot <- td_plot %>%
    dplyr::mutate(
      p_x = pmax(
        ifelse(is.finite(conf.high), conf.high, NA_real_),
        ifelse(is.finite(x_adjusted), x_adjusted, NA_real_),
        ifelse(is.finite(x_raw), x_raw, NA_real_),
        ifelse(is.finite(x_overall), x_overall, NA_real_),
        ifelse(is.finite(x_low), x_low, NA_real_),
        ifelse(is.finite(x_mid), x_mid, NA_real_),
        ifelse(is.finite(x_high), x_high, NA_real_),
        na.rm = TRUE
      ) + xpad
    )

  df_ci <- td_plot %>% dplyr::filter(is.finite(conf.low), is.finite(conf.high), is.finite(x_adjusted))
  df_est <- td_plot %>% dplyr::filter(is.finite(x_adjusted))
  df_raw <- td_plot %>% dplyr::filter(is.finite(x_raw))
  df_overall <- td_plot %>% dplyr::filter(is.finite(x_overall))
  df_low <- td_plot %>% dplyr::filter(is.finite(x_low))
  df_mid <- td_plot %>% dplyr::filter(is.finite(x_mid))
  df_high <- td_plot %>% dplyr::filter(is.finite(x_high))
  df_flip <- td_plot %>% dplyr::filter(!is.na(flip), flip == TRUE, is.finite(x_adjusted))

  p <- ggplot2::ggplot(td_plot, ggplot2::aes(y = term_pretty)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = COL$zero, linewidth = 0.8) +

    ggplot2::geom_segment(
      data = df_ci,
      ggplot2::aes(x = conf.low, xend = conf.high, y = term_pretty, yend = term_pretty),
      inherit.aes = FALSE,
      linewidth = 0.8, color = COL$muted, alpha = 0.95
    ) +

    ggplot2::geom_point(
      data = df_est,
      ggplot2::aes(x = x_adjusted, y = term_pretty, shape = "Selected estimate", fill = "Selected estimate"),
      inherit.aes = FALSE,
      size = 3.0, stroke = 0.8, color = COL$ink
    ) +

    ggplot2::geom_point(
      data = df_raw,
      ggplot2::aes(x = x_raw, y = term_pretty, shape = "Raw β", fill = "Raw β"),
      inherit.aes = FALSE,
      size = 2.8, stroke = 0.8, color = COL$ink
    ) +

    ggplot2::geom_point(
      data = df_overall,
      ggplot2::aes(x = x_overall, y = term_pretty, shape = "Overall", fill = "Overall"),
      inherit.aes = FALSE,
      size = 2.8, stroke = 0.8, color = COL$ink
    ) +

    ggplot2::geom_point(
      data = df_low,
      ggplot2::aes(x = x_low, y = term_pretty, shape = "Low", fill = "Low"),
      inherit.aes = FALSE,
      size = 2.7, stroke = 0.8, color = COL$ink
    ) +

    ggplot2::geom_point(
      data = df_mid,
      ggplot2::aes(x = x_mid, y = term_pretty, shape = "Mid", fill = "Mid"),
      inherit.aes = FALSE,
      size = 2.7, stroke = 0.8, color = COL$ink
    ) +

    ggplot2::geom_point(
      data = df_high,
      ggplot2::aes(x = x_high, y = term_pretty, shape = "High", fill = "High"),
      inherit.aes = FALSE,
      size = 2.7, stroke = 0.8, color = COL$ink
    ) +

    ggplot2::geom_point(
      data = df_flip,
      ggplot2::aes(x = x_adjusted, y = term_pretty),
      inherit.aes = FALSE,
      shape = 8, size = 4, color = "red"
    ) +

    ggplot2::geom_text(
      data = td_plot %>% dplyr::filter(p_label != ""),
      ggplot2::aes(x = p_x, y = term_pretty, label = p_label),
      inherit.aes = FALSE,
      hjust = 0, size = 3.2, color = COL$ink
    ) +

    ggplot2::scale_shape_manual(
      values = c(
        "Selected estimate" = 21,
        "Raw β" = 23,
        "Overall" = 22,
        "Low" = 24,
        "Mid" = 21,
        "High" = 25
      ),
      breaks = c("Selected estimate", "Raw β", "Overall", "Low", "Mid", "High")
    ) +

    ggplot2::scale_fill_manual(
      values = c(
        "Selected estimate" = COL$point,
        "Raw β" = "white",
        "Overall" = COL$overall,
        "Low" = COL$low,
        "Mid" = COL$mid,
        "High" = COL$high
      ),
      breaks = c("Selected estimate", "Raw β", "Overall", "Low", "Mid", "High")
    ) +

    ggplot2::labs(
      title = "Fixed-effect coefficients",
      subtitle = paste0(model_name, " • interaction terms included as rows"),
      x = "Adjusted estimate / marginal effect (95% CI where available)",
      y = NULL,
      shape = NULL,
      fill = NULL
    ) +

    ggplot2::coord_cartesian(
      xlim = c(min(xr) - 0.15 * diff(xr), max(td_plot$p_x, na.rm = TRUE) + 0.05 * diff(xr)),
      clip = "off"
    ) +

    theme_modern() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = COL$tint_dark, color = NA),
      plot.margin = ggplot2::margin(8, 52, 8, 8),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.text = ggplot2::element_text(size = 10)
    ) +

    ggplot2::guides(
      shape = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 3.5)),
      fill  = ggplot2::guide_legend(
        nrow = 1,
        override.aes = list(shape = c(21, 23, 22, 24, 21, 25), size = 3.5)
      )
    )

  safe_ggsave(
    file.path(out_dir, "coef_fixed.png"),
    p,
    width = 11.8,
    height = max(5.2, 0.46 * nrow(td_plot) + 2.6),
    dpi = 340
  )

  safe_write_csv(td_plot, file.path(out_dir, "coef_fixed_table.csv"))

  p
}
