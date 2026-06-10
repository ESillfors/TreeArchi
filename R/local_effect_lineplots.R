theme_effect <- function(base_size = 19) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans", face = "bold", color = "black"),
      plot.background  = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey84", linewidth = 0.65),
      plot.title = ggplot2::element_text(face = "bold", size = 27, color = "black"),
      strip.text = ggplot2::element_text(face = "bold", size = 18, color = "black"),
      axis.title = ggplot2::element_text(size = 21, face = "bold", color = "black"),
      axis.text  = ggplot2::element_text(size = 17, face = "bold", color = "black"),
      legend.position = "none"
    )
}

save_local_effect_lineplots <- function(local_effects,
                                        out_dir,
                                        response_label = NULL) {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(local_effects) || nrow(local_effects) == 0) {
    return(invisible(NULL))
  }

  response_label <- response_label %||% "response"

  plot_df <- local_effects |>
    dplyr::filter(!is.na(partner_var)) |>
    dplyr::mutate(
      focal_class = factor(focal_class, levels = c("low", "medium", "high")),
      partner_class = factor(partner_class, levels = c("low", "medium", "high")),
      focal_clean = short_pretty(pretty_term(focal_var, pretty_names)),
      partner_clean = short_pretty(pretty_term(partner_var, pretty_names)),
      facet_label = factor(
        paste0(partner_clean, ": ", partner_class),
        levels = paste0(unique(partner_clean)[1], ": ", c("low", "medium", "high"))
      )
    )

  if (nrow(plot_df) == 0) {
    message("No interaction local effects found. Skipping local effect lineplots.")
    return(invisible(NULL))
  }

  combos <- plot_df |>
    dplyr::distinct(focal_var, partner_var, focal_clean, partner_clean)

  plots <- list()

  for (i in seq_len(nrow(combos))) {

    focal_i <- combos$focal_var[i]
    partner_i <- combos$partner_var[i]

    df_i <- plot_df |>
      dplyr::filter(
        focal_var == focal_i,
        partner_var == partner_i
      ) |>
      dplyr::mutate(
        facet_label = factor(
          paste0(partner_clean, ": ", partner_class),
          levels = paste0(unique(partner_clean)[1], ": ", c("low", "medium", "high"))
        )
      )

    focal_lab <- combos$focal_clean[i]

    p <- ggplot2::ggplot(
      df_i,
      ggplot2::aes(
        x = as.numeric(focal_class),
        y = local_beta,
        group = 1
      )
    ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "grey55",
        linewidth = 1.0
      ) +
      ggplot2::geom_line(
        linewidth = 2.0,
        color = "#2F5D46",
        alpha = 0.90,
        lineend = "butt",
        linejoin = "mitre"
      ) +
      ggplot2::geom_point(
        size = 4.2,
        color = "black",
        fill = "#5E8A6A",
        shape = 21,
        stroke = 1.0
      ) +
      ggplot2::facet_wrap(~ facet_label, nrow = 1) +
      ggplot2::scale_x_continuous(
        breaks = c(1, 2, 3),
        labels = c("Low", "Medium", "High"),
        expand = ggplot2::expansion(mult = c(0.08, 0.08))
      ) +
      ggplot2::labs(
        title = paste0(focal_lab, " effect to ", short_pretty(response_label)),
        x = paste0(focal_lab, " class"),
        y = paste0("Local marginal effect to ", short_pretty(response_label))
      ) +
      theme_effect(base_size = 19)

    out_name <- paste0(
      "EFFECT_",
      safe_filename(focal_i, 40),
      "_BY_",
      safe_filename(partner_i, 40)
    )

    safe_write_csv(
      df_i,
      file.path(out_dir, paste0(out_name, ".csv"))
    )

    safe_ggsave(file.path(out_dir, paste0(out_name, ".png")), p, width = 16.5, height = 7.2, dpi = 600)
    safe_ggsave(file.path(out_dir, paste0(out_name, ".pdf")), p, width = 16.5, height = 7.2, dpi = 600)
    safe_ggsave(file.path(out_dir, paste0(out_name, ".svg")), p, width = 16.5, height = 7.2, dpi = 600)

    plots[[out_name]] <- p
  }

  invisible(plots)
}
