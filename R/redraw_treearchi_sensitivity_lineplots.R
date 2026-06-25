#' Redraw TreeArchi sensitivity lineplots
#'
#' @param sensitivity_dir Sensitivity output directory.
#' @param effects_file Optional collected effects CSV.
#' @param out_dir Optional output directory.
#' @param top_cases Optional data frame from select_treearchi_sensitivity_top_cases().
#' @param output_level Either "focused" or "full".
#'
#' @return Invisibly returns output directory.
#' @export
redraw_treearchi_sensitivity_lineplots <- function(sensitivity_dir,
                                                   effects_file = NULL,
                                                   out_dir = NULL,
                                                   top_cases = NULL,
                                                   output_level = c("focused", "full")) {

  output_level <- match.arg(output_level)

  if (!dir.exists(sensitivity_dir)) {
    stop("sensitivity_dir does not exist.", call. = FALSE)
  }

  if (is.null(effects_file)) {
    effects_file <- file.path(sensitivity_dir, "ALL_SENSITIVITY_LOCAL_EFFECTS.csv")
  }

  if (!file.exists(effects_file)) {
    collect_treearchi_sensitivity_effects(sensitivity_dir, effects_file)
  }

  if (is.null(out_dir)) {
    out_dir <- file.path(
      sensitivity_dir,
      ifelse(
        output_level == "focused",
        "SENSITIVITY_LINEPLOTS_FOCUSED",
        "SENSITIVITY_LINEPLOTS_FULL"
      )
    )
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  dat <- utils::read.csv(effects_file, stringsAsFactors = FALSE)

  needed <- c(
    "target_group",
    "trim_var",
    "trim_tail",
    "trim_prop",
    "group_level",
    "effect_type",
    "focal_var",
    "partner_var",
    "focal_class",
    "partner_class",
    "local_beta"
  )

  missing <- setdiff(needed, names(dat))

  if (length(missing) > 0) {
    stop(
      "Missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  dat <- dat[dat$effect_type == "interaction", , drop = FALSE]
  dat <- dat[!is.na(dat$partner_var) & dat$partner_var != "", , drop = FALSE]

  if (nrow(dat) == 0) {
    stop("No interaction local effects found.", call. = FALSE)
  }

  dat$target_group <- as.character(dat$target_group)
  dat$trim_var <- as.character(dat$trim_var)
  dat$trim_tail <- as.character(dat$trim_tail)
  dat$group_level <- as.character(dat$group_level)
  dat$focal_var <- as.character(dat$focal_var)
  dat$partner_var <- as.character(dat$partner_var)
  dat$trim_prop <- suppressWarnings(as.numeric(dat$trim_prop))
  dat$local_beta <- suppressWarnings(as.numeric(dat$local_beta))

  dat <- dat[
    is.finite(dat$trim_prop) &
      is.finite(dat$local_beta),
    ,
    drop = FALSE
  ]

  if (nrow(dat) == 0) {
    stop("No finite sensitivity lineplot rows found.", call. = FALSE)
  }

  # ------------------------------------------------------------
  # Optional focused filtering
  # ------------------------------------------------------------

  if (!is.null(top_cases) && output_level == "focused") {

    needed_top <- c(
      "target_group",
      "trim_var",
      "trim_tail",
      "focal_var",
      "partner_var"
    )

    missing_top <- setdiff(needed_top, names(top_cases))

    if (length(missing_top) > 0) {
      stop(
        "top_cases is missing required column(s): ",
        paste(missing_top, collapse = ", "),
        call. = FALSE
      )
    }

    key_dat <- paste(
      dat$target_group,
      dat$trim_var,
      dat$trim_tail,
      dat$focal_var,
      dat$partner_var,
      sep = "___"
    )

    key_top <- paste(
      as.character(top_cases$target_group),
      as.character(top_cases$trim_var),
      as.character(top_cases$trim_tail),
      as.character(top_cases$focal_var),
      as.character(top_cases$partner_var),
      sep = "___"
    )

    dat <- dat[key_dat %in% key_top, , drop = FALSE]

    if (nrow(dat) == 0) {
      stop("No rows left after filtering to top_cases.", call. = FALSE)
    }
  }

  dat$trim_label <- ifelse(
    dat$trim_prop == 0,
    "Original (0%)",
    paste0(round(dat$trim_prop * 100), "% removed")
  )

  combos <- unique(dat[, c(
    "target_group",
    "trim_var",
    "trim_tail",
    "focal_var",
    "partner_var"
  )])

  index_rows <- list()

  for (i in seq_len(nrow(combos))) {

    cc <- combos[i, ]

    dd <- dat[
      dat$target_group == cc$target_group &
        dat$trim_var == cc$trim_var &
        dat$trim_tail == cc$trim_tail &
        dat$focal_var == cc$focal_var &
        dat$partner_var == cc$partner_var,
      ,
      drop = FALSE
    ]

    if (nrow(dd) == 0) next

    trim_order <- dd[order(dd$trim_prop), c("trim_prop", "trim_label")]
    trim_order <- trim_order[!duplicated(trim_order$trim_label), "trim_label"]

    dd$focal_class <- factor(dd$focal_class, levels = c("low", "medium", "high"))
    dd$partner_class <- factor(dd$partner_class, levels = c("low", "medium", "high"))
    dd$trim_label <- factor(dd$trim_label, levels = trim_order)

    dd$linewidth_val <- ifelse(
      dd$group_level == cc$target_group,
      1.35,
      0.75
    )

    dd$alpha_val <- ifelse(
      dd$group_level == cc$target_group,
      1.0,
      0.65
    )

    p <- ggplot2::ggplot(
      dd,
      ggplot2::aes(
        x = focal_class,
        y = local_beta,
        group = group_level,
        color = group_level,
        linewidth = linewidth_val,
        alpha = alpha_val
      )
    ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "grey60"
      ) +
      ggplot2::geom_line(lineend = "round") +
      ggplot2::geom_point(size = 2.3) +
      ggplot2::facet_grid(partner_class ~ trim_label) +
      ggplot2::scale_linewidth_identity() +
      ggplot2::scale_alpha_identity() +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "#E4F0E4", color = NA),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(color = "#D9E6DB", linewidth = 0.35),
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 10),
        strip.text = ggplot2::element_text(face = "bold", size = 11),
        legend.position = "right",
        plot.caption = ggplot2::element_text(size = 10, hjust = 0),
        plot.margin = ggplot2::margin(t = 12, r = 25, b = 25, l = 12)
      ) +
      ggplot2::labs(
        title = paste0("Sensitivity of effect of ", cc$focal_var),
        subtitle = paste0(
          "Target group: ", cc$target_group,
          " | trim variable: ", cc$trim_var,
          " | tail: ", cc$trim_tail,
          " | output: ", output_level
        ),
        x = paste0("Level of ", cc$focal_var),
        y = "Estimated local effect",
        color = "Group",
        caption = paste0(
          "Rows: ", cc$partner_var,
          " level. Columns: trimming proportion. Thicker line = target group."
        )
      )

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    fname <- paste0(
      "TG",
      safe_filename(cc$target_group),
      "_",
      safe_filename(cc$trim_var),
      "_",
      safe_filename(cc$trim_tail),
      "_",
      safe_filename(cc$focal_var),
      "_BY_",
      safe_filename(cc$partner_var),
      ".png"
    )

    out_path <- file.path(out_dir, fname)

    ggplot2::ggsave(
      filename = out_path,
      plot = p,
      width = 16,
      height = 8.5,
      dpi = 320,
      bg = "white",
      limitsize = FALSE
    )

    index_rows[[length(index_rows) + 1]] <- data.frame(
      target_group = cc$target_group,
      trim_var = cc$trim_var,
      trim_tail = cc$trim_tail,
      focal_var = cc$focal_var,
      partner_var = cc$partner_var,
      output_level = output_level,
      figure_png = out_path,
      stringsAsFactors = FALSE
    )
  }

  if (length(index_rows) == 0) {
    stop("No sensitivity lineplots were created.", call. = FALSE)
  }

  index_tbl <- do.call(rbind, index_rows)

  utils::write.csv(
    index_tbl,
    file.path(out_dir, "INDEX_SENSITIVITY_LINEPLOTS.csv"),
    row.names = FALSE
  )

  base::writeLines(
    c(
      "TreeArchi sensitivity lineplots",
      paste0("Created: ", Sys.time()),
      paste0("Sensitivity directory: ", sensitivity_dir),
      paste0("Effects file: ", effects_file),
      paste0("Output directory: ", out_dir),
      paste0("Output level: ", output_level),
      "",
      "Output level interpretation:",
      "- focused: lineplots are drawn only for top-ranked sensitivity cases if top_cases is provided.",
      "- full: lineplots are drawn for all sensitivity combinations in the effects file.",
      "",
      "Index:",
      "INDEX_SENSITIVITY_LINEPLOTS.csv"
    ),
    con = file.path(out_dir, "RUN_INFO.txt")
  )

  message("Sensitivity lineplots saved to: ", out_dir)

  invisible(out_dir)
}
