#' Select top TreeArchi sensitivity cases
#'
#' Selects the strongest sensitivity cases from collected group-level
#' local effects. This is intended to be used after
#' collect_treearchi_sensitivity_effects().
#'
#' @param effects_tbl Data frame returned by collect_treearchi_sensitivity_effects().
#' @param effects_file Optional path to ALL_SENSITIVITY_LOCAL_EFFECTS.csv.
#' @param out_dir Output directory for top-case tables.
#' @param top_n Number of top cases per target group x interaction.
#' @param exclude_trim_vars Trim variables excluded from ranking.
#'
#' @return A data frame of selected top sensitivity cases.
#' @export
select_treearchi_sensitivity_top_cases <- function(
    effects_tbl = NULL,
    effects_file = NULL,
    out_dir,
    top_n = 5,
    exclude_trim_vars = c("projected_area_m2", "cd_raw")
) {

  if (is.null(effects_tbl)) {
    if (is.null(effects_file) || !file.exists(effects_file)) {
      stop("Provide either effects_tbl or an existing effects_file.", call. = FALSE)
    }

    effects_tbl <- utils::read.csv(
      effects_file,
      stringsAsFactors = FALSE
    )
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

  missing <- setdiff(needed, names(effects_tbl))

  if (length(missing) > 0) {
    stop(
      "Missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  clean_text <- function(x) {
    x <- as.character(x)
    x <- tolower(x)
    x <- gsub("^log_", "", x)
    x <- gsub("_", " ", x)
    x <- gsub("\\.", " ", x)
    x <- gsub("\\([^\\)]*\\)", "", x)
    x <- gsub("\\[[^\\]]*\\]", "", x)
    x <- gsub("\\s+", " ", x)
    trimws(x)
  }

  is_forbidden_trim_var <- function(x) {
    x_clean <- clean_text(x)
    x_raw <- as.character(x)

    x_raw %in% exclude_trim_vars |
      grepl("projected", x_clean) |
      grepl("^cd raw$|^cd$|crown diameter", x_clean)
  }

  dat <- effects_tbl

  dat$target_group <- as.character(dat$target_group)
  dat$trim_var <- as.character(dat$trim_var)
  dat$trim_tail <- as.character(dat$trim_tail)
  dat$group_level <- as.character(dat$group_level)
  dat$effect_type <- as.character(dat$effect_type)
  dat$focal_var <- as.character(dat$focal_var)
  dat$partner_var <- as.character(dat$partner_var)
  dat$focal_class <- as.character(dat$focal_class)
  dat$partner_class <- as.character(dat$partner_class)
  dat$trim_prop <- suppressWarnings(as.numeric(dat$trim_prop))
  dat$local_beta <- suppressWarnings(as.numeric(dat$local_beta))

  dat <- dat[
    dat$effect_type == "interaction" &
      !is.na(dat$partner_var) &
      dat$partner_var != "" &
      is.finite(dat$local_beta) &
      is.finite(dat$trim_prop) &
      !is_forbidden_trim_var(dat$trim_var),
    ,
    drop = FALSE
  ]

  if (nrow(dat) == 0) {
    stop("No usable interaction sensitivity effects found.", call. = FALSE)
  }

  dat$interaction <- paste0(dat$focal_var, "_BY_", dat$partner_var)

  key_cols <- c(
    "target_group",
    "trim_var",
    "trim_tail",
    "interaction",
    "focal_var",
    "partner_var",
    "group_level",
    "focal_class",
    "partner_class"
  )

  baseline <- dat[
    dat$trim_prop == 0,
    c(key_cols, "local_beta"),
    drop = FALSE
  ]

  names(baseline)[names(baseline) == "local_beta"] <- "baseline_beta"

  dat2 <- merge(
    dat,
    baseline,
    by = key_cols,
    all.x = TRUE
  )

  dat2$delta_from_original <- dat2$local_beta - dat2$baseline_beta

  dat2 <- dat2[
    dat2$trim_prop > 0 &
      is.finite(dat2$delta_from_original),
    ,
    drop = FALSE
  ]

  if (nrow(dat2) == 0) {
    stop(
      "No non-zero trimming cases with matching original baseline were found.",
      call. = FALSE
    )
  }

  # Rank one trim variable + tail scenario at a time.
  # The score describes how strongly the target-group local effect changes
  # across trimming levels compared with the original 0% run.
  target_rows <- dat2[
    dat2$group_level == dat2$target_group,
    ,
    drop = FALSE
  ]

  summary_tbl <- stats::aggregate(
    abs(target_rows$delta_from_original),
    by = list(
      target_group = target_rows$target_group,
      interaction = target_rows$interaction,
      focal_var = target_rows$focal_var,
      partner_var = target_rows$partner_var,
      trim_var = target_rows$trim_var,
      trim_tail = target_rows$trim_tail
    ),
    FUN = function(z) {
      mean(z, na.rm = TRUE)
    }
  )

  names(summary_tbl)[names(summary_tbl) == "x"] <- "mean_abs_target_change"

  max_tbl <- stats::aggregate(
    abs(target_rows$delta_from_original),
    by = list(
      target_group = target_rows$target_group,
      interaction = target_rows$interaction,
      focal_var = target_rows$focal_var,
      partner_var = target_rows$partner_var,
      trim_var = target_rows$trim_var,
      trim_tail = target_rows$trim_tail
    ),
    FUN = function(z) {
      max(z, na.rm = TRUE)
    }
  )

  names(max_tbl)[names(max_tbl) == "x"] <- "max_abs_target_change"

  summary_tbl <- merge(
    summary_tbl,
    max_tbl,
    by = c(
      "target_group",
      "interaction",
      "focal_var",
      "partner_var",
      "trim_var",
      "trim_tail"
    ),
    all.x = TRUE
  )

  summary_tbl <- summary_tbl[order(
    summary_tbl$target_group,
    summary_tbl$interaction,
    -summary_tbl$mean_abs_target_change,
    -summary_tbl$max_abs_target_change
  ), ]

  # Same trim variable only once per tail per target group x interaction.
  summary_tbl <- summary_tbl[
    !duplicated(summary_tbl[, c(
      "target_group",
      "interaction",
      "trim_var",
      "trim_tail"
    )]),
    ,
    drop = FALSE
  ]

  split_id <- paste(summary_tbl$target_group, summary_tbl$interaction, sep = "___")
  split_tbl <- split(summary_tbl, split_id)

  top_tbl <- do.call(
    rbind,
    lapply(split_tbl, function(x) {
      x <- x[order(
        -x$mean_abs_target_change,
        -x$max_abs_target_change
      ), , drop = FALSE]

      utils::head(x, top_n)
    })
  )

  rownames(top_tbl) <- NULL

  top_tbl$top_rank <- ave(
    -top_tbl$mean_abs_target_change,
    top_tbl$target_group,
    top_tbl$interaction,
    FUN = function(z) rank(z, ties.method = "first")
  )

  top_tbl <- top_tbl[order(
    top_tbl$target_group,
    top_tbl$interaction,
    top_tbl$top_rank
  ), ]

  utils::write.csv(
    summary_tbl,
    file.path(out_dir, "SENSITIVITY_CASE_RANKING_ALL.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    top_tbl,
    file.path(out_dir, paste0("TOP", top_n, "_SENSITIVITY_CASES.csv")),
    row.names = FALSE
  )

  # ------------------------------------------------------------
  # Top-case barplot
  # ------------------------------------------------------------

  plot_tbl <- top_tbl

  plot_tbl$case_label <- paste0(
    "TG", plot_tbl$target_group,
    " | ",
    plot_tbl$trim_var,
    " | ",
    plot_tbl$trim_tail,
    " | ",
    plot_tbl$focal_var,
    " × ",
    plot_tbl$partner_var
  )

  plot_tbl$case_label <- stats::reorder(
    plot_tbl$case_label,
    plot_tbl$mean_abs_target_change
  )

  p_top <- ggplot2::ggplot(
    plot_tbl,
    ggplot2::aes(
      x = mean_abs_target_change,
      y = case_label,
      fill = trim_tail
    )
  ) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::facet_wrap(
      ~ interaction,
      scales = "free_y",
      ncol = 1
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(size = 11),
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      axis.text.y = ggplot2::element_text(size = 9),
      legend.position = "bottom"
    ) +
    ggplot2::labs(
      title = paste0("Top ", top_n, " TreeArchi sensitivity cases"),
      subtitle = "Cases ranked by mean absolute change in the target-group local effect",
      x = "Mean absolute target-group change",
      y = NULL,
      fill = "Trim tail"
    )

  ggplot2::ggsave(
    filename = file.path(out_dir, paste0("TOP", top_n, "_SENSITIVITY_CASES.png")),
    plot = p_top,
    width = 13,
    height = max(6, 0.42 * nrow(plot_tbl)),
    dpi = 320,
    bg = "white",
    limitsize = FALSE
  )

  message("Top sensitivity cases saved to: ", out_dir)

  top_tbl
}
