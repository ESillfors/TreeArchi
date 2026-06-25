#' Select TreeArchi convergence cases
#'
#' Finds trimming scenarios where a target group becomes most similar to a
#' user-defined reference profile. If several reference groups are provided,
#' the reference profile is calculated as the mean local-effect profile across
#' those groups.
#'
#' Ranking is primarily based on the smallest distance_after, meaning the
#' trimmed target-group profile is closest to the reference profile.
#' Convergence is also reported as distance_before - distance_after.
#'
#' @param effects_tbl Data frame returned by collect_treearchi_sensitivity_effects().
#' @param effects_file Optional path to ALL_SENSITIVITY_LOCAL_EFFECTS.csv.
#' @param out_dir Output directory for convergence tables and plot.
#' @param target_group Target group being trimmed.
#' @param reference_group One or more reference groups used to form the reference profile.
#' @param top_n Number of top convergence cases per interaction.
#' @param exclude_trim_vars Trim variables excluded from ranking.
#'
#' @return A data frame of top convergence cases.
#' @export
select_treearchi_convergence_cases <- function(
    effects_tbl = NULL,
    effects_file = NULL,
    out_dir,
    target_group,
    reference_group,
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

  target_group <- as.character(target_group)
  reference_group <- as.character(reference_group)

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

  if (length(target_group) != 1) {
    stop("target_group must contain exactly one group.", call. = FALSE)
  }

  if (!(target_group %in% unique(dat$group_level))) {
    stop("target_group was not found in group_level.", call. = FALSE)
  }

  missing_ref <- setdiff(reference_group, unique(dat$group_level))

  if (length(missing_ref) > 0) {
    stop(
      "reference_group value(s) not found in group_level: ",
      paste(missing_ref, collapse = ", "),
      call. = FALSE
    )
  }

  dat$cell_id <- paste(
    dat$focal_class,
    dat$partner_class,
    sep = "_"
  )

  dat$interaction <- paste0(
    dat$focal_var,
    "_BY_",
    dat$partner_var
  )

  target_dat <- dat[
    dat$group_level == target_group,
    ,
    drop = FALSE
  ]

  ref_dat <- dat[
    dat$group_level %in% reference_group &
      dat$trim_prop == 0,
    ,
    drop = FALSE
  ]

  if (nrow(target_dat) == 0 || nrow(ref_dat) == 0) {
    stop("No usable target or reference profile rows found.", call. = FALSE)
  }

  # ------------------------------------------------------------
  # Reference profile
  # ------------------------------------------------------------

  ref_profile <- stats::aggregate(
    ref_dat$local_beta,
    by = list(
      interaction = ref_dat$interaction,
      focal_var = ref_dat$focal_var,
      partner_var = ref_dat$partner_var,
      cell_id = ref_dat$cell_id
    ),
    FUN = function(z) mean(z, na.rm = TRUE)
  )

  names(ref_profile)[names(ref_profile) == "x"] <- "reference_beta"

  ref_profile$n_reference_groups <- length(unique(reference_group))
  ref_profile$reference_group <- paste(reference_group, collapse = ", ")

  scenarios <- unique(
    target_dat[, c(
      "trim_var",
      "trim_tail",
      "trim_prop",
      "interaction",
      "focal_var",
      "partner_var"
    )]
  )

  scenarios <- scenarios[
    scenarios$trim_prop > 0,
    ,
    drop = FALSE
  ]

  if (nrow(scenarios) == 0) {
    stop("No non-zero trimming scenarios found for target_group.", call. = FALSE)
  }

  result_rows <- list()
  rr <- 0

  for (ii in seq_len(nrow(scenarios))) {

    ss <- scenarios[ii, ]

    target_after <- target_dat[
      target_dat$trim_var == ss$trim_var &
        target_dat$trim_tail == ss$trim_tail &
        target_dat$trim_prop == ss$trim_prop &
        target_dat$focal_var == ss$focal_var &
        target_dat$partner_var == ss$partner_var,
      ,
      drop = FALSE
    ]

    target_before <- target_dat[
      target_dat$trim_prop == 0 &
        target_dat$focal_var == ss$focal_var &
        target_dat$partner_var == ss$partner_var,
      ,
      drop = FALSE
    ]

    ref_i <- ref_profile[
      ref_profile$focal_var == ss$focal_var &
        ref_profile$partner_var == ss$partner_var,
      ,
      drop = FALSE
    ]

    if (nrow(target_after) == 0 ||
        nrow(target_before) == 0 ||
        nrow(ref_i) == 0) {
      next
    }

    before_tbl <- merge(
      target_before[, c("cell_id", "local_beta")],
      ref_i[, c("cell_id", "reference_beta")],
      by = "cell_id"
    )

    after_tbl <- merge(
      target_after[, c("cell_id", "local_beta")],
      ref_i[, c("cell_id", "reference_beta")],
      by = "cell_id"
    )

    if (nrow(before_tbl) == 0 || nrow(after_tbl) == 0) next

    distance_before <- sqrt(
      mean(
        (before_tbl$local_beta - before_tbl$reference_beta)^2,
        na.rm = TRUE
      )
    )

    distance_after <- sqrt(
      mean(
        (after_tbl$local_beta - after_tbl$reference_beta)^2,
        na.rm = TRUE
      )
    )

    convergence <- distance_before - distance_after

    rr <- rr + 1

    result_rows[[rr]] <- data.frame(
      target_group = target_group,
      reference_group = paste(reference_group, collapse = ", "),
      n_reference_groups = length(reference_group),
      trim_var = ss$trim_var,
      trim_tail = ss$trim_tail,
      trim_prop = ss$trim_prop,
      interaction = ss$interaction,
      focal_var = ss$focal_var,
      partner_var = ss$partner_var,
      distance_before = distance_before,
      distance_after = distance_after,
      convergence = convergence,
      stringsAsFactors = FALSE
    )
  }

  if (length(result_rows) == 0) {
    stop("No convergence results could be calculated.", call. = FALSE)
  }

  ranking_all <- do.call(rbind, result_rows)

  # ------------------------------------------------------------
  # Rank convergence cases
  # ------------------------------------------------------------

  ranking_all <- ranking_all[
    order(
      ranking_all$interaction,
      ranking_all$distance_after,
      -ranking_all$convergence
    ),
    ,
    drop = FALSE
  ]

  dedup_key <- paste(
    ranking_all$interaction,
    ranking_all$trim_var,
    ranking_all$trim_tail,
    sep = "___"
  )

  ranking_dedup <- ranking_all[!duplicated(dedup_key), , drop = FALSE]

  split_by_interaction <- split(
    ranking_dedup,
    ranking_dedup$interaction
  )

  top_tbl <- do.call(
    rbind,
    lapply(split_by_interaction, function(x) {
      x <- x[
        order(
          x$distance_after,
          -x$convergence
        ),
        ,
        drop = FALSE
      ]

      utils::head(x, top_n)
    })
  )

  rownames(top_tbl) <- NULL

  top_tbl <- top_tbl[
    order(
      top_tbl$interaction,
      top_tbl$distance_after,
      -top_tbl$convergence
    ),
    ,
    drop = FALSE
  ]

  top_tbl$convergence_rank <- ave(
    top_tbl$distance_after,
    top_tbl$interaction,
    FUN = function(z) rank(z, ties.method = "first")
  )

  utils::write.csv(
    ref_profile,
    file.path(out_dir, "CONVERGENCE_REFERENCE_PROFILE.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    ranking_all,
    file.path(out_dir, "CONVERGENCE_RANKING_ALL.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    ranking_dedup,
    file.path(out_dir, "CONVERGENCE_RANKING_DEDUP.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    top_tbl,
    file.path(out_dir, paste0("TOP", top_n, "_CONVERGENCE_CASES.csv")),
    row.names = FALSE
  )

  # ------------------------------------------------------------
  # Barplot: top convergence cases per interaction
  # ------------------------------------------------------------

  plot_tbl <- top_tbl

  plot_tbl$trim_percent <- round(plot_tbl$trim_prop * 100)

  plot_tbl$interaction_label <- paste0(
    pretty_var_name(plot_tbl$focal_var),
    " × ",
    pretty_var_name(plot_tbl$partner_var)
  )

  plot_tbl$case_label <- make_trim_case_label(
    trim_var = plot_tbl$trim_var,
    trim_tail = plot_tbl$trim_tail,
    trim_prop = plot_tbl$trim_prop
  )

  plot_tbl$case_label_panel <- paste0(
    plot_tbl$interaction_label,
    "___",
    plot_tbl$case_label
  )

  plot_tbl$case_label_panel <- stats::reorder(
    plot_tbl$case_label_panel,
    plot_tbl$distance_after,
    FUN = function(x) -mean(x, na.rm = TRUE)
  )

  p_top <- ggplot2::ggplot(
    plot_tbl,
    ggplot2::aes(
      x = distance_after,
      y = case_label_panel,
      fill = convergence
    )
  ) +
    ggplot2::geom_col(width = 0.70) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = paste0(
          "d = ",
          round(distance_after, 3),
          " | Δ = ",
          round(convergence, 3)
        )
      ),
      hjust = -0.05,
      size = 3.3,
      fontface = "bold"
    ) +
    ggplot2::facet_wrap(
      ~ interaction_label,
      scales = "free_y",
      ncol = 1
    ) +
    ggplot2::scale_y_discrete(
      labels = function(x) sub("^.*___", "", x)
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#B22234",
      mid = "grey82",
      high = "#18B7BE",
      midpoint = 0
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.25))
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 17),
      plot.subtitle = ggplot2::element_text(size = 12),
      strip.text = ggplot2::element_text(face = "bold", size = 12),
      axis.text.y = ggplot2::element_text(size = 10, face = "bold"),
      axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(
      title = paste0("Top ", top_n, " convergence cases per interaction"),
      subtitle = paste0(
        "Target group: ",
        target_group,
        " → reference profile: ",
        paste(reference_group, collapse = ", "),
        ". Smaller distance means closer match."
      ),
      x = "Distance to reference profile after trimming",
      y = NULL,
      fill = "Convergence\n(before − after)"
    )

  ggplot2::ggsave(
    filename = file.path(out_dir, paste0("TOP", top_n, "_CONVERGENCE_CASES.png")),
    plot = p_top,
    width = 15,
    height = max(8, 2.8 * length(unique(plot_tbl$interaction_label))),
    dpi = 320,
    bg = "white",
    limitsize = FALSE
  )

  base::writeLines(
    c(
      "TreeArchi convergence ranking",
      paste0("Created: ", Sys.time()),
      paste0("Target group: ", target_group),
      paste0("Reference group(s): ", paste(reference_group, collapse = ", ")),
      paste0("Reference profile: mean local-effect profile across ", length(reference_group), " reference group(s)."),
      "",
      "Ranking:",
      "- Primary ranking is based on the smallest distance_after.",
      "- distance_after is the RMS distance between the trimmed target-group local-effect profile and the reference profile.",
      "- convergence = distance_before - distance_after.",
      "- Positive convergence means the target group moved closer to the reference profile after trimming.",
      "",
      "Output files:",
      "CONVERGENCE_REFERENCE_PROFILE.csv",
      "CONVERGENCE_RANKING_ALL.csv",
      "CONVERGENCE_RANKING_DEDUP.csv",
      paste0("TOP", top_n, "_CONVERGENCE_CASES.csv"),
      paste0("TOP", top_n, "_CONVERGENCE_CASES.png")
    ),
    con = file.path(out_dir, "CONVERGENCE_RUN_INFO.txt")
  )

  message("Convergence ranking saved to: ", out_dir)

  top_tbl
}
