#' Collect TreeArchi sensitivity local effects
#'
#' @param sensitivity_dir Sensitivity output directory.
#' @param out_file Optional output CSV path.
#'
#' @return Data frame of collected local effects.
#' @export
collect_treearchi_sensitivity_effects <- function(
    sensitivity_dir,
    out_file = NULL
) {

  if (!dir.exists(sensitivity_dir)) {
    stop("sensitivity_dir does not exist.", call. = FALSE)
  }

  effect_files <- list.files(
    sensitivity_dir,
    pattern = "LOCAL_EFFECTS.*\\.csv$|ALL_LOCAL_EFFECTS.*\\.csv$|ALL_ANALYTIC_MARGINAL_EFFECTS.*\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  effect_files <- effect_files[
    !grepl("ALL_SENSITIVITY_LOCAL_EFFECTS", basename(effect_files))
  ]

  read_one <- function(path) {

    x <- tryCatch(
      utils::read.csv(path, stringsAsFactors = FALSE),
      error = function(e) NULL
    )

    if (is.null(x)) return(NULL)

    needed <- c("focal_var", "local_beta")

    if (!all(needed %in% names(x))) {
      return(NULL)
    }

    x$source_file <- path

    # --------------------------------------------------
    # Extract TGshorter_branch_len_l10 style folder
    # --------------------------------------------------

    path_parts <- strsplit(
      normalizePath(path, winslash = "/", mustWork = FALSE),
      "/"
    )[[1]]

    tg_dirs <- path_parts[
      grepl("^TG(shorter|taller)_", path_parts)
    ]

    if (length(tg_dirs) >= 1) {

      base_path <- tg_dirs[1]

      parts <- strsplit(base_path, "_")[[1]]

      target_group_value <- sub("^TG", "", parts[1])

      trim_tail_prop <- parts[length(parts)]

      trim_tail_code <- substr(trim_tail_prop, 1, 1)

      trim_tail_value <- ifelse(
        trim_tail_code == "l",
        "lower",
        "upper"
      )

      trim_prop_value <- as.numeric(
        sub("^[lu]", "", trim_tail_prop)
      ) / 100

      trim_var_value <- paste(
        parts[2:(length(parts) - 1)],
        collapse = "_"
      )

    } else {

      target_group_value <- NA_character_
      trim_var_value <- NA_character_
      trim_tail_value <- NA_character_
      trim_prop_value <- NA_real_
    }

    x$target_group <- target_group_value
    x$trim_var <- trim_var_value
    x$trim_tail <- trim_tail_value
    x$trim_prop <- trim_prop_value

    x
  }

  lst <- lapply(effect_files, read_one)

  lst <- lst[
    !vapply(lst, is.null, logical(1))
  ]

  if (length(lst) == 0) {
    stop(
      "No usable local effect CSV files found.",
      call. = FALSE
    )
  }

  all_names <- unique(
    unlist(
      lapply(lst, names)
    )
  )

  align_cols <- function(x) {

    missing <- setdiff(
      all_names,
      names(x)
    )

    for (m in missing) {
      x[[m]] <- NA
    }

    x[, all_names, drop = FALSE]
  }

  out <- do.call(
    rbind,
    lapply(lst, align_cols)
  )

  if (is.null(out_file)) {
    out_file <- file.path(
      sensitivity_dir,
      "ALL_SENSITIVITY_LOCAL_EFFECTS.csv"
    )
  }

  utils::write.csv(
    out,
    out_file,
    row.names = FALSE
  )

  message(
    "Collected local effects: ",
    nrow(out)
  )

  message(
    "Saved to: ",
    out_file
  )

  out
}
