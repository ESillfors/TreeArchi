# ------------------------------------------------------------
# TreeArchi label helpers
# ------------------------------------------------------------

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

pretty_names <- c(

  # Shape traits
  sba_degrees = "Stem Branch angle (°)",
  ba2 = "Second Branch angle (°)",
  bar = "Branch angle ratio",
  cdhr = "Crown spread",
  rvr = "Relative volume ratio",
  dhr = "DBH / height",
  sbd = "Stem–branch distance (m)",
  ch_raw = "Crown depth (m)",
  csh_raw = "Crown start height (m)",

  # Size traits
  sbr = "Stem–branch radius (m)",
  sbl = "Stem–branch length (m)",
  clvr_m.2 = "Branch slenderness (m²)",
  branch_len = "Total branch length (m)",
  branch_vol_m3 = "Total Branch volume (m³)",
  cd_raw = "Crown diameter (m)",
  projected_area_m2 = "Projected crown area (m²)",
  alpha_volume_m3 = "Alpha volume (m³)",
  dbh_m = "DBH (m)",
  tree_height_m = "Tree height (m)",
  tree_vol_m3 = "Tree volume (m³)",
  trunk_vol_m3 = "Trunk volume (m³)",
  base_vol_0_10 = "Base volume (0–10%)",
  AGB_TLS = "AGB (TLS)",

  # Common identifiers / grouping variables
  tls_id = "TLS ID",
  plot = "Plot",
  plot_factor = "Plot",
  genus = "Genus",
  species = "Species",
  family = "Family"
)

label_term <- function(term, labels_map = pretty_names) {
  vapply(
    as.character(term),
    function(tt) {
      base <- if (!is.null(labels_map) && tt %in% names(labels_map)) unname(labels_map[tt]) else tt

      if (tt == "alpha_volume_m3") {
        base <- paste0(base, " [3D]")
      }

      if (tt == "projected_area_m2") {
        base <- paste0(base, " [2D]")
      }

      base
    },
    character(1)
  )
}

pretty_var_name <- function(x, labels_map = pretty_names) {
  vapply(
    as.character(x),
    function(xx) {
      label_term(xx, labels_map)
    },
    character(1)
  )
}

wrap_words <- function(x, n_words = 6) {
  if (is.na(x) || !nzchar(x)) return(x)

  w <- unlist(strsplit(x, "\\s+"))

  if (length(w) <= n_words) return(x)

  idx <- seq(n_words + 1, length(w), by = n_words)
  w[idx] <- paste0("\n", w[idx])

  paste(w, collapse = " ")
}

safe_p <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(
      p < 0.001,
      "<0.001",
      sub("^0", "", formatC(p, digits = 3, format = "f"))
    )
  )
}

sign_to_arrow <- function(est) {
  ifelse(
    is.na(est),
    "(?)",
    ifelse(
      est > 0,
      "-> (+)",
      ifelse(est < 0, "-> (-)", "-> (0)")
    )
  )
}

fmt_beta <- function(x) {
  ifelse(is.finite(x), sprintf("%.3f", x), "NA")
}

pretty_term <- function(term, labels_map = pretty_names) {
  vapply(
    term,
    function(tt) {

      if (tt %in% names(labels_map)) {
        base <- unname(labels_map[tt])

        if (tt == "alpha_volume_m3") {
          base <- paste0(base, " [3D]")
        }

        if (tt == "projected_area_m2") {
          base <- paste0(base, " [2D]")
        }

        return(base)
      }

      if (grepl("^I\\(.+\\^2\\)$", tt)) {
        base <- gsub("^I\\((.+)\\^2\\)$", "\\1", tt)

        base_pretty <- if (base %in% names(labels_map)) {
          unname(labels_map[base])
        } else {
          base
        }

        if (base == "alpha_volume_m3") {
          base_pretty <- paste0(base_pretty, " [3D]")
        }

        if (base == "projected_area_m2") {
          base_pretty <- paste0(base_pretty, " [2D]")
        }

        return(paste0(base_pretty, " [Sq]"))
      }

      if (grepl(":", tt, fixed = TRUE)) {
        parts <- unlist(strsplit(tt, ":", fixed = TRUE))
        parts_pretty <- vapply(
          parts,
          function(z) pretty_term(z, labels_map),
          character(1)
        )

        return(paste(parts_pretty, collapse = " × "))
      }

      tt
    },
    character(1)
  )
}

make_path_label <- function(pred_term, estimate, response_label, labels_map = pretty_names) {
  pred_pretty <- pretty_term(pred_term, labels_map)[1]

  pred_pretty <- gsub(" \\[3D\\]", "", pred_pretty)
  pred_pretty <- gsub(" \\[2D\\]", "", pred_pretty)

  resp_clean <- gsub(" \\[3D\\]", "", response_label)
  resp_clean <- gsub(" \\[2D\\]", "", resp_clean)

  paste0(pred_pretty, " ", sign_to_arrow(estimate), " ", resp_clean)
}

short_pretty <- function(x) {
  x <- gsub(" \\[3D\\]", "", x)
  x <- gsub(" \\[2D\\]", "", x)
  x
}

trim_tail_label <- function(trim_tail) {
  trim_tail <- as.character(trim_tail)

  ifelse(
    trim_tail == "lower",
    "Lowest",
    ifelse(trim_tail == "upper", "Highest", trim_tail)
  )
}

make_trim_case_label <- function(trim_var, trim_tail, trim_prop, labels_map = pretty_names) {
  trim_percent <- round(as.numeric(trim_prop) * 100)

  paste0(
    trim_tail_label(trim_tail),
    " ",
    trim_percent,
    "% removed: ",
    pretty_var_name(trim_var, labels_map)
  )
}
