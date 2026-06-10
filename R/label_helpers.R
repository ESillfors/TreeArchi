label_term <- function(term, labels_map = pretty_names) {
  if (term == "alpha_volume_m3") return(paste0(labels_map[[term]] %||% term, " [3D]"))
  if (term == "projected_area_m2") return(paste0(labels_map[[term]] %||% term, " [2D]"))
  labels_map[[term]] %||% term
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
    is.na(p), NA_character_,
    ifelse(p < 0.001, "<0.001", sub("^0", "", formatC(p, digits = 3, format = "f")))
  )
}

sign_to_arrow <- function(est) {
  ifelse(is.na(est), "(?)", ifelse(est > 0, "-> (+)", ifelse(est < 0, "-> (-)", "-> (0)")))
}

fmt_beta <- function(x) {
  ifelse(is.finite(x), sprintf("%.3f", x), "NA")
}

pretty_term <- function(term, labels_map = pretty_names) {
  vapply(term, function(tt) {
    if (tt %in% names(labels_map)) {
      base <- unname(labels_map[tt])
      if (tt == "alpha_volume_m3") base <- paste0(base, " [3D]")
      if (tt == "projected_area_m2") base <- paste0(base, " [2D]")
      return(base)
    }

    if (grepl("^I\\(.+\\^2\\)$", tt)) {
      base <- gsub("^I\\((.+)\\^2\\)$", "\\1", tt)
      base_pretty <- if (base %in% names(labels_map)) unname(labels_map[base]) else base
      if (base == "alpha_volume_m3") base_pretty <- paste0(base_pretty, " [3D]")
      if (base == "projected_area_m2") base_pretty <- paste0(base_pretty, " [2D]")
      return(paste0(base_pretty, " [Sq]"))
    }

    if (grepl(":", tt, fixed = TRUE)) {
      parts <- unlist(strsplit(tt, ":", fixed = TRUE))
      parts_pretty <- vapply(parts, function(z) pretty_term(z, labels_map), character(1))
      return(paste(parts_pretty, collapse = " × "))
    }

    tt
  }, character(1))
}

make_path_label <- function(pred_term, estimate, response_label, labels_map = pretty_names) {
  pred_pretty <- pretty_term(pred_term, labels_map)[1]
  pred_pretty <- gsub(" \\[3D\\]", "", pred_pretty)
  pred_pretty <- gsub(" \\[2D\\]", "", pred_pretty)
  resp_clean  <- gsub(" \\[3D\\]", "", response_label)
  resp_clean  <- gsub(" \\[2D\\]", "", response_label)
  paste0(pred_pretty, " ", sign_to_arrow(estimate), " ", resp_clean)
}

short_pretty <- function(x) {
  x <- gsub(" \\[3D\\]", "", x)
  x <- gsub(" \\[2D\\]", "", x)
  x
}
