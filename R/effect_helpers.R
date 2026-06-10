flip_test <- function(model, dat, response, x_vars) {
  fixed_terms <- attr(stats::terms(lme4::nobars(stats::formula(model))), "term.labels")
  fixed_terms <- fixed_terms[fixed_terms != ""]
  fixed_terms <- fixed_terms[!grepl(":", fixed_terms, fixed = TRUE)]
  fixed_terms <- fixed_terms[!grepl("^I\\(.+\\^2\\)$", fixed_terms)]
  fixed_terms <- intersect(fixed_terms, x_vars)
  fixed_terms <- fixed_terms[fixed_terms %in% names(dat)]
  fixed_terms <- fixed_terms[vapply(fixed_terms, function(v) is.numeric(dat[[v]]), logical(1))]

  if (length(fixed_terms) == 0) return(tibble::tibble())

  td <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)

  out <- lapply(fixed_terms, function(pred) {
    dd <- dat %>%
      dplyr::select(dplyr::all_of(c(response, pred))) %>%
      tidyr::drop_na() %>%
      dplyr::filter(is.finite(.data[[response]]), is.finite(.data[[pred]]))

    if (nrow(dd) < 10) return(NULL)

    m_raw <- stats::lm(dd[[response]] ~ dd[[pred]])
    raw_slope <- unname(stats::coef(m_raw)[2])

    row_beta <- td %>% dplyr::filter(term == pred)
    adj_slope <- if (nrow(row_beta) > 0) unname(row_beta$estimate[1]) else NA_real_

    tibble::tibble(
      predictor = pred,
      raw_slope = raw_slope,
      adj_slope = adj_slope,
      flip = ifelse(is.finite(raw_slope) && is.finite(adj_slope), sign(raw_slope) != sign(adj_slope), NA)
    )
  })

  dplyr::bind_rows(out)
}

get_fixef_safe <- function(model, term) {
  cf <- lme4::fixef(model)
  if (term %in% names(cf)) return(unname(cf[[term]]))
  NA_real_
}

find_quad_term_exact <- function(var, coef_names) {
  cand <- paste0("I(", var, "^2)")
  if (cand %in% coef_names) return(cand)
  NA_character_
}

get_interaction_rows_for_focal <- function(model, focal_var) {
  td <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  td <- td %>% dplyr::filter(term != "(Intercept)")

  int_rows <- td %>%
    dplyr::filter(grepl(":", term, fixed = TRUE)) %>%
    dplyr::mutate(
      part1 = vapply(strsplit(term, ":", fixed = TRUE), `[`, character(1), 1),
      part2 = vapply(strsplit(term, ":", fixed = TRUE), `[`, character(1), 2)
    ) %>%
    dplyr::filter(part1 == focal_var | part2 == focal_var) %>%
    dplyr::mutate(
      partner_var = ifelse(part1 == focal_var, part2, part1)
    )

  int_rows
}

compute_marginal_effect_vector <- function(model, dat, focal_var) {
  cf <- lme4::fixef(model)
  coef_names <- names(cf)

  if (!focal_var %in% names(dat)) return(rep(NA_real_, nrow(dat)))
  if (!is.numeric(dat[[focal_var]])) return(rep(NA_real_, nrow(dat)))

  out <- rep(0, nrow(dat))

  beta_main <- if (focal_var %in% coef_names) unname(cf[[focal_var]]) else 0
  out <- out + beta_main

  int_rows <- get_interaction_rows_for_focal(model, focal_var)
  if (nrow(int_rows) > 0) {
    for (i in seq_len(nrow(int_rows))) {
      partner <- int_rows$partner_var[i]
      beta_int <- int_rows$estimate[i]
      if (partner %in% names(dat) && is.numeric(dat[[partner]])) {
        x_partner <- dat[[partner]]
        x_partner[!is.finite(x_partner)] <- NA_real_
        out <- out + beta_int * x_partner
      }
    }
  }

  quad_term <- find_quad_term_exact(focal_var, coef_names)
  if (!is.na(quad_term)) {
    beta_quad <- unname(cf[[quad_term]])
    xk <- dat[[focal_var]]
    xk[!is.finite(xk)] <- NA_real_
    out <- out + 2 * beta_quad * xk
  }

  out
}

make_3class_means <- function(x) {
  out_class <- rep(NA_character_, length(x))
  ok <- is.finite(x)

  if (sum(ok) < 3) {
    return(list(
      class = factor(out_class, levels = c("low", "medium", "high")),
      means = c(low = NA_real_, medium = NA_real_, high = NA_real_)
    ))
  }

  rank_grp <- dplyr::ntile(x[ok], 3)
  out_class[ok] <- c("low", "medium", "high")[rank_grp]
  out_class <- factor(out_class, levels = c("low", "medium", "high"))

  means <- tapply(x, out_class, mean, na.rm = TRUE)
  means <- c(
    low = unname(means["low"]),
    medium = unname(means["medium"]),
    high = unname(means["high"])
  )

  list(class = out_class, means = means)
}

compute_local_effects_table <- function(model, dat) {
  coef_names <- names(lme4::fixef(model))

  fixed_terms <- attr(stats::terms(lme4::nobars(stats::formula(model))), "term.labels")
  fixed_terms <- fixed_terms[fixed_terms != ""]
  fixed_terms <- fixed_terms[!grepl("^I\\(.+\\^2\\)$", fixed_terms)]

  main_terms_formula <- fixed_terms[!grepl(":", fixed_terms, fixed = TRUE)]
  int_parts <- coef_names[grepl(":", coef_names, fixed = TRUE)]
  int_parts <- unique(unlist(strsplit(int_parts, ":", fixed = TRUE)))
  quad_parts <- coef_names[grepl("^I\\(.+\\^2\\)$", coef_names)]
  quad_parts <- sub("^I\\((.+)\\^2\\)$", "\\1", quad_parts)

  basis_terms <- unique(c(main_terms_formula, int_parts, quad_parts))
  basis_terms <- basis_terms[basis_terms %in% names(dat)]
  basis_terms <- basis_terms[vapply(basis_terms, function(v) is.numeric(dat[[v]]), logical(1))]

  all_out <- list()

  for (focal_var in basis_terms) {
    beta_main <- get_fixef_safe(model, focal_var)
    if (!is.finite(beta_main)) beta_main <- 0

    focal_info  <- make_3class_means(dat[[focal_var]])
    focal_means <- focal_info$means

    quad_term <- find_quad_term_exact(focal_var, coef_names)
    beta_quad <- if (!is.na(quad_term)) get_fixef_safe(model, quad_term) else 0
    if (!is.finite(beta_quad)) beta_quad <- 0

    marg_vec <- compute_marginal_effect_vector(model, dat, focal_var)
    overall_local_beta <- mean(marg_vec, na.rm = TRUE)

    int_rows <- get_interaction_rows_for_focal(model, focal_var)

    if (nrow(int_rows) == 0) {
      out <- tibble::tibble(
        focal_var = focal_var,
        focal_class = c("low", "medium", "high"),
        partner_var = NA_character_,
        partner_class = NA_character_,
        focal_mean = unname(focal_means[c("low", "medium", "high")]),
        partner_mean = NA_real_,
        beta_main = beta_main,
        beta_int = NA_real_,
        beta_quad = beta_quad,
        local_beta = beta_main + 2 * beta_quad * unname(focal_means[c("low", "medium", "high")]),
        overall_local_beta = overall_local_beta
      )
      all_out[[length(all_out) + 1]] <- out
    }

    if (nrow(int_rows) > 0) {
      for (ii in seq_len(nrow(int_rows))) {
        partner_var <- int_rows$partner_var[ii]
        beta_int <- int_rows$estimate[ii]

        if (!partner_var %in% names(dat)) next
        if (!is.numeric(dat[[partner_var]])) next

        partner_info  <- make_3class_means(dat[[partner_var]])
        partner_means <- partner_info$means

        grid_df <- expand.grid(
          focal_class = c("low", "medium", "high"),
          partner_class = c("low", "medium", "high"),
          stringsAsFactors = FALSE
        ) %>%
          tibble::as_tibble() %>%
          dplyr::mutate(
            focal_var = focal_var,
            partner_var = partner_var,
            focal_mean = unname(focal_means[focal_class]),
            partner_mean = unname(partner_means[partner_class]),
            beta_main = beta_main,
            beta_int = beta_int,
            beta_quad = beta_quad,
            local_beta = beta_main + beta_int * partner_mean + 2 * beta_quad * focal_mean,
            overall_local_beta = overall_local_beta
          )

        all_out[[length(all_out) + 1]] <- grid_df
      }
    }
  }

  dplyr::bind_rows(all_out)
}

summarise_local_effects <- function(local_df) {
  if (is.null(local_df) || nrow(local_df) == 0) {
    return(tibble::tibble())
  }

  focals <- unique(local_df$focal_var)
  out <- list()

  for (focal in focals) {
    d0 <- local_df %>% dplyr::filter(focal_var == focal)
    if (nrow(d0) == 0) next

    overall_val <- unique(d0$overall_local_beta)
    overall_val <- overall_val[is.finite(overall_val)]
    overall_val <- if (length(overall_val) > 0) overall_val[1] else NA_real_

    d_self <- d0 %>%
      dplyr::mutate(
        partner_rank = dplyr::case_when(
          is.na(partner_class) ~ 1L,
          partner_class == "medium" ~ 1L,
          TRUE ~ 2L
        )
      ) %>%
      dplyr::arrange(partner_rank) %>%
      dplyr::group_by(focal_class) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

    low_val  <- if ("low" %in% d_self$focal_class) d_self$local_beta[d_self$focal_class == "low"][1] else NA_real_
    mid_val  <- if ("medium" %in% d_self$focal_class) d_self$local_beta[d_self$focal_class == "medium"][1] else NA_real_
    high_val <- if ("high" %in% d_self$focal_class) d_self$local_beta[d_self$focal_class == "high"][1] else NA_real_

    out[[length(out) + 1]] <- tibble::tibble(
      term = focal,
      overall_local_beta = overall_val,
      low_local_beta = low_val,
      mid_local_beta = mid_val,
      high_local_beta = high_val
    )
  }

  dplyr::bind_rows(out)
}
