`%>%` <- magrittr::`%>%`

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
})

fit_lmer <- function(formula_obj, dat) {
  lme4::lmer(
    formula_obj,
    data = dat,
    REML = FALSE,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5),
      check.scaleX = "warning"
    )
  )
}

make_formula <- function(resp, fixed_terms, rand_terms) {
  fixed_terms <- fixed_terms[fixed_terms != ""]
  fixed_part <- if (length(fixed_terms) == 0) "1" else paste(fixed_terms, collapse = " + ")
  stats::as.formula(paste(resp, "~", fixed_part, "+", rand_terms))
}

mk_sq <- function(x) {
  paste0("I(", x, "^2)")
}

build_interaction_terms <- function(x_vars) {
  if (length(x_vars) < 2) return(character(0))
  combn(x_vars, 2, FUN = function(z) paste(z, collapse = ":"), simplify = TRUE)
}

build_int_quad_terms <- function(x_vars) {
  lin_terms  <- x_vars
  int_terms  <- build_interaction_terms(x_vars)
  quad_terms <- mk_sq(x_vars)
  unique(c(lin_terms, int_terms, quad_terms))
}

backward_select_lmm <- function(dat, resp, forced, candidates, rand_terms,
                                criterion = c("BIC"),
                                independent = FALSE,
                                tol = 1e-9) {

  criterion <- match.arg(criterion)

  current_terms <- unique(c(forced, candidates))
  current_form  <- make_formula(resp, current_terms, rand_terms)
  current_mod   <- fit_lmer(current_form, dat)
  current_ic_v  <- stats::BIC(current_mod)

  history <- list()
  step <- 0
  improved <- TRUE

  while (improved) {
    improved <- FALSE
    step <- step + 1

    drop_candidates <- setdiff(current_terms, forced)
    if (length(drop_candidates) == 0) break

    trials <- purrr::map_dfr(drop_candidates, function(term_to_drop) {
      trial_terms <- setdiff(current_terms, term_to_drop)

      if (!independent) {
        ints <- trial_terms[grepl(":", trial_terms, fixed = TRUE)]
        if (length(ints) > 0) {
          mains_needed <- unique(unlist(strsplit(ints, ":", fixed = TRUE)))
          if (term_to_drop %in% mains_needed) {
            return(tibble::tibble(drop = term_to_drop, ok = FALSE, crit = NA_real_))
          }
        }
      }

      trial_form <- make_formula(resp, trial_terms, rand_terms)
      trial_mod  <- tryCatch(fit_lmer(trial_form, dat), error = function(e) NULL)
      if (is.null(trial_mod)) {
        return(tibble::tibble(drop = term_to_drop, ok = FALSE, crit = NA_real_))
      }

      crit_v <- stats::BIC(trial_mod)
      tibble::tibble(drop = term_to_drop, ok = TRUE, crit = crit_v, model = list(trial_mod))
    })

    trials_ok <- dplyr::arrange(
      dplyr::filter(trials, ok),
      crit
    )
    if (nrow(trials_ok) == 0) break

    best <- dplyr::slice(trials_ok, 1)

    if (is.finite(best$crit) && best$crit < (current_ic_v - tol)) {
      term_drop     <- best$drop
      current_terms <- setdiff(current_terms, term_drop)
      current_mod   <- best$model[[1]]
      current_ic_v  <- best$crit

      history[[length(history) + 1]] <- tibble::tibble(
        step = step,
        dropped = term_drop,
        criterion = "BIC",
        crit_value = current_ic_v,
        terms_remaining = paste(current_terms, collapse = " + ")
      )

      improved <- TRUE
    }
  }

  hist_df <- if (length(history) == 0) tibble::tibble() else dplyr::bind_rows(history)

  list(
    model = current_mod,
    terms = current_terms,
    history = hist_df,
    BIC = current_ic_v
  )
}
