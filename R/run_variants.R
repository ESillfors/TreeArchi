run_pure_bic <- function(dat0, response, forced_var, x_vars, rand_terms) {
  form_full <- make_formula(
    response,
    fixed_terms = unique(c(forced_var, setdiff(x_vars, forced_var))),
    rand_terms = rand_terms
  )
  m_full_pure <- fit_lmer(form_full, dat0)

  sel_pure <- backward_select_lmm(
    dat0,
    response,
    forced = forced_var,
    candidates = setdiff(x_vars, forced_var),
    rand_terms = rand_terms,
    independent = FALSE
  )

  list(
    full_model = m_full_pure,
    selected = sel_pure
  )
}

run_interaction_bic <- function(datC, response, forced_var, x_vars, rand_terms) {
  form_int_full <- stats::as.formula(
    paste0(response, " ~ (", paste(x_vars, collapse = " + "), ")^2 + ", rand_terms)
  )
  m_full_int <- fit_lmer(form_int_full, datC)

  fixed_now <- attr(stats::terms(lme4::nobars(stats::formula(m_full_int))), "term.labels")
  fixed_now <- fixed_now[fixed_now != ""]
  candidates_int <- setdiff(fixed_now, forced_var)

  sel_int <- backward_select_lmm(
    datC,
    response,
    forced = forced_var,
    candidates = candidates_int,
    rand_terms = rand_terms,
    independent = TRUE
  )

  list(
    full_model = m_full_int,
    selected = sel_int
  )
}

run_interaction_quadratic_bic <- function(datQ, response, forced_var, x_vars, rand_terms) {
  lin_terms  <- x_vars
  int_terms  <- combn(x_vars, 2, FUN = function(z) paste(z, collapse = ":"), simplify = TRUE)
  quad_terms <- mk_sq(x_vars)

  fixed_terms_iq <- unique(c(lin_terms, int_terms, quad_terms))

  form_iq_full <- make_formula(response, fixed_terms = fixed_terms_iq, rand_terms = rand_terms)
  m_full_iq <- fit_lmer(form_iq_full, datQ)

  fixed_now_iq <- attr(stats::terms(lme4::nobars(stats::formula(m_full_iq))), "term.labels")
  fixed_now_iq <- fixed_now_iq[fixed_now_iq != ""]
  candidates_iq <- setdiff(fixed_now_iq, forced_var)

  sel_iq <- backward_select_lmm(
    datQ,
    response,
    forced = forced_var,
    candidates = candidates_iq,
    rand_terms = rand_terms,
    independent = TRUE
  )

  list(
    full_model = m_full_iq,
    selected = sel_iq
  )
}
