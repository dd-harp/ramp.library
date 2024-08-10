
#' @title Dawn, day, dusk, night model for the blood feeding rate
#' @description Implements [F_f] for a dddn model
#' @inheritParams ramp.xds::F_f
#' @return a [numeric] vector of length `nPatches`
#' @export
F_f.dddn <- function(t, vars, f_par){
  with(vars, f_par, t(matrix(f0, 4, nPatches)))
}

#' @title Dawn, day, dusk, night model for the human fraction
#' @description Implements [F_q] for a dddn model
#' @inheritParams ramp.xds::F_q
#' @return a [numeric] vector of length `nPatches`
#' @export
F_q.dddn <- function(t, vars, q_par){
  with(vars, q_par, t(matrix(q0, 4, nPatches)))
}

#' @title Dawn, day, dusk, night model for the human fraction
#' @description Implements [F_p] for a dddn model
#' @inheritParams ramp.xds::F_p
#' @return a [numeric] vector of length `nPatches`
#' @export
F_p.dddn <- function(t, vars, p_par){
  with(vars, p_par, t(matrix(p0, 4, nPatches)))
}

#' @title Dawn, day, dusk, night model for the human fraction
#' @description Implements [F_sigma] for a dddn model
#' @inheritParams ramp.xds::F_sigma
#' @return a [numeric] vector of length `nPatches`
#' @export
F_sigma.dddn <- function(t, vars, sigma_par){
  with(vars, sigma_par, t(matrix(sigma0, 4, nPatches)))
}

#' @title Dawn, day, dusk, night model for the human fraction
#' @description Implements [F_nu] for a dddn model
#' @inheritParams ramp.xds::F_nu
#' @return a [numeric] vector of length `nPatches`
#' @export
F_nu.dddn <- function(t, vars, nu_par){
  with(vars, nu_par, t(matrix(nu0, 4, nPatches)))
}
