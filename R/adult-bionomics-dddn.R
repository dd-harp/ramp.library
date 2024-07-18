
#' @title Dawn, day, dusk, night model for the blood feeding rate
#' @description Implements [F_f] for a dddn model
#' @inheritParams ramp.xds::F_f
#' @return a [numeric] vector of length `nPatches`
#' @export
F_f.dddn <- function(t, MYZpar){
  with(MYZpar, t(matrix(f0, 4, nPatches)))
}

#' @title Dawn, day, dusk, night model for the human fraction
#' @description Implements [F_q] for a dddn model
#' @inheritParams ramp.xds::F_q
#' @return a [numeric] vector of length `nPatches`
#' @export
F_q.dddn <- function(t, MYZpar){
  with(MYZpar, t(matrix(q0, 4, nPatches)))
}

#' @title Dawn, day, dusk, night model for the human fraction
#' @description Implements [F_p] for a dddn model
#' @inheritParams ramp.xds::F_p
#' @return a [numeric] vector of length `nPatches`
#' @export
F_p.dddn <- function(t, MYZpar){
  with(MYZpar, t(matrix(p0, 4, nPatches)))
}

#' @title Dawn, day, dusk, night model for the human fraction
#' @description Implements [F_sigma] for a dddn model
#' @inheritParams ramp.xds::F_sigma
#' @return a [numeric] vector of length `nPatches`
#' @export
F_sigma.dddn <- function(t, MYZpar){
  with(MYZpar, t(matrix(sigma0, 4, nPatches)))
}

#' @title Dawn, day, dusk, night model for the human fraction
#' @description Implements [F_nu] for a dddn model
#' @inheritParams ramp.xds::F_nu
#' @return a [numeric] vector of length `nPatches`
#' @export
F_nu.dddn <- function(t, MYZpar){
  with(MYZpar, t(matrix(nu0, 4, nPatches)))
}
