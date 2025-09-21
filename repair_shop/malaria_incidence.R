
#' @title Malaria Incidence from Exposure
#' @description Compute the probability of
#' getting malaria as a fraction of new
#' infections.
#' @param y state vector
#' @param Xpars an `Xmod` object
#' @param mod_ince_h  an `F_incidence_h` model object
#' @return a [numeric] vector of length `nStrata`
#' @export
F_incidence_h <- function(y, Xpars, mod_ince_h) {
  UseMethod("F_incidence_h", mod_ince_h)
}

#' @title Malaria Incidence from Infection
#' @description Compute the rate that malaria
#' cases arise from infected states
#' @param y state vector
#' @param Xpars an `Xmod` object
#' @param mod_ince_x  an `F_incidence_x` model object
#' @return a [numeric] vector of length `nStrata`
#' @export
F_incidence_x <- function(y, Xpars, mod_ince_x) {
  UseMethod("F_incidence_x", mod_ince_x)
}

#' @title Setup Incidence from Exposure
#' @description Set up an `F_incidence_h` model object.
#' @param mod_ince_h is a string that defines the S3 class
#' @param malaria_opts a set of options to override the defaults
#' @param nStrata  the number of population strata
#' @description an `F_incidence_h` model object
#' @export
setup_incidence_h <- function(mod_ince_h, malaria_opts, nStrata) {
  class(mod_ince_h) <- mod_ince_h
  UseMethod("setup_incidence_h", mod_ince_h)
}



#' @title Setup Malaria from Infection
#' @description Set up a `F_incidence_x` model object.
#' @param mod_ince_x is a string that defines the S3 class
#' @param malaria_opts a set of options to override the defaults
#' @param nStrata  the number of population strata
#' @description an `F_incidence_x` model object
#' @export
setup_incidence_x <- function(mod_ince_x, malaria_opts, nStrata) {
  class(mod_ince_x) <- mod_ince_x
  UseMethod("setup_incidence_x", mod_ince_x)
}

#' @title Malaria from Exposure
#' @description Computes the probability of
#' getting malaria from a new infection
#' @inheritParams F_incidence_h
#' @return a [numeric] vector of length `nStrata`
#' @export
F_incidence_h.p <- function(y, Xpars, mod_ince_h) {
  return(mod_ince_h$p)
}

#' @title Setup Malaria from Exposure, Basic
#' @description Set up an `F_incidence_h.p` model object
#' @inheritParams setup_incidence_h
#' @description an `F_incidence_h` model object
#' @export
setup_incidence_h.p <- function(mod_ince_h, malaria_opts, nStrata) {
  pars <- list()
  class(pars) <- "p"
  pars$p <- checkIt(malaria_opts$p, nStrata)
  return(pars)
}

#' @title Basic model for
#' @description Implements a basic model for
#' incidence from infection
#' @inheritParams F_incidence_x
#' @return a [numeric] vector of length `nStrata`
#' @export
F_incidence_x.rate <- function(y, Xpars, mod_ince_x){
  return(mod_ince_x$rate)
}

#' @title Setup Malaria from Infection, Basic
#' @description Set up a basic model for malaria
#' incidence from infection, a constant
#' rate of malaria per infection.
#' @inheritParams setup_incidence_x
#' @description an `F_incidence_x` model object
#' @export
setup_incidence_x.rate <- function(mod_ince_x, malaria_opts, nStrata) {
  pars <- list()
  class(pars) <- "rate"
  pars$rate <- checkIt(malaria_opts$rate, nStrata)
  return(pars)
}

#' @title Malaria Incidents, per infection
#' @description A model for incidence
#' assuming a *bad types* model
#' @inheritParams F_incidence_h
#' @return a [numeric] vector of length `nStrata`
#' @export
F_incidence_h.protect <- function(y, Xpars, mod_ince_h) {with(mod_ince_h,{
  vh <- y[Xpars$ix$vh_ix]
  p <- p0*exp(-vh/Np)
  return(p)
})}

#' @title Setup an `F_incidence_h.protect` model object
#' @description Set up a model object
#' for [F_incidence_h.protect]
#' @inheritParams setup_incidence_h
#' @return an `F_incidence_h.b` model object
#' @export
setup_incidence_h.protect <- function(mod_ince_h, malaria_opts,  nStrata) {
  ince_pars <- with(malaria_opts, setup_incidence_h_protect(p0, Np, nStrata))
  return(ince_pars)
}

#' @title Setup an `F_incidence_h.protect` model object
#' @description Set up a F_incidence_h model
#' with immunity of the form \deqn{p_0 e^{-v_h/N_p}}
#' @param p0 probability of malaria, per infection
#' @param Np the number of effective types
#' @param nStrata the number of population strata
#' @return an `F_incidence_h.protect` model object
#' @export
setup_incidence_h_protect <- function(p0, Np, nStrata) {
  p_pars <- list()
  class(p_pars) <- "protect"
  p_pars$p0 <- checkIt(p0, nStrata)
  p_pars$Np <- checkIt(Np, nStrata)
  return(p_pars)
}

#' @title Malaria Incidents, per infection
#' @description A model for incidence
#' assuming a *bad types* model
#' @inheritParams F_incidence_x
#' @return a [numeric] vector of length `nStrata`
#' @export
F_incidence_x.protect <- function(y, Xpars, mod_ince_x) {with(mod_ince_x,{
  vh <- y[Xpars$ix$vh_ix]
  rate <- rate0*exp(-vh/Np)
  return(rate)
})}

#' @title Setup an `F_incidence_x.protect` model object
#' @description Set up a model object
#' for [F_incidence_x.protect]
#' @inheritParams setup_incidence_x
#' @return an `F_incidence_h.b` model object
#' @export
setup_incidence_x.protect <- function(mod_ince_x, malaria_opts,  nStrata) {
  ince_pars <- with(malaria_opts, setup_incidence_x_protect(rate0, Np, nStrata))
  return(ince_pars)
}

#' @title Setup an `F_incidence_h.protect` model object
#' @description Set up a F_incidence_h model
#' with immunity of the form \deqn{p_0 e^{-v_h/N_p}}
#' @param rate0 probability of malaria, per infection
#' @param Np the number of effective types
#' @param nStrata the number of population strata
#' @return an `F_incidence_x.protect` model object
#' @export
setup_incidence_x_protect <- function(rate0, Np, nStrata) {
  p_pars <- list()
  class(p_pars) <- "protect"
  p_pars$rate0 <- checkIt(rate0, nStrata)
  p_pars$Np <- checkIt(Np, nStrata)
  return(p_pars)
}
