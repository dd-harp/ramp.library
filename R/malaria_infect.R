
#' @title Infections per Infectious Bite
#' @description Compute the probability of
#' infection, per infectious bite. Implements
#' models of naturally-acquired or vaccine-induced
#' pre-erythrocytic immunity
#' @param y state vector
#' @param Xpars an `Xmod` object
#' @return a [numeric] vector of length `nStrata`
#' @export
F_infect <- function(y, Xpars) {
  UseMethod("F_infect", Xpars$b_mod)
}

#' @title Infections per Infectious Bite
#' @description The basic model for the
#' risk of infection, per infectious bite:
#' a constant probability of
#' infection, per infectious bite
#' @inheritParams F_infect
#' @return a [numeric] vector of length `nStrata`
#' @export
F_infect.b <- function(y, Xpars) {
  return(Xpars$b_mod$b)
}

#' @title Infections per Infectious Bite
#' @description The basic model for the
#' risk of infection, per infectious bite:
#' a constant probability of
#' infection, per infectious bite
#' @inheritParams F_infect
#' @return a [numeric] vector of length `nStrata`
#' @export
F_infect.im <- function(y, Xpars){
  vh <- y[Xpars$ix$vh_ix]
  with(Xpars$b_mod,F_infect_im(vh, b0, Nb))
}

#' @title Infections per Infectious Bite
#' @description The basic model for the
#' risk of infection, per infectious bite:
#' a constant probability of
#' infection, per infectious bite
#' @param vh cumulative exposure
#' @param b0 infection, per infectious bite, no immunity
#' @param Nb effective number of infections
#' @return a [numeric] vector of length `nStrata`
#' @export
F_infect_im <- function(vh, b0, Nb) {
  return(b0*exp(-vh/Nb))
}


#' @title Setup an `F_infect` model object
#' @description Set up an `F_infect` model object
#' @param mod_b is a string that defines the S3 class
#' @param b_opts  a set of options to override the defaults
#' @param Xpars an `Xmod` object
#' @param nStrata  the number of population strata
#' @return an `Xmod` model object
#' @export
setup_F_infect <- function(mod_b, b_opts, Xpars, nStrata) {
  class(mod_b) <- mod_b
  UseMethod("setup_F_infect", mod_b)
}

#' @title Setup an `F_infect.b` model object
#' @description Set up a model object
#' for [F_infect.b]
#' @inheritParams setup_F_infect
#' @return an `F_infect.b` model object
#' @export
setup_F_infect.b <- function(mod_b, b_opts, Xpars, nStrata) {
  b_pars <- list()
  class(b_pars) <- "b"
  b_pars$b <- checkIt(b_opts$b, nStrata)
  Xpars$b_mod <- b_pars
  return(Xpars)
}

#' @title Setup an `F_infect.b` model object
#' @description Set up a model object
#' for [F_infect.b]
#' @inheritParams setup_F_infect
#' @return an `F_infect.b` model object
#' @export
setup_F_infect.im <- function(mod_b, b_opts, Xpars, nStrata) {
  b_pars <- with(b_opts, setup_F_infect_im(b0, Nb, nStrata))
  Xpars$b_mod <- b_pars
  return(Xpars)
}

#' @title Setup an `F_infect.b` model object
#' @description Set up a F_infect model
#' with immunity of the form \deqn{b_0 e^{-v_h/N_b}}
#' @param b0 probability of infection per infectious bite, no immunity
#' @param Nb the number of effective types
#' @param nStrata the number of population strata
#' @return an `F_infect.b` model object
#' @export
setup_F_infect_im <- function(b0, Nb, nStrata) {
  b_pars <- list()
  class(b_pars) <- "im"
  b_pars$b0 <- checkIt(b0, nStrata)
  b_pars$Nb <- checkIt(Nb, nStrata)
  return(b_pars)
}
