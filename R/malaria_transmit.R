
#' @title Infectiousness
#' @description This method computes infectiousness,
#' the probability a mosquito would become infected
#' after blood feeding on a human, by strata
#' @param y state vector
#' @param Xpars an `Xmod` object
#' @return a [numeric] vector of length `nStrata`
#' @export
F_transmit <- function(y, Xpars) {
  UseMethod("F_transmit", Xpars$c_mod)
}

#' @title Infectiousness, basic
#' @description Return a constant value, \eqn{c},
#' for the probability a mosquito would become infected
#' after blood feeding on a human, by strata
#' @inheritParams F_transmit
#' @return a [numeric] vector of length `nStrata`
#' @export
F_transmit.c <- function(y, Xpars) {
  return(Xpars$c_mod$c)
}

#' @title Setup Infectiousness
#' @description Set up an `F_transmit` model object
#' @param mod_c is a string that defines the S3 class
#' @param c_opts  a set of options to override the defaults
#' @param Xpars an `Xmod` object
#' @param nStrata  the number of population strata
#' @return an `Xmod` model object
#' @export
setup_F_transmit <- function(mod_c, c_opts, Xpars, nStrata) {
  class(mod_c) <- mod_c
  UseMethod("setup_F_transmit", mod_c)
}

#' @title Set up default F_transmit
#' @description Implements [F_transmit] for the SIP model.
#' @inheritParams setup_F_transmit
#' @return an `Xmod` model object
#' @export
setup_F_transmit.c <- function(mod_c, c_opts, Xpars, nStrata) {
  c_pars <- list()
  class(c_pars) <- "c"
  c_pars$c <- checkIt(c_opts$c, nStrata)
  Xpars$c_mod <- c_pars
  return(Xpars)
}

#' @title A model for detection
#' @description The basic model for the
#' risk of transmission, per bloodmeal:
#' a constant probability of
#' transmission, per bloodmeal
#' @inheritParams F_transmit
#' @return a [numeric] vector of length `nStrata`
#' @export
F_transmit.dim <- function(y, Xpars) {with(Xpars$c_mod,{
  vh <- y[Xpars$ix$vh_ix]
  moi <- y[Xpars$ix$moi_ix]
  aoi <- y[Xpars$ix$aoi_ix]
  im_es = exp(-vh/Nc)
  c <- c0*(1-exp(-moi*exp(-im_es*aoi/pda)))
  return(c)
})}

#' @title Setup an `F_transmit.dim` model object
#' @description Set up a model object
#' for [F_transmit.c]
#' @inheritParams setup_F_transmit
#' @return an `F_transmit.c` model object
#' @export
setup_F_transmit.dim <- function(mod_c, c_opts, Xpars, nStrata) {
  c_pars <- with(c_opts, setup_F_transmit_dim(c0, pda, Nc, nStrata))
  Xpars$c_mod <- c_pars
  return(Xpars)
}

#' @title Setup an `F_transmit.dim` model object
#' @description Set up a F_transmit model
#' with immunity of the form \deqn{c_0 e^{-v_h/N_c}}
#' @param c0 probability of transmition per transmitious cite, no immunity
#' @param pda a model for detection by aoi, motivated by a model for parasite densities and detection
#' @param Nc a model for protection by exposure
#' @param nStrata the number of population strata
#' @return an `F_transmit.c` model object
#' @export
setup_F_transmit_dim <- function(c0, pda, Nc, nStrata) {
  c_pars <- list()
  class(c_pars) <- "dim"
  c_pars$c0 <- checkIt(c0, nStrata)
  c_pars$pda <- checkIt(pda, nStrata)
  c_pars$Nc <- checkIt(Nc, nStrata)
  return(c_pars)
}
