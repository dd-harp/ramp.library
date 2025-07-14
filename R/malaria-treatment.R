
#' @title Treat and Cure
#' @description Define a model for the
#' fraction of infections arising that
#' are treated and cured. The model should
#' align with a model for malaria incidence.
#' @param y state vector
#' @param Xpars an `Xmod` object
#' @param mod_treat an `F_treat` model object
#' @return a [numeric] vector of length `nStrata`
#' @export
F_treat <- function(y, Xpars, mod_treat) {
  UseMethod("F_treat", mod_treat)
}

#' @title Set up F_treat
#' @description Set up model objects for F_treat
#' @param mod_treat is a string that defines the S3 class
#' @param treat_opts  a set of options to override the defaults
#' @param nStrata  the number of population strata
#' @return an `Xmod` object
#' @export
setup_F_treat <- function(mod_treat, treat_opts, nStrata) {
  class(mod_treat) <- mod_treat
  UseMethod("setup_F_treat", mod_treat)
}

#' @title Basic Treat and Cure
#' @description Implements [F_treat] for the SIP model.
#' @inheritParams F_treat
#' @return a [numeric] vector of length `nStrata`
#' @export
F_treat.p <- function(y, Xpars, mod_treat) {
  return(with(Xpars$mod_treat, c(p_treat)))
}

#' @title Set up a basic F_treat
#' @description Set up model objects for a basic F_treat model
#' @inheritParams setup_F_treat
#' @return a [numeric] vector of length `nStrata`
#' @export
setup_F_treat.p <- function(mod_treat, treat_opts, nStrata) {
  treat_mod <- list()
  class(treat_mod) <- "p"
  treat_mod$p_treat <-  checkIt(treat_opts$p_treat, nStrata)
  Xpars$mod_treat <- treat_mod
  return(Xpars)
}

#' @title Basic Treat and Cure
#' @description Implements [F_treat] for the SIP model.
#' @inheritParams F_treat
#' @return a [numeric] vector of length `nStrata`
#' @export
F_treat.p3 <- function(y, Xpars, mod_treat) {
  return(with(mod_treat, c(p_sev, p_mod, p_mild)))
}

#' @title Set up a basic F_treat
#' @description Set up model objects for a basic F_treat model
#' @inheritParams setup_F_treat
#' @return a [F_treat] model object
#' @export
setup_F_treat.p3 <- function(mod_treat, treat_opts, nStrata) {
  treat_mod <- list()
  class(treat_mod) <- "p3"
  treat_mod$p_sev  <-  checkIt(treat_opts$p_sev, nStrata)
  treat_mod$p_mod  <-  checkIt(treat_opts$p_mod, nStrata)
  treat_mod$p_mild <-  checkIt(treat_opts$p_mild, nStrata)
  return(treat_mod)
}

