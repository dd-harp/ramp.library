
#' @title Natural Parasite Clearance
#' @description Compute natural parasite clearance
#' in a population
#' @param y state vector
#' @param pars an **`xds`** object
#' @return a [numeric] vector of length `nStrata`
#' @export
F_clear <- function(y, Xpars) {
  UseMethod("F_clear", Xpars$r_mod)
}

#' @title Natural Parasite Clearance
#' @description Implements [F_clear] for a constant
#' rate of clearance for simple infections
#' @inheritParams F_clear
#' @return a [numeric] vector of length `nStrata`
#' @export
F_clear.r <- function(y, Xpars) {
  return(Xpars$r_mod$r)
}

#' @title Setup an `F_clear` model object
#' @description Setup an `F_clear` model object
#' @param mod_clear is a string that defines the S3 class
#' @param r_opts  a set of options to override the defaults
#' @param Xpars an `Xmod` object
#' @param nStrata  the number of population strata
#' @return an `Xmod` object
#' @export
setup_F_clear <- function(mod_clear, r_opts, Xpars, nStrata) {
  class(mod_clear) <- mod_clear
  UseMethod("setup_F_clear", mod_clear)
}

#' @title Set up default F_clear
#' @description Implements [F_clear] for the SIP model.
#' @inheritParams setup_F_clear
#' @return an `Xmod` object
#' @export
setup_F_clear.r <- function(mod_clear, r_opts, Xpars, nStrata) {
  r_pars <- list()
  class(r_pars) <- "r"
  r_pars$r <- checkIt(r_opts$r, nStrata)
  Xpars$r_mod <- r_pars
  return(Xpars)
}
