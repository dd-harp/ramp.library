
#' @title Set up a model for detection
#' @description A model could have models
#' for detection by different methods, so detection
#' models must be dispatched by passing their
#' model objects.
#' @param y state vector
#' @param Xpars an `Xmod` object
#' @param mod_detect an `F_detect` model object
#' @return a [numeric] vector of length `nStrata`
#' @export
F_detect_y <- function(y, Xpars, mod_detect) {
  UseMethod("F_detect_y", mod_detect)
}

#' @title Set up a model for detection
#' @description A model could have models
#' for detection by different methods, so `F_detect`
#' must pass the specific detection model.
#' @param vars variables as a named list
#' @param Xpars an **X** model object
#' @param mod_detect an detection model object
#' @return a [numeric] vector of length `nStrata`
#' @export
F_detect_vars <- function(vars, Xpars, mod_detect) {
  UseMethod("F_detect_vars", mod_detect)
}

#' @title Set up default detect
#' @description Set up an `F_detect` model object
#' @param mod_detect a string, the S3 `class` for `F_detect`
#' @param detect_opts  options to set up the model
#' @param nStrata  the number of population strata
#' @return an `F_detect` model object
#' @export
setup_F_detect <- function(mod_detect, detect_opts, nStrata) {
  class(mod_detect) <- mod_detect
  UseMethod("setup_F_detect", mod_detect)
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements the linear model for detection
#' @inheritParams F_detect_y
#' @return a [numeric] vector of length `nStrata`
#' @export
F_detect_y.d <- function(y, Xpars, mod_detect) {
  return(mod_detect$d)
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements the linear model for detection
#' @inheritParams F_detect_vars
#' @return a [numeric] vector of length `nStrata`
#' @export
F_detect_vars.d <- function(vars, Xpars, mod_detect) {
  return(mod_detect$d)
}

#' @title Set up a basic model for detection
#' @description Set up a `F_detect.d` model object,
#' a constant fraction detected
#' @inheritParams setup_F_detect
#' @return an `F_detect` model object
#' @export
setup_F_detect.d <- function(mod_detect, detect_opts, nStrata) {
  detect_mod <- list()
  class(detect_mod) <- "d"
  detect_mod$d <- checkIt(detect_opts$d, nStrata)
  return(detect_mod)
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements a model for detection motivated
#' by parasite densities
#' @inheritParams F_detect_y
#' @return a [numeric] vector of length `nStrata`
#' @export
F_detect_y.mav <- function(y, Xpars, mod_detect) {with(mod_detect,{
  aoi = pars$Xpars[[i]]$ix$aoi_ix
  moi = pars$Xpars[[i]]$ix$moi_ix
  v0 = pars$Xpars[[i]]$ix$v0_ix
  max_detect <- d0*exp(-mod_detect$imm_rt*v0)
  dt <- max_detect*exp(-pd_aoi*aoi*moi)
  return(dt)
})}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements a model for detection motivated by
#' parasite densities and detection
#' @inheritParams F_detect_vars
#' @return a [numeric] vector of length `nStrata`
#' @export
F_detect_vars.mav <- function(vars, Xpars, mod_detect) {
  dt <- with(vars, with(mod_detect, F_detect_mav(vh, moi, aoi, d0, imm_rt, pd_aoi)))
  return(dt)
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_detect_vars.mav] for the SIP model.
#' @param vh cumulative exposure
#' @param moi mean multiplicity of infection
#' @param aoi mean age of infection
#' @param d0 maximum probability of detection
#' @param imm_rt effect modivication by
#' @param pd_aoi parasite densities by aoi
#' @return a [numeric] vector of length `nStrata`
#' @export
F_detect_mav <- function(vh, moi, aoi, d0, imm_rt, pd_aoi) {
  d_max <- d0*exp(-imm_rt*vh)
  d_moi <- d_max*exp(-pd_aoi*aoi*moi)
  return(d_moi)
}

#' @title Set up a basic model for detection
#' @description Set up a `F_detect.q` model object,
#' a constant fraction detected
#' @inheritParams setup_F_detect
#' @return an `F_detect` model object
#' @export
setup_F_detect.mav <- function(mod_detect, detect_opts, nStrata){
  detect_mod <- list()
  class(detect_mod) <- "exp"
  detect_mod$d0 <- checkIt(detect_opts$d0, nStrata)
  detect_mod$imm_rt <- checkIt(detect_opts$imm_rt, nStrata)
  detect_mod$pd_aoi <- checkIt(detect_opts$pd_aoi, nStrata)
  return(detect_mod)
}


