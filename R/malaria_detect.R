
#' @title Set up a model for detection
#' @description A model could have models
#' for detection by different methods, so `F_detect`
#' must pass the specific detection model.
#' @param y state vector
#' @param Xpars an `Xmod` object
#' @param mod_detect an `F_detect` model object
#' @return a [numeric] vector of length `nStrata`
#' @export
F_detect <- function(y, Xpars, mod_detect) {
  UseMethod("F_detect", mod_detect)
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [detect] for the SIP model.
#' @inheritParams F_detect
#' @return a [numeric] vector of length `nStrata`
#' @export
F_detect.d <- function(y, Xpars, mod_detect) {
  return(mod_detect$q)
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
#' @description Implements [detect] for the SIP model.
#' @inheritParams F_detect
#' @return a [numeric] vector of length `nStrata`
#' @export
F_detect.exp <- function(y, Xpars, mod_detect) {
  aoi = pars$Xpar[[i]]$ix$aoi_ix
  moi = pars$Xpar[[i]]$ix$moi_ix
  v0 = pars$Xpar[[i]]$ix$v0_ix
  max_detect <- d0*exp(-mod_detect$imm_rt*v0)
  dt <- max_detect*exp(-pd_aoi*aoi*moi)
  return(dt)
}

#' @title Set up a basic model for detection
#' @description Set up a `F_detect.q` model object,
#' a constant fraction detected
#' @inheritParams setup_F_detect
#' @return an `F_detect` model object
#' @export
setup_F_detect.q <- function(mod_detect, detect_opts, nStrata) {
  detect_mod <- list()
  class(detect_mod) <- "exp"
  detect_mod$d0 <- checkIt(detect_opts$d0, nStrata)
  detect_mod$imm_rt <- checkIt(detect_opts$imm_rt, nStrata)
  detect_mod$pd_aoi <- checkIt(detect_opts$pd_aoi, nStrata)
  return(detect_mod)
}

#' @title Detection
#' @description The basic model for the
#' risk of transmission, per bloodmeal:
#' a constant probability of
#' transmission, per bloodmeal
#' @inheritParams F_detect
#' @return a [numeric] vector of length `nStrata`
#' @export
F_detect.dim <- function(y, Xpars, mod_detect) {with(mod_detect,{
  vh <- y[Xpars$ix$vh_ix]
  moi <- y[Xpars$ix$moi_ix]
  aoi <- y[Xpars$ix$aoi_ix]
  im_es = exp(-vh/Nd)
  d <- d0*(1-exp(-moi*exp(-im_is*aoi/pda)))
  return(c)
})}

#' @title Setup an `F_detect.dim` model object
#' @description Set up a model object
#' for [F_detect.dim]
#' @inheritParams F_detect
#' @return an `F_detect` model object
#' @export
setup_F_detect.dim <- function(mod_d, d_opts, Xpars, nStrata) {
  d_pars <- with(d_opts, setup_F_detect_dim(d0, pda, nStrata))
  Xpars$d_mod <- d_pars
  return(Xpars)
}

#' @title Setup an `F_detect.dim` model object
#' @description Set up a F_detect model
#' with immunity of the form \deqn{d_0(1- e^{-\alpha}}
#' @param c0 probability of detection per detectious cite, no immunity
#' @param pda a model for detection by aoi, motivated by a model for parasite densities and detection
#' @param nStrata the number of population strata
#' @return an `F_detect.c` model object
#' @export
setup_F_detect_dim <- function(c0, pda, nStrata) {
  d_pars <- list()
  class(d_pars) <- "dim"
  d_pars$d0 <- checkIt(d0, nStrata)
  d_pars$pda <- checkIt(pda, nStrata)
  d_pars$Nd <- checkIt(Nd, nStrata)
  return(d_pars)
}
