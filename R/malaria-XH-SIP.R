# specialized methods for the human SIP model



#' @title The `SIP` module for the **XH** component
#' @description
#' Implements the **XH** component using a Susceptible-Infected-chemoProtected
#' (`SIP`) compartmental model of malaria infection dynamics.
#' Treatment clears infection and moves individuals into \eqn{P}, where they
#' are protected from infection; they re-enter susceptibility when chemoprotection wanes.
#'
#' This implementation includes:
#' + parameters to model treating and curing infections;
#' + functions to model mass treatment;
#' + functions to model human demographic change;
#' + support for progressive stratification.
#'
#' @section State Variables:
#' \describe{
#'   \item{`H`}{total human (or host) population density}
#'   \item{`I`}{density of infectious humans}
#'   \item{`P`}{density of humans under chemoprotection}
#'   \item{`S=H-I-P`}{density of susceptible humans or hosts}
#' }
#'
#' @inheritSection ramp.xds::SIS Demography
#'
#' @section Dynamical Term:
#' \describe{
#'   \item{\eqn{h} or `foi`}{The force of infection is computed and stored by [Exposure] as `FoI`}
#' }
#' In the `SIP` module, a constant fraction of infectious bites causes an
#' infection: the parameter \eqn{b}.
#'
#' @section Parameters:
#' The parameters returned by [get_XH_pars.SIP] and modifiable by
#' [change_XH_pars.SIP] describe infection, treatment, and chemoprotection.
#' \describe{
#'   \item{`q`}{test positivity (default `0.8`)}
#'   \item{`b`}{transmission probability from mosquito to human}
#'   \item{`c`}{transmission probability from human to mosquito}
#'   \item{`r`}{clearance rate for infections}
#'   \item{`rho` or \eqn{\rho}}{fraction of incident infections treated and cured, entering \eqn{P}}
#'   \item{`xi` or \eqn{\xi}}{background treatment rate for susceptible and infectious individuals}
#'   \item{`sigma` or \eqn{\sigma}}{additional treatment rate for infectious individuals}
#'   \item{`eta`}{rate of chemoprotection waning}
#' }
#'
#' @section Infection Dynamics:
#' Infections clear naturally at rate \eqn{r}, returning individuals to
#' susceptibility. Treatment instead moves individuals into chemoprotection.
#' A fraction \eqn{\rho} of incident infections \eqn{hS} are treated immediately.
#' Define \eqn{\xi_t = \xi + mda(t)} and \eqn{\sigma_t = \sigma + msat(t)}.
#' \deqn{
#' \begin{array}{rl}
#' dI/dt &= (1-\rho)hS - (r + \xi_t + \sigma_t)I - \mathbf{D}(I) \\
#' dP/dt &= \rho hS + \xi_t(S+I) + \sigma_t I - \eta P - \mathbf{D}(P) \\
#' \end{array}}
#' where \eqn{h} is the force of infection and \eqn{S = H - I - P}.
#' The mortality function uses the same model for \eqn{I} and \eqn{P} as for
#' \eqn{H}. Susceptible density is redundant and is computed after solving.
#'
#' @section Prevalence:
#' [F_pr_obs.SIP] returns observed prevalence \eqn{x_q=qI/H}.
#'
#' True prevalence is \eqn{x = I/H}. The light microscopy, RDT, and PCR
#' prevalence functions currently also return \eqn{I/H}, without diagnostic
#' sensitivity parameters.
#'
#' @section Infectiousness:
#' The probability that a mosquito becomes infected after feeding on an
#' infectious host is \eqn{c}. Infectious density is \eqn{cI}, and net
#' infectiousness is \eqn{cI/H}.
#'
#' @inheritSection ramp.xds::SIS Mass Treatment
#'
#' @note In `SIP`, mass drug administration treats susceptible and infectious
#' individuals, while mass screen and treat acts on infectious individuals.
#' Treatment moves individuals into chemoprotection \eqn{P}, not susceptibility.
#' Background treatment rates
#' `xi` and `sigma` remain active when the mass-treatment functions are zero.
#'
#' @seealso
#' + Basic Setup: [make_XH_obj_SIP], [setup_XH_obj.SIP]
#' + Dynamics: [dXHdt.SIP]
#' + Parameters: [get_XH_pars.SIP], [change_XH_pars.SIP]
#' + Variables: [make_XH_inits_SIP], [get_XH_inits], [change_XH_inits.SIP], [setup_XH_inits.SIP]
#' + [SIS]
#'
#' @name SIP
#' @rdname SIP
NULL

#' @title Compute the derivatives for `SIP` (**XH**)
#' @description The derivatives for the `SIP` module are:
#' \deqn{\begin{array}{rl}
#' dH/dt &= \mathbf{B}(t,y)-\mathbf{D}(t,H) \\
#' dI/dt &= (1-\rho)hS-(r+\xi_t+\sigma_t)I-\mathbf{D}(t,I) \\
#' dP/dt &= \rho hS+\xi_t(S+I)+\sigma_t I-\eta P-\mathbf{D}(t,P)
#' \end{array}}
#' Here \eqn{h} is `foi`, \eqn{S=H-I-P},
#' \eqn{\xi_t=\xi+mda(t)}, and \eqn{\sigma_t=\sigma+msat(t)}.
#' The term \eqn{\mathbf{B}(t,y)} is returned by [ramp.xds::Births], and
#' \eqn{\mathbf{D}(t,X)} is returned by [ramp.xds::Deaths] for population \eqn{X}.
#' @inheritSection SIP Demography
#' @inheritSection SIP Infection Dynamics
#' @inheritSection SIP Mass Treatment
#'
#' @seealso [SIP], [make_XH_obj_SIP]
#' @inheritParams ramp.xds::dXHdt
#' @return a [numeric] vector
#' @keywords internal
#' @export
dXHdt.SIP <- function(t, y, xds_obj, i){

  foi <- xds_obj$terms$FoI[[i]]

  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {

      xi_t = xi + mda(t)
      sigma_t = sigma + msat(t)

      dH <- Births(t, y, xds_obj, i) - Deaths(t, H, xds_obj, i)
      dI <- (1-rho)*foi*S - (r+xi_t+sigma_t)*I - Deaths(t, I, xds_obj, i)
      dP <- rho*foi*S + xi_t*(S+I) + sigma_t*I - eta*P - Deaths(t, P, xds_obj, i)

      return(c(dH, dI, dP))
    })
  })
}

#' @title Update states for `SIP` (**XH**)
#' @description Implements [Update_XHt] for the `SIP` model.
#' @inheritParams ramp.xds::Update_XHt
#' @return a [numeric] vector
#' @keywords internal
#' @seealso [SIP]
#' @export
Update_XHt.SIP <- function(t, y, xds_obj, i){

  attack <- xds_obj$AR[[i]]

  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {

      St <- (1-attack)*(S+r*I) + eta*P - xi*S
      It <- (1-r)*I + attack*(1-rho)*(S+r*I) - xi*I
      Pt <- xi*(S+I) + attack*rho*(S+r*I) + (1-eta)*P

      St <- dHdt(t, St, xds_obj$Hpar[[i]]) + Births(t, y, xds_obj, i)
      It <- dHdt(t, It, xds_obj$Hpar[[i]])
      Pt <- dHdt(t, Pt, xds_obj$Hpar[[i]])

      return(c(St, It, Pt))
    })
  })
}


#' @title Get the state variables for `SIP` (**XH**)
#' @description Extract the state variables for the `SIP` module
#' from `y` and return them by name in a list for use during model solving.
#' Also compute susceptible density `S = H - I - P`.
#' @inheritParams ramp.xds::get_XH_vars
#' @return a named [list] of state variables and module-specific derived quantities
#' @keywords internal
#' @seealso [SIP]
#' @export
#' @note Variable indices are stored at `xds_obj$XH_obj[[i]]$ix`.
get_XH_vars.SIP <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
    H = y[H_ix]
    I = y[I_ix]
    P = y[P_ix]
    S = H-I-P
    return(list(S=S,I=I,P=P,H=H))})
}

#' @title Size of effective infectious human population
#' @description Implements [F_I] for the `SIP` model.
#' @inheritParams ramp.xds::F_I
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIP]
#' @export
F_I.SIP <- function(t, y, xds_obj, i) {
  I = y[xds_obj$XH_obj[[i]]$ix$I_ix]
  X = with(xds_obj$XH_obj[[i]], c*I)
  return(X)
}

#' @title Compute population density
#' @description Implements [F_H] for the `SIP` model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIP]
#' @export
F_H.SIP <- function(t, y, xds_obj, i){
  with(get_XH_vars(y, xds_obj, i), return(H))
}


#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_infectivity] for the `SIP` model.
#' @inheritParams ramp.xds::F_infectivity
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIP]
#' @export
F_infectivity.SIP <- function(y, xds_obj,i) {
  with(xds_obj$XH_obj[[i]], b)
}


#' @title Make a module object for `SIP` (**XH**)
#' @description Create an `SIP` module object and initialize its parameters.
#' Inspect parameters with [get_XH_pars.SIP] and change supported values
#' with [change_XH_pars.SIP]. Values in `options` override parameter arguments.
#' @seealso [get_XH_pars.SIP], [change_XH_pars.SIP]
#' @inheritSection SIP Parameters
#' @param nStrata the number of population strata
#' @param options a [list] that could override defaults
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @param rho fraction of incident infections treated and cured, entering P
#' @param eta chemoprotection waning rate
#' @param xi background treatment rate for susceptible and infectious individuals
#' @param sigma additional treatment rate for infectious individuals
#' @importFrom ramp.xds F_zero checkIt
#' @return an **XH** model object
#' @keywords internal
#' @seealso [SIP]
#' @param q test positivity; a numeric scalar or input of length `nStrata`
#' @export
make_XH_obj_SIP = function(nStrata, options=list(),
                           b=0.55, r=1/180, c=0.15,
                           rho =.1, eta=1/25,
                           xi =1/365, sigma = 1/365, q=0.8){
  with(options,{
    XH_obj = list()
    class(XH_obj) <- c("SIP")

    XH_obj$b = checkIt(b, nStrata)
    XH_obj$c = checkIt(c, nStrata)
    XH_obj$r = checkIt(r, nStrata)
    XH_obj$rho = checkIt(rho, nStrata)
    XH_obj$eta = checkIt(eta, nStrata)
    XH_obj$xi = checkIt(xi, nStrata)
    XH_obj$sigma = checkIt(sigma, nStrata)


    XH_obj$q = checkIt(q, nStrata)

    return(XH_obj)
  })}

#' @title Set up the module object for `SIP` (**XH**)
#' @description
#'
#' Set up an `SIP` module object for the **XH** component and the \eqn{i^{th}} host species.
#'
#' + Set the module name to `SIP`.
#' + Create the module object with [make_XH_obj_SIP], passing `nStrata[i]` and `options`.
#' + Set initial values using [ramp.xds::setup_XH_inits], passing `options` to override defaults.
#' + Initialize the time-spent matrix using [ramp.xds::setup_timespent.setup].
#' + Initialize mass-treatment functions using [ramp.xds::setup_mass_treatment].
#' + Initialize birth and death rates using [ramp.xds::setup_demog_obj].
#' + Initialize blood search weights using [ramp.xds::setup_blood_search_weights.setup].
#' + Initialize time away using [ramp.xds::setup_time_away.setup].
#' + Initialize travel EIR using [ramp.xds::setup_travel_eir.no_travel].
#' + Return the updated **`xds`** model object.
#'
#' @seealso [SIP], [make_XH_obj_SIP]
#'
#' @note The module object is stored at `xds_obj$XH_obj[[i]]`
#' @inheritSection SIP State Variables
#' @inheritSection SIP Parameters
#'
#' @inheritParams ramp.xds::setup_XH_obj
#' @return an **`xds`** model object
#' @keywords internal
#' @export
setup_XH_obj.SIP = function(XHname, residence, xds_obj, options=list(), i=1){
  xds_obj$XHname = "SIP"
  xds_obj$XH_obj[[i]] = make_XH_obj_SIP(xds_obj$nStrata[i], options)
  xds_obj <- setup_XH_inits(xds_obj, options=options, i=i)
  xds_obj <- setup_timespent("setup", xds_obj, list(residence=residence), i)
  xds_obj <- setup_mass_treatment(xds_obj, i=i)
  xds_obj <- setup_demog_obj(xds_obj, i)
  xds_obj <- setup_blood_search_weights("setup", xds_obj, i=i)
  xds_obj <- setup_time_away("setup", xds_obj, i=i)
  xds_obj <- setup_travel_eir("no_travel", xds_obj, i=i)
  return(xds_obj)
}


#' @title Get parameters for `SIP` (**XH**)
#' @seealso [make_XH_obj_SIP], [change_XH_pars.SIP]
#' @description Return the parameters for the `SIP` module as a named list.
#' @note Parameter values for the \eqn{i^{th}} host are
#' stored as `xds_obj$XH_obj[[i]]`.
#' @inheritParams ramp.xds::get_XH_pars
#' @inheritSection SIP Parameters
#' @return a named [list] of parameter values
#' @seealso [SIP], [make_XH_obj_SIP]
#' @keywords internal
#' @export
get_XH_pars.SIP <- function(xds_obj, i=1) {
  with(xds_obj$XH_obj[[i]],list(q=q, b=b, c=c, r=r, rho=rho, eta=eta, xi=xi, sigma=sigma))
}



#' @title Change parameters for `SIP` (**XH**)
#' @seealso [make_XH_obj_SIP], [get_XH_pars.SIP]
#' @description Change parameter values named in `options` for the `SIP`
#' module, leaving other values unchanged.
#' @inheritParams ramp.xds::change_XH_pars
#' @inheritSection SIP Parameters
#' @return an **`xds`** model object
#' @keywords internal
#' @seealso [SIP]
#' @export
change_XH_pars.SIP <- function(xds_obj, options=list(), i=1) {
  nHabitats <- xds_obj$nHabitats
  with(xds_obj$XH_obj[[i]], with(options,{
    xds_obj$XH_obj[[i]]$b <- b
    xds_obj$XH_obj[[i]]$c <- c
    xds_obj$XH_obj[[i]]$r <- r
    xds_obj$XH_obj[[i]]$rho <- rho
    xds_obj$XH_obj[[i]]$eta <- eta
    xds_obj$XH_obj[[i]]$xi <- xi
    xds_obj$XH_obj[[i]]$sigma <- sigma
    xds_obj$XH_obj[[i]]$mda <- mda
    xds_obj$XH_obj[[i]]$msat <- msat
    xds_obj$XH_obj[[i]]$q <- q
    return(xds_obj)
  }))}

#' @title Make initial values for `SIP` (**XH**)
#' @description Create a named list of initial values for the `SIP` module.
#' Values in `options` override the defaults.
#' @details Initial values must be numeric scalars or have length `nStrata`.
#' Scalars are expanded to that length.
#' @inheritSection SIP State Variables
#' @param nStrata the number of population strata
#' @param options a [list] that could override defaults
#' @param H the initial human population density
#' @param I the initial number of infectious humans
#' @param P the initial number of humans under chemoprotection
#' @return a named [list] of initial values
#' @keywords internal
#' @seealso [SIP]
#' @export
#' @seealso [setup_XH_inits.SIP], [change_XH_inits.SIP]
make_XH_inits_SIP = function(nStrata, options = list(),
                             H=1000, I=1, P=0){
with(options,{
  H = checkIt(H, nStrata)
  I = checkIt(I, nStrata)
  P = checkIt(P, nStrata)
  return(list(H=H, I=I, P=P))
})}

#' @title Set up initial values for `SIP` (**XH**)
#' @description Set up initial values for the `SIP` module using
#' [make_XH_inits_SIP]. Values in `options` override the defaults.
#' @inheritParams ramp.xds::setup_XH_inits
#' @inheritSection SIP State Variables
#' @return an **`xds`** model object
#' @keywords internal
#' @seealso [SIP]
#' @export
#' @seealso [change_XH_inits.SIP]
#' @note Initial values are stored at `xds_obj$XH_obj[[i]]$inits`.
#' @seealso [ramp.xds::get_XH_inits]
setup_XH_inits.SIP = function(xds_obj, options=list(), i=1){
  xds_obj$XH_obj[[i]]$inits = make_XH_inits_SIP(xds_obj$nStrata[i], options)
  return(xds_obj)
}


#' @title Change initial values for `SIP` (**XH**)
#' @description Change initial values named in `options` for the `SIP`
#' module, leaving other values unchanged.
#' @inheritParams ramp.xds::change_XH_inits
#' @inheritSection SIP State Variables
#' @note This method changes `I` and `P`; it leaves `H` unchanged.
#' @return an **`xds`** model object
#' @keywords internal
#' @seealso [SIP]
#' @export
#' @seealso [make_XH_inits_SIP], [setup_XH_inits.SIP]
#' @note Stored initial values for this module are at
#' `xds_obj$XH_obj[[i]]$inits`.
change_XH_inits.SIP <- function(xds_obj, options=list(), i=1) {
  with(get_XH_inits(xds_obj, i), with(options,{
    xds_obj$XH_obj[[i]]$inits$I = I
    xds_obj$XH_obj[[i]]$inits$P = P
    return(xds_obj)
  }))}

#' @title Set up variable indices for `SIP` (**XH**)
#' @description Assign consecutive variable indices for the `SIP` module,
#' starting after `xds_obj$max_ix`, and update the maximum index.
#' Allocate blocks in the order `H_ix`, `I_ix`, `P_ix`.
#' Each block has length `nStrata[i]`.
#' @inheritParams ramp.xds::setup_XH_ix
#' @return an **`xds`** model object
#' @importFrom utils tail
#' @keywords internal
#' @seealso [SIP]
#' @export
#' @note Indices are stored at `xds_obj$XH_obj[[i]]$ix`.
#' The maximum index is stored at `xds_obj$max_ix`.
setup_XH_ix.SIP <- function(xds_obj, i) {with(xds_obj,{

  H_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(H_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(I_ix, 1)

  P_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(P_ix, 1)

  xds_obj$max_ix = max_ix
  xds_obj$XH_obj[[i]]$ix = list(H_ix=H_ix, I_ix=I_ix, P_ix=P_ix)
  return(xds_obj)
})}

#' @title Parse outputs for `SIP` (**XH**)
#'
#' @description Extract the model state outputs for the `SIP` module
#' from `outputs` and return them by name in a list.
#' Also compute `S = H-I-P`.
#' @note Variable indices are stored at `xds_obj$XH_obj[[i]]$ix`.
#'
#' @inheritParams ramp.xds::parse_XH_orbits
#'
#' @return a [list]
#' @keywords internal
#' @seealso [SIP]
#' @export
parse_XH_orbits.SIP <- function(outputs, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
    H = outputs[,H_ix]
    I = outputs[,I_ix]
    P = outputs[,P_ix]
    S = H-I-P
    vars <- list(S=S, I=I, P=P, H=H)
    return(vars)
})}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_prevalence] for the `SIP` model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIP]
#' @export
F_prevalence.SIP <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}


#' @title Compute observed prevalence
#' @description Return observed prevalence for the `SIP` module as
#' `q` times the true prevalence returned by [ramp.xds::F_prevalence].
#' @inheritParams ramp.xds::F_pr_obs
#' @return a [numeric] vector of length `nStrata`
#' @seealso [SIP]
#' @keywords internal
#' @export
F_pr_obs.SIP <- function(vars, XH_obj) {
  XH_obj$q * ramp.xds::F_prevalence(vars, XH_obj)
}

#' @title Compute the net infectiousness
#' @description Implements [F_ni] for the `SIP` model.
#' @inheritParams ramp.xds::F_ni
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIP]
#' @export
F_ni.SIP <- function(vars, XH_obj) {
  ni = with(vars, XH_obj$c*I/H)
  return(ni)
}


#' @title Compute the HTC for the `SIP` model
#' @description Implements [get_HTC] for the `SIP` model with demography.
#' @inheritParams ramp.xds::get_HTC
#' @return a [numeric] vector
#' @keywords internal
#' @seealso [SIP]
#' @export
get_HTC.SIP <- function(xds_obj, i) {
  with(xds_obj$XH_obj[[i]],
       return((1-rho)*b/(r+xi)*xi/(eta+xi))
  )
}

#' @title Compute observed prevalence by light microscopy
#' @description Implements [F_prevalence] for the `SIP` model.
#' @inheritParams ramp.xds::F_pfpr_by_lm
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIP]
#' @export
F_pfpr_by_lm.SIP <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute observed prevalence by RDT
#' @description Implements [F_prevalence] for the `SIP` model.
#' @inheritParams ramp.xds::F_pfpr_by_rdt
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIP]
#' @export
F_pfpr_by_rdt.SIP <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute observed prevalence by PCR
#' @description Implements [F_prevalence] for the `SIP` model.
#' @inheritParams ramp.xds::F_pfpr_by_pcr
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIP]
#' @export
F_pfpr_by_pcr.SIP <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Check `SIP` (**XH**)
#'
#' @inheritParams ramp.xds::check_XH
#'
#' @return an **`xds`** model object
#' @keywords internal
#' @seealso [SIP]
#' @export
check_XH.SIP = function(xds_obj, i){
  return(xds_obj)
}

#' @title Get the skill set for `SIP` (**XH**)
#'
#' @description Return the capabilities and compatibilities of the `SIP` module as a list.
#'
#' @inheritParams ramp.xds::get_skillset_XH
#'
#' @return the **XH** module skill set, as a list
#' @keywords internal
#'
#' @seealso [SIP]
#' @export
get_skillset_XH.SIP = function(xds_obj,i){
  skills =   list(
    dts         = FALSE,
    xde         = TRUE,
    delay       = FALSE,
    demography  = TRUE,
    malaria     = TRUE,
    pr_obs      = TRUE,
    pf_lm       = TRUE,
    pf_rdt      = TRUE,
    pf_pcr      = TRUE,
    mda         = TRUE,
    msat        = TRUE
  )
  return(skills)
}
#' @title Compute steady states for `SIP` (**XH**)
#' @description Compute steady states given daily FoI `foi` and fixed host density `H`.
#' @note Delegate to [steady_state_X_SIP], which dispatches on `xds_obj$xds`.
#' @inheritParams ramp.xds::steady_state_X
#' @return a named [list] containing `H`, `I`, and `P` for `xde`, or a named
#' [numeric] vector containing `S`, `I`, and `P` for `dts`
#' @keywords internal
#' @export
steady_state_X.SIP = function(foi, H, xds_obj, i=1){
  steady_state_X_SIP(foi, H, xds_obj, i)
}

#' @rdname steady_state_X.SIP
#' @export
steady_state_X_SIP = function(foi, H, xds_obj, i=1){
  UseMethod("steady_state_X_SIP", xds_obj$xds)
}

#' @title Compute steady states for `SIP` (**XH**)
#' @description Compute steady-state infection densities for the `SIP` module
#' given daily FoI `foi` and fixed host population density `H` for the
#' differential-equation system.
#' @inheritParams steady_state_X_SIP
#' @return a named [list] with elements `H`, `I`, and `P`
#' @keywords internal
#' @seealso [SIP]
#' @export
steady_state_X_SIP.xde = function(foi, H, xds_obj, i=1){
  with(xds_obj$XH_obj[[i]],{
    Ieq = (foi*H*eta*(1-rho))/((foi+r+xi+sigma)*(eta+xi) +foi*((r-eta)*rho+sigma))
    Peq  = (H*xi*(foi+r+xi+sigma) + (foi*H*(r*rho + sigma)))/((foi+r+xi+sigma)*(eta+xi) +foi*((r-eta)*rho+sigma))
    return(list(H=H, I=as.vector(Ieq), P = as.vector(Peq)))
})}

#' @title Compute steady states for `SIP` (**XH**)
#' @description Compute steady-state infection densities for the `SIP` module
#' given daily FoI `foi` and fixed host population density `H` for the
#' discrete-time system.
#' @inheritParams steady_state_X_SIP
#' @return a named [numeric] vector containing `S`, `I`, and `P`
#' @keywords internal
#' @seealso [SIP]
#' @export
steady_state_X_SIP.dts = function(foi, H, xds_obj, i=1){
  ar = exp(-foi)
  with(xds_obj$XH_obj[[i]],{
    Iteq =(ar*H*eta*(1-rho))/((eta+xi)*(ar*(1-r)+(r+xi)+ar*(eta*(r-1)+r+xi)))
    Pteq = (H*(xi*(ar+r*(1-ar)+xi)+ar*(1+xi)*rho))/((eta+xi)*(ar*(1-r)+(r+xi)+ar*(eta*(r-1)+r+xi)))
    Steq = H-Iteq-Pteq
    return(c(S=Steq, I=Iteq,P =Pteq))
})}

#' Plot the density of infected individuals for the `SIP` model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @keywords internal
#' @seealso [SIP]
#' @export
xds_plot_X.SIP = function(xds_obj, i=1, clrs=c("darkblue", "darkred", "darkgreen"), llty=1, add=FALSE){
  XH = xds_obj$outputs$orbits$XH[[i]]
  time = xds_obj$outputs$time

  if(add==FALSE)
    plot(time, 0*time, type = "n", ylim = c(0, max(XH$H)),
         ylab = "# Infected", xlab = "Time")

  add_lines_X_SIP(time, XH, xds_obj$nStrata[i], clrs, llty)
}


#' Add lines for the density of infected individuals for the `SIP` model
#'
#' @param time time points for the observations
#' @param XH a list with the outputs of parse_XH_orbits.SIP
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @keywords internal
#' @seealso [SIP]
#' @export
add_lines_X_SIP = function(time, XH, nStrata, clrs=c("darkblue", "darkred", "darkgreen"), llty=1){
  if (length(llty)< nStrata) llty = rep(llty, nStrata)
  with(XH,{
    if(nStrata == 1){
      lines(time, S, col=clrs[1], lty = llty)
      lines(time, I, col=clrs[2], lty = llty)
      lines(time, P, col=clrs[3], lty = llty)
    } else {
      for(i in 1:nStrata)
        lines(time, S[,i], col=clrs[1], lty = llty[i])
      lines(time, I[,i], col=clrs[2], lty = llty[i])
      lines(time, P[,i], col=clrs[3], lty = llty[i])
    }})
}
