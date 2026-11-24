# specialized methods for the human SEISd model (delay SEIS)
# TODO: Review change_XH_inits.SEISd repairs: remove the extra species index
# inside inits and update H, E, I without adding derived S. Confirm cases
# remains unchanged, partial updates work, and other species are preserved.

#' @title The `SEISd` module for the **XH** component
#' @description
#' Implements the **XH** component using a delay Susceptible-Exposed-Infectious-Susceptible
#' (`SEISd`) model of human infection dynamics. This is a delay differential equation
#' version of the `SEIS` model, in which the exposed
#' period is modelled as a fixed delay \eqn{\nu} rather than an exponential rate.
#'
#' @section State Variables:
#' \describe{
#'   \item{`H`}{total human (or host) population density}
#'   \item{`E`}{density of exposed (infected but not yet infectious) humans}
#'   \item{`I`}{density of infectious humans}
#'   \item{`S=H-E-I`}{density of susceptible humans or hosts}
#'   \item{`cases`}{cumulative exposures used to track delayed progression}
#' }
#' This implementation includes human demographic change, progressive
#' stratification, and parameters to model detection by diagnostic methods.
#' It also includes parameters for treating and curing infections and functions
#' to model mass treatment.
#'
#' @inheritSection ramp.xds::SIS Demography
#'
#' @section Dynamical Term:
#' \describe{
#'   \item{\eqn{h} or `foi`}{The force of infection is computed and stored by [Exposure] as `FoI`}
#' }
#' In the `SEISd` module, there is no effect of pre-erythrocytic immunity, so
#' a constant fraction of infectious bites causes an infection: the parameter \eqn{b}.
#'
#' @section Parameters:
#' The parameters returned by [get_XH_pars.SEISd] and modifiable by
#' [change_XH_pars.SEISd] describe the dynamics or detection by diagnostic methods.
#' \describe{
#'   \item{`q`}{test positivity (default `0.8`)}
#'   \item{`b`}{transmission probability from mosquito to human}
#'   \item{`c`}{transmission probability from human to mosquito}
#'   \item{`nu`}{intrinsic incubation period (fixed delay, in days)}
#'   \item{`r`}{clearance rate for infections}
#'   \item{`rho` or \eqn{\rho}}{the fraction of incident infections that are treated and cured}
#'   \item{`sigma` or \eqn{\sigma}}{the rate that prevalent infections are treated and cured}
#'   \item{`q_lm`}{test positivity by light microscopy}
#'   \item{`q_rdt`}{test positivity by RDT}
#'   \item{`q_pcr`}{test positivity by PCR}
#' }
#'
#' @section Infection Dynamics:
#' Infections clear naturally at the per-capita rate \eqn{r}, after which
#' individuals are susceptible to infection.
#' Incident infections are delayed progression \eqn{C(t)}, not new exposures
#' \eqn{hS}. A fraction \eqn{\rho} of incident infections are treated and
#' cured, returning directly to susceptibility. Prevalent infections are treated
#' and cured at rate \eqn{\sigma}; mass treatment adds clearance at rate
#' \eqn{\Xi(t) = mda(t) + msat(t)}. Treatment does not change the exposed-state
#' equation or exposure history, and does not add a \eqn{\rho h} clearance term.
#' \deqn{
#' \begin{array}{rl}
#' dE/dt &= hS - C(t) - \mathbf{D}(E) \\
#' dI/dt &= (1-\rho)C(t) - (r+\sigma+\Xi(t))I - \mathbf{D}(I) \\
#' d\,\mathrm{cases}/dt &= hS \\
#' \end{array}}
#' where \eqn{h} is the force of infection, \eqn{S = H - E - I}, and
#' \eqn{C(t)} is zero for \eqn{t \leq \nu}; otherwise it is the
#' rate of exposure \eqn{h(t-\nu)S(t-\nu)} recorded \eqn{\nu} days ago
#' (computed via [deSolve::lagderiv]).
#' The mortality function uses the same model for \eqn{E} and \eqn{I} as for
#' \eqn{H}. Susceptible density is redundant and is computed after solving.
#'
#' @inheritSection SEIS Prevalence
#' @inheritSection SEIS Infectiousness
#'
#' @inheritSection ramp.xds::SIS Mass Treatment
#'
#' @note In `SEISd`, mass treatment clears infectious individuals, returning
#' them to susceptibility; exposed individuals and exposure history are unaffected.
#'
#' @seealso
#' + Basic Setup: [make_XH_obj_SEISd], [setup_XH_obj.SEISd]
#' + Dynamics: [dXHdt.SEISd]
#' + Parameters: [get_XH_pars.SEISd], [change_XH_pars.SEISd]
#' + Variables: [make_XH_inits_SEISd], [get_XH_inits], [change_XH_inits.SEISd], [setup_XH_inits.SEISd]
#' + [SEIS], [SIS]
#'
#' @name SEISd
#' @rdname SEISd
NULL

#' @title Compute the derivatives for `SEISd` (**XH**)
#' @inheritSection SEISd Demography
#' @inheritSection SEISd Infection Dynamics
#'
#' @description The derivatives for the `SEISd` module are:
#' \deqn{\begin{array}{rl}
#' dH/dt &= \mathbf{B}(t,y)-\mathbf{D}(t,H) \\
#' dE/dt &= hS-C_\nu(t)-\mathbf{D}(t,E) \\
#' dI/dt &= (1-\rho)C_\nu(t)-(r+\sigma+\Xi(t))I-\mathbf{D}(t,I) \\
#' d\,\mathrm{cases}/dt &= hS
#' \end{array}}
#' Here \eqn{h} is `foi`, \eqn{S=H-E-I}, and
#' \eqn{\Xi(t)=mda(t)+msat(t)}.
#' The delayed progression rate \eqn{C_\nu(t)} is zero for \eqn{t\leq\nu};
#' otherwise it is the derivative of `cases` at \eqn{t-\nu},
#' obtained with [deSolve::lagderiv].
#' The term \eqn{\mathbf{B}(t,y)} is returned by [ramp.xds::Births], and
#' \eqn{\mathbf{D}(t,X)} is returned by [ramp.xds::Deaths] for population \eqn{X}.
#' @inheritParams ramp.xds::dXHdt
#' @return a [numeric] vector
#' @importFrom deSolve lagvalue
#' @importFrom deSolve lagderiv
#' @keywords internal
#' @seealso [SEISd]
#' @export
dXHdt.SEISd <- function(t, y, xds_obj, i) {

  foi <- xds_obj$terms$FoI[[i]]

  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {

      if (t <= nu) {
        cases = 0
      } else {
        cases = lagderiv(t=t-nu, nr=ix$cases_ix)
      }

      dH <- Births(t, y, xds_obj, i) - Deaths(t, H, xds_obj, i)
      dE <- foi*S - cases - Deaths(t, E, xds_obj, i)
      dI <- (1-rho)*cases - (r+sigma)*I - Deaths(t, I, xds_obj, i)
      dI <- dI - mda(t)*I - msat(t)*I
      dcases <- foi*S

      return(c(dH, dE, dI, dcases))
    })
  })
}




#' @title Get the state variables for `SEISd` (**XH**)
#' @description Extract the state variables for the `SEISd` module
#' from `y` and return them by name in a list for use during model solving.
#' Also compute susceptible density `S = H - E - I` and extract `cases`.
#' @inheritParams ramp.xds::get_XH_vars
#' @return a named [list] of state variables and module-specific derived quantities
#' @keywords internal
#' @seealso [SEISd]
#' @export
#' @note Variable indices are stored at `xds_obj$XH_obj[[i]]$ix`.
get_XH_vars.SEISd <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
    H = y[H_ix]
    E = y[E_ix]
    I = y[I_ix]
    cases = y[cases_ix]
    S = H-E-I
    return(list(S=S,E=E,I=I,H=H,cases=cases))})
}

# specialized methods for the human SEISd model

#' @title Size of effective infectious human population
#' @description Implements [F_I] for the `SEISd` model.
#' @inheritParams ramp.xds::F_I
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SEISd]
#' @export
F_I.SEISd <- function(t, y, xds_obj, i) {
  I = y[xds_obj$XH_obj[[i]]$ix$I_ix]
  X = with(xds_obj$XH_obj[[i]], c*I)
  return(X)
}

#' @title Compute population density
#' @description Implements [F_H] for the `SEISd` model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SEISd]
#' @export
F_H.SEISd <- function(t, y, xds_obj, i){
  with(get_XH_vars(y, xds_obj, i), return(H))
}


#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_infectivity] for the `SEISd` model.
#' @inheritParams ramp.xds::F_infectivity
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SEISd]
#' @export
F_infectivity.SEISd <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]], b)
}

#' @title Make a module object for `SEISd` (**XH**)
#' @description Create an `SEISd` module object and initialize its parameters.
#' Inspect parameters with [get_XH_pars.SEISd] and change supported values
#' with [change_XH_pars.SEISd]. Values in `options` override parameter arguments.
#' @seealso [get_XH_pars.SEISd], [change_XH_pars.SEISd]
#' @inheritSection SEISd Parameters
#' @param nStrata is the number of population strata
#' @param options a [list] that could override defaults
#' @param b transmission probability (efficiency) from mosquito to human
#' @param r recovery rate
#' @param nu intrinsic incubation period (fixed delay, in days)
#' @param c transmission probability (efficiency) from human to mosquito
#' @param q_lm test positivity by light microscopy
#' @param q_rdt test positivity by RDT
#' @param q_pcr test positivity by PCR
#' @param rho the fraction of incident infections that are treated and cured
#' @param sigma the rate prevalent infections are treated and cured
#' @return an **XH** model object
#' @keywords internal
#' @seealso [SEISd]
#' @param q test positivity; a numeric scalar or input of length `nStrata`
#' @export
make_XH_obj_SEISd = function(nStrata, options=list(),
                             b=0.55, r=1/180, nu=14, c=0.15,
                             q_lm=0.8, q_rdt=0.8, q_pcr=0.9,
                             rho=0, sigma=0, q=0.8){
  with(options,{
    XH_obj = list()
    class(XH_obj) <- "SEISd"

    XH_obj$b = checkIt(b, nStrata)
    XH_obj$c = checkIt(c, nStrata)
    XH_obj$r = checkIt(r, nStrata)
    XH_obj$nu = checkIt(nu, nStrata)
    XH_obj$rho = checkIt(rho, nStrata)
    XH_obj$sigma = checkIt(sigma, nStrata)
    XH_obj$q_lm = checkIt(q_lm, nStrata)
    XH_obj$q_rdt = checkIt(q_rdt, nStrata)
    XH_obj$q_pcr = checkIt(q_pcr, nStrata)

    XH_obj$q = checkIt(q, nStrata)

    return(XH_obj)
  })}

#' @title Set up the module object for `SEISd` (**XH**)
#' @description
#'
#' Set up an `SEISd` module object for the **XH** component and the \eqn{i^{th}} host species.
#'
#' + Set the module name to `SEISd`.
#' + Create the module object with [make_XH_obj_SEISd], passing `nStrata[i]` and `options`.
#' + Set initial values using [ramp.xds::setup_XH_inits], passing `options` to override defaults.
#' + Initialize mass-treatment functions using [ramp.xds::setup_mass_treatment].
#' + Initialize the time-spent matrix using [ramp.xds::setup_timespent.setup].
#' + Initialize birth and death rates using [ramp.xds::setup_demog_obj].
#' + Initialize blood search weights using [ramp.xds::setup_blood_search_weights.setup].
#' + Initialize time away using [ramp.xds::setup_time_away.setup].
#' + Initialize travel EIR using [ramp.xds::setup_travel_eir.no_travel].
#' + Return the updated **`xds`** model object.
#'
#' @seealso [SEISd], [make_XH_obj_SEISd]
#'
#' @note The module object is stored at `xds_obj$XH_obj[[i]]`
#' @inheritSection SEISd State Variables
#' @inheritSection SEISd Parameters
#'
#' @inheritParams ramp.xds::setup_XH_obj
#' @return an **`xds`** model object
#' @keywords internal
#' @export
setup_XH_obj.SEISd = function(XHname, residence, xds_obj, options=list(), i=1){
  xds_obj$XHname = "SEISd"
  xds_obj$XH_obj[[i]] = make_XH_obj_SEISd(xds_obj$nStrata[i], options)
  xds_obj <- setup_XH_inits(xds_obj, options=options, i=i)
  xds_obj <- setup_mass_treatment(xds_obj, i=i)
  xds_obj <- setup_timespent("setup", xds_obj, list(residence=residence), i)
  xds_obj <- setup_demog_obj(xds_obj, i)
  xds_obj <- setup_blood_search_weights("setup", xds_obj, i=i)
  xds_obj <- setup_time_away("setup", xds_obj, i=i)
  xds_obj <- setup_travel_eir("no_travel", xds_obj, i=i)
  return(xds_obj)
}


#' @title Get parameters for `SEISd` (**XH**)
#' @seealso [make_XH_obj_SEISd], [change_XH_pars.SEISd]
#' @description Return the parameters for the `SEISd` module as a named list.
#' @note Parameter values for the \eqn{i^{th}} host are
#' stored as `xds_obj$XH_obj[[i]]`.
#' @inheritParams ramp.xds::get_XH_pars
#' @inheritSection SEISd Parameters
#' @return a named [list] of parameter values
#' @seealso [SEISd], [make_XH_obj_SEISd]
#' @keywords internal
#' @export
get_XH_pars.SEISd <- function(xds_obj, i=1) {
  with(xds_obj$XH_obj[[i]],list(q=q, b=b, c=c, r=r, nu=nu,
                                q_lm=q_lm, q_rdt=q_rdt, q_pcr=q_pcr,
                                rho=rho, sigma=sigma))
}


#' @title Change parameters for `SEISd` (**XH**)
#' @seealso [make_XH_obj_SEISd], [get_XH_pars.SEISd]
#' @description Change parameter values named in `options` for the `SEISd`
#' module, leaving other values unchanged.
#' @inheritParams ramp.xds::change_XH_pars
#' @inheritSection SEISd Parameters
#' @return an **`xds`** model object
#' @keywords internal
#' @seealso [SEISd]
#' @export
change_XH_pars.SEISd <- function(xds_obj, options=list(), i=1) {
  with(xds_obj$XH_obj[[i]], with(options,{
    xds_obj$XH_obj[[i]]$b <- b
    xds_obj$XH_obj[[i]]$c <- c
    xds_obj$XH_obj[[i]]$r <- r
    xds_obj$XH_obj[[i]]$nu <- nu
    xds_obj$XH_obj[[i]]$rho <- rho
    xds_obj$XH_obj[[i]]$sigma <- sigma
    xds_obj$XH_obj[[i]]$q_lm <- q_lm
    xds_obj$XH_obj[[i]]$q_rdt <- q_rdt
    xds_obj$XH_obj[[i]]$q_pcr <- q_pcr
    xds_obj$XH_obj[[i]]$q <- q
    return(xds_obj)
  }))}

#' @title Make initial values for `SEISd` (**XH**)
#' @description Create a named list of initial values for the `SEISd` module.
#' Values in `options` override the defaults.
#' @details Initial values must be numeric scalars or have length `nStrata`.
#' Scalars are expanded to that length.
#' Initialize `cases` to zero in each stratum. This accessory value is constructed
#' internally rather than independently set through `options`.
#' @inheritSection SEISd State Variables
#' @param nStrata the number of strata in the model
#' @param options a [list] to override defaults
#' @param H the initial human population density
#' @param E the initial number of exposed humans
#' @param I the initial number of infectious humans
#' @return a named [list] of initial values
#' @keywords internal
#' @seealso [SEISd]
#' @export
#' @seealso [setup_XH_inits.SEISd], [change_XH_inits.SEISd]
make_XH_inits_SEISd = function(nStrata, options = list(), H=1000, E=0, I=1){
  with(options,{
    H = checkIt(H, nStrata)
    E = checkIt(E, nStrata)
    I = checkIt(I, nStrata)
    cases = rep(0, nStrata)
    return(list(H=H, E=E, I=I, cases=cases))
})}

#' @title Set up initial values for `SEISd` (**XH**)
#' @description Set up initial values for the `SEISd` module using
#' [make_XH_inits_SEISd]. Values in `options` override the defaults.
#' @inheritParams ramp.xds::setup_XH_inits
#' @inheritSection SEISd State Variables
#' @return an **`xds`** model object
#' @keywords internal
#' @seealso [SEISd]
#' @export
#' @seealso [change_XH_inits.SEISd]
#' @note Initial values are stored at `xds_obj$XH_obj[[i]]$inits`.
#' @seealso [ramp.xds::get_XH_inits]
setup_XH_inits.SEISd = function(xds_obj, options=list(), i=1){
  xds_obj$XH_obj[[i]]$inits = make_XH_inits_SEISd(xds_obj$nStrata[i], options)
  return(xds_obj)
}


#' @title Change initial values for `SEISd` (**XH**)
#' @description Change initial values named in `options` for the `SEISd`
#' module, leaving other values unchanged.
#' Supported values are `H`, `E`, and `I`; `cases` is left unchanged.
#' @inheritParams ramp.xds::change_XH_inits
#' @inheritSection SEISd State Variables
#' @return an **`xds`** model object
#' @keywords internal
#' @seealso [SEISd]
#' @export
#' @seealso [make_XH_inits_SEISd], [setup_XH_inits.SEISd]
#' @note Stored initial values for this module are at
#' `xds_obj$XH_obj[[i]]$inits`.
change_XH_inits.SEISd <- function(xds_obj, options=list(), i=1) {
  with(xds_obj$XH_obj[[i]]$inits,
    with(options,{
      xds_obj$XH_obj[[i]]$inits$H = H
      xds_obj$XH_obj[[i]]$inits$E = E
      xds_obj$XH_obj[[i]]$inits$I = I
      return(xds_obj)
}))}

#' @title Set up variable indices for `SEISd` (**XH**)
#' @description Assign consecutive variable indices for the `SEISd` module,
#' starting after `xds_obj$max_ix`, and update the maximum index.
#' Allocate blocks in the order `H_ix`, `E_ix`, `I_ix`, `cases_ix`.
#' Each block has length `nStrata[i]`.
#' @inheritParams ramp.xds::setup_XH_ix
#' @return an **`xds`** model object
#' @importFrom utils tail
#' @keywords internal
#' @seealso [SEISd]
#' @export
#' @note Indices are stored at `xds_obj$XH_obj[[i]]$ix`.
#' The maximum index is stored at `xds_obj$max_ix`.
setup_XH_ix.SEISd <- function(xds_obj, i) {with(xds_obj,{

  H_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(H_ix, 1)

  E_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(E_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(I_ix, 1)

  cases_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(cases_ix, 1)

  xds_obj$max_ix = max_ix
  xds_obj$XH_obj[[i]]$ix = list(H_ix=H_ix, E_ix=E_ix, I_ix=I_ix,
                                cases_ix=cases_ix)
  return(xds_obj)
})}


#' @title Parse outputs for `SEISd` (**XH**)
#' @description Extract the model state outputs for the `SEISd` module
#' from `outputs` and return them by name in a list.
#' Also compute `S = H-E-I`, `ni` using [ramp.xds::F_ni], and `pr`
#' using [ramp.xds::F_prevalence]. The `cases` state is not returned.
#' @note Variable indices are stored at `xds_obj$XH_obj[[i]]$ix`.
#' @inheritParams ramp.xds::parse_XH_orbits
#' @return a [list]
#' @keywords internal
#' @seealso [SEISd]
#' @export
parse_XH_orbits.SEISd <- function(outputs, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix, {
    H = outputs[,H_ix]
    E = outputs[,E_ix]
    I = outputs[,I_ix]
    S = H-E-I
    vars <- list(S=S, E=E, I=I, H=H)
    vars$ni <- F_ni(vars, xds_obj$XH_obj[[i]])
    vars$pr <- F_prevalence(vars, xds_obj$XH_obj[[i]])
    return(vars)
})}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_prevalence] for the `SEISd` model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SEISd]
#' @export
F_prevalence.SEISd <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute observed prevalence
#' @description Return observed prevalence for the `SEISd` module as
#' `q` times the true prevalence returned by [ramp.xds::F_prevalence].
#' @inheritParams ramp.xds::F_pr_obs
#' @return a [numeric] vector of length `nStrata`
#' @seealso [SEISd]
#' @keywords internal
#' @export
F_pr_obs.SEISd <- function(vars, XH_obj) {
  XH_obj$q * ramp.xds::F_prevalence(vars, XH_obj)
}

#' @title Compute the NI
#' @description Implements [F_ni] for the `SEISd` model.
#' @inheritParams ramp.xds::F_ni
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SEISd]
#' @export
F_ni.SEISd <- function(vars, XH_obj) {
  ni = with(vars, XH_obj$c*I/H)
  return(ni)
}

#' @title Compute the HTC for the `SEISd` model
#' @description Implements [get_HTC] for the `SEISd` model with demography.
#' @inheritParams ramp.xds::get_HTC
#' @return a [numeric] vector
#' @keywords internal
#' @seealso [SEISd]
#' @export
get_HTC.SEISd <- function(xds_obj, i) {
  with(xds_obj$XH_obj[[i]],
       return(c/r)
  )
}

#' @title Compute observed prevalence by light microscopy
#' @description Implements [F_prevalence] for the `SEISd` model.
#' @inheritParams ramp.xds::F_pfpr_by_lm
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SEISd]
#' @export
F_pfpr_by_lm.SEISd <- function(vars, XH_obj) {
  pr = with(XH_obj, with(vars, q_lm*I/H))
  return(pr)
}

#' @title Compute observed prevalence by RDT
#' @description Implements [F_prevalence] for the `SEISd` model.
#' @inheritParams ramp.xds::F_pfpr_by_rdt
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SEISd]
#' @export
F_pfpr_by_rdt.SEISd <- function(vars, XH_obj) {
  pr = with(XH_obj, with(vars, q_rdt*I/H))
  return(pr)
}

#' @title Compute observed prevalence by PCR
#' @description Implements [F_prevalence] for the `SEISd` model.
#' @inheritParams ramp.xds::F_pfpr_by_pcr
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SEISd]
#' @export
F_pfpr_by_pcr.SEISd <- function(vars, XH_obj) {
  pr = with(XH_obj, with(vars, q_pcr*I/H))
  return(pr)
}

#' @title Check `SEISd` (**XH**)
#'
#' @inheritParams ramp.xds::check_XH
#'
#' @return an **`xds`** model object
#' @keywords internal
#' @seealso [SEISd]
#' @export
check_XH.SEISd = function(xds_obj, i){
  return(xds_obj)
}

#' @title Get the skill set for `SEISd` (**XH**)
#'
#' @description Return the capabilities and compatibilities of the `SEISd` module as a list.
#'
#' @inheritParams ramp.xds::get_skillset_XH
#'
#' @return the **XH** module skill set, as a list
#' @keywords internal
#'
#' @seealso [SEISd]
#' @export
get_skillset_XH.SEISd = function(xds_obj,i){
  skills =   list(
    dts         = FALSE,
    xde         = TRUE,
    delay       = TRUE,
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


#' @title Compute steady states for `SEISd` (**XH**)
#' @description Compute steady-state infection densities for the `SEISd` module
#' given daily FoI `foi` and fixed host population density `H` for the
#' delay differential-equation system.
#' @note Assumes no demographic change or mass treatment and constant exposure
#' history after the delay. Includes `rho` and `sigma`.
#' @inheritParams ramp.xds::steady_state_X
#' @return a named [numeric] vector containing `S`, `E`, `I`, and `H`
#' @keywords internal
#' @seealso [SEISd]
#' @export
steady_state_X.SEISd = function(foi, H, xds_obj, i=1){
  with(xds_obj$XH_obj[[i]],{
    Seq = H/(1+foi*nu + (1-rho)*foi/(r+sigma))
    Eeq = foi*Seq*nu
    Ieq = (1-rho)*foi*Seq/(r+sigma)
    return(c(S=Seq, E = Eeq, I=Ieq, H=H))
})}

#' Plot the density of infected individuals for the `SEISd` model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @keywords internal
#' @seealso [SEISd]
#' @export
xds_plot_X.SEISd = function(xds_obj, i=1, clrs=c("darkblue","darkred"), llty=1, add=FALSE){
  XH = xds_obj$outputs$orbits$XH[[i]]
  time = xds_obj$outputs$time

  if(add==FALSE)
    plot(time, 0*time, type = "n", ylim = c(0, max(XH$H)),
         ylab = "# Infected", xlab = "Time")

  add_lines_X_SIS(time, XH, xds_obj$nStrata[i], clrs, llty)
}

#' Add lines for the density of infected individuals for the `SEISd` model
#'
#' @param time time points for the observations
#' @param XH a list with the outputs of parse_XH_orbits_SEISd
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @keywords internal
#' @seealso [SEISd]
#' @export
add_lines_X_SEISd = function(time, XH, nStrata, clrs=c("darkblue","darkred"), llty=1){
  if (length(llty)< nStrata) llty = rep(llty, nStrata)
  with(XH,{
    if(nStrata == 1){
      lines(time, S+E, col=clrs[1], lty = llty)
      lines(time, I, col=clrs[2], lty = llty)
    } else {
      for(i in 1:nStrata)
        lines(time, S[,i] + E[,i], col=clrs[1], lty = llty[i])
      lines(time, I[,i], col=clrs[2], lty = llty[i])
    }})
}
