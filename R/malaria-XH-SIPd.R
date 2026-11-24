# specialized methods for the human SIPd model


# TODO: Review steady-state dispatch repair: separate SIPd's formerly duplicate
# definitions into xde/dts methods; dts converts foi to ar using exp(-foi).
# TODO: Review parse_XH_orbits.SIPd repair: read indices from XH_obj[[i]]$ix
# instead of the legacy X_obj path.
# TODO: Review change_XH_inits.SIPd repairs: remove the extra species index
# inside inits and update stored H, I, P instead of derived S. Confirm H
# update semantics, preservation of cfoi, partial updates, and species isolation.

#' @title The `SIPd` module for the **XH** component
#' @description
#' Implements the **XH** component using a delay Susceptible-Infected-chemoProtected
#' (`SIPd`) model of malaria infection dynamics. This is the delay differential
#' equation version of the `SIP` model, in which the duration of chemoprotection
#' \eqn{\eta} is modelled as a fixed delay.
#' Treatment clears infection and moves individuals into \eqn{P}, where they
#' are protected from infection until chemoprotection ends.
#' This implementation includes treatment parameters, human demographic change,
#' and support for progressive stratification.
#'
#' @section State Variables:
#' \describe{
#'   \item{`H`}{total human (or host) population density}
#'   \item{`I`}{density of infectious humans}
#'   \item{`P`}{density of humans under chemoprotection}
#'   \item{`S=H-I-P`}{density of susceptible humans or hosts}
#'   \item{`treated`}{cumulative treatment history used for delayed outflow from \eqn{P}}
#' }
#' The delay derivative tracks treatment history. The initializer currently
#' returns this auxiliary initial value under the name `cfoi`.
#'
#' @inheritSection ramp.xds::SIS Demography
#'
#' @section Dynamical Term:
#' \describe{
#'   \item{\eqn{h} or `foi`}{The force of infection is computed and stored by [Exposure] as `FoI`}
#' }
#' In the `SIPd` module, a constant fraction of infectious bites causes an
#' infection: the parameter \eqn{b}.
#'
#' @section Parameters:
#' The parameters returned by [get_XH_pars.SIPd] and modifiable by
#' [change_XH_pars.SIPd] describe infection, treatment, and chemoprotection.
#' \describe{
#'   \item{`q`}{test positivity (default `0.8`)}
#'   \item{`b`}{transmission probability from mosquito to human}
#'   \item{`c`}{transmission probability from human to mosquito}
#'   \item{`r`}{clearance rate for infections}
#'   \item{`rho` or \eqn{\rho}}{fraction of incident infections treated and cured, entering \eqn{P}}
#'   \item{`xi` or \eqn{\xi}}{background treatment rate for susceptible and infectious individuals}
#'   \item{`eta`}{chemoprotection duration (fixed delay, in days)}
#' }
#'
#' @section Infection Dynamics:
#' Infections clear naturally at rate \eqn{r}, returning individuals to
#' susceptibility. Treatment moves individuals into chemoprotection.
#' A fraction \eqn{\rho} of incident infections \eqn{hS} are treated immediately.
#' \deqn{
#' \begin{array}{rl}
#' dI/dt &= (1-\rho)hS - (r + \xi)I - \mathbf{D}(I) \\
#' dP/dt &= \rho hS + \xi(S+I) - T(t) - \mathbf{D}(P) \\
#' d\,\mathrm{treated}/dt &= \rho hS + \xi(S+I) \\
#' \end{array}}
#' where \eqn{S=H-I-P}. Delayed outflow \eqn{T(t)} is zero for
#' \eqn{t \leq \eta}; otherwise it is the treatment inflow recorded
#' \eqn{\eta} days ago, computed via [deSolve::lagderiv].
#' The mortality function uses the same model for \eqn{I} and \eqn{P} as for \eqn{H}.
#'
#' @section Prevalence:
#' [F_pr_obs.SIPd] returns observed prevalence \eqn{x_q=qI/H}.
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
#' @section Mass Treatment:
#' This implementation does not currently initialize or use `mda(t)` or
#' `msat(t)`. Treatment is modelled through the parameters `rho` and `xi`;
#' unlike `SIP`, it has no `sigma` parameter.
#'
#' @seealso
#' + Basic Setup: [make_XH_obj_SIPd], [setup_XH_obj.SIPd]
#' + Dynamics: [dXHdt.SIPd]
#' + Parameters: [get_XH_pars.SIPd], [change_XH_pars.SIPd]
#' + Variables: [make_XH_inits_SIPd], [get_XH_inits], [change_XH_inits.SIPd], [setup_XH_inits.SIPd]
#' + [SIP], [SIS]
#'
#' @name SIPd
#' @rdname SIPd
NULL

#' @title Compute the derivatives for `SIPd` (**XH**)
#' @description The derivatives for the `SIPd` module are:
#' \deqn{\begin{array}{rl}
#' dH/dt &= \mathbf{B}(t,y)-\mathbf{D}(t,H) \\
#' dI/dt &= (1-\rho)hS-(r+\xi)I-\mathbf{D}(t,I) \\
#' dP/dt &= \rho hS+\xi(S+I)-T_\eta(t)-\mathbf{D}(t,P) \\
#' d\,\mathrm{treated}/dt &= (\rho h+\xi)S+\xi I
#' \end{array}}
#' Here \eqn{h} is `foi` and \eqn{S} is the susceptible density
#' returned by [get_XH_vars.SIPd].
#' The delayed treatment rate \eqn{T_\eta(t)} is zero for \eqn{t\leq\eta};
#' otherwise it is the derivative of `treated` at \eqn{t-\eta},
#' obtained with [deSolve::lagderiv].
#' The term \eqn{\mathbf{B}(t,y)} is returned by [ramp.xds::Births], and
#' \eqn{\mathbf{D}(t,X)} is returned by [ramp.xds::Deaths] for population \eqn{X}.
#' @inheritSection SIPd Demography
#' @inheritSection SIPd Infection Dynamics
#'
#' @seealso [SIPd], [make_XH_obj_SIPd]
#' @inheritParams ramp.xds::dXHdt
#' @return a [numeric] vector
#' @keywords internal
#' @export
dXHdt.SIPd <- function(t, y, xds_obj, i){

  foi <- xds_obj$terms$FoI[[i]]

  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {
      if (t <= eta) {
        treated_eta = 0
      } else {
        treated_eta = lagderiv(t=t-eta, nr=xds_obj$XH_obj[[i]]$ix$treated_ix)
      }

      dH <- Births(t, y, xds_obj, i) - Deaths(t, H, xds_obj, i)
      dI <- (1-rho)*foi*S - (r+xi)*I - Deaths(t, I, xds_obj, i)
      dP <- rho*foi*S + xi*(S+I) - treated_eta - Deaths(t, P, xds_obj, i)
      dtreated <- (rho*foi+xi)*S + xi*I

      return(c(dH, dI, dP, dtreated))
    })
  })
}

#' @title Update states for `SIPd` (**XH**)
#' @description Implements [Update_XHt] for the `SIPd` model.
#' @inheritParams ramp.xds::Update_XHt
#' @return a [numeric] vector
#' @keywords internal
#' @seealso [SIPd]
#' @export
Update_XHt.SIPd <- function(t, y, xds_obj, i){

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



#' @title Get the state variables for `SIPd` (**XH**)
#' @description Extract the state variables for the `SIPd` module
#' from `y` and return them by name in a list for use during model solving.
#' Also compute total host density `H = S + I + P`.
#' @inheritParams ramp.xds::get_XH_vars
#' @return a named [list] of state variables and module-specific derived quantities
#' @keywords internal
#' @seealso [SIPd]
#' @export
#' @note Variable indices are stored at `xds_obj$XH_obj[[i]]$ix`.
get_XH_vars.SIPd <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
    S = y[S_ix]
    I = y[I_ix]
    P = y[P_ix]
    H = S+I+P
    return(list(S=S,I=I,P=P,H=H))})
}

#' @title Size of effective infectious human population
#' @description Implements [F_I] for the `SIPd` model.
#' @inheritParams ramp.xds::F_I
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIPd]
#' @export
F_I.SIPd <- function(t, y, xds_obj, i) {
  I = y[xds_obj$XH_obj[[i]]$ix$I_ix]
  X = with(xds_obj$XH_obj[[i]], c*I)
  return(X)
}

#' @title Compute population density
#' @description Implements [F_H] for the `SIPd` model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIPd]
#' @export
F_H.SIPd <- function(t, y, xds_obj, i){
  with(get_XH_vars(y, xds_obj, i), return(H))
}


#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_infectivity] for the `SIPd` model.
#' @inheritParams ramp.xds::F_infectivity
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIPd]
#' @export
F_infectivity.SIPd <- function(y, xds_obj,i) {
  with(xds_obj$XH_obj[[i]], b)
}


#' @title Make a module object for `SIPd` (**XH**)
#' @description Create an `SIPd` module object and initialize its parameters.
#' Inspect parameters with [get_XH_pars.SIPd] and change supported values
#' with [change_XH_pars.SIPd]. Values in `options` override parameter arguments.
#' @seealso [get_XH_pars.SIPd], [change_XH_pars.SIPd]
#' @inheritSection SIPd Parameters
#' @param nStrata the number of population strata
#' @param options a [list] that could override defaults
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @param rho fraction of incident infections treated and cured, entering P
#' @param eta chemoprotection duration (fixed delay, in days)
#' @param xi background treatment rate for susceptible and infectious individuals
#' @return an **XH** model object
#' @keywords internal
#' @seealso [SIPd]
#' @param q test positivity; a numeric scalar or input of length `nStrata`
#' @export
make_XH_obj_SIPd = function(nStrata, options=list(),
                            b=0.55, r=1/180, c=0.15,
                            rho=.1, eta=25, xi=1/365, q=0.8){
  with(options,{
    XH_obj = list()
    class(XH_obj) <- c("SIPd")

    XH_obj$b = checkIt(b, nStrata)
    XH_obj$c = checkIt(c, nStrata)
    XH_obj$r = checkIt(r, nStrata)
    XH_obj$rho = checkIt(rho, nStrata)
    XH_obj$eta = checkIt(eta, nStrata)
    XH_obj$xi = checkIt(xi, nStrata)



    XH_obj$q = checkIt(q, nStrata)

    return(XH_obj)
  })}

#' @title Set up the module object for `SIPd` (**XH**)
#' @description
#'
#' Set up an `SIPd` module object for the **XH** component and the \eqn{i^{th}} host species.
#'
#' + Set the module name to `SIPd`.
#' + Create the module object with [make_XH_obj_SIPd], passing `nStrata[i]` and `options`.
#' + Set initial values using [ramp.xds::setup_XH_inits], passing `options` to override defaults.
#' + Initialize the time-spent matrix using [ramp.xds::setup_timespent.setup].
#' + Initialize birth and death rates using [ramp.xds::setup_demog_obj].
#' + Initialize blood search weights using [ramp.xds::setup_blood_search_weights.setup].
#' + Initialize time away using [ramp.xds::setup_time_away.setup].
#' + Initialize travel EIR using [ramp.xds::setup_travel_eir.no_travel].
#' + Return the updated **`xds`** model object.
#'
#' @seealso [SIPd], [make_XH_obj_SIPd]
#'
#' @note The module object is stored at `xds_obj$XH_obj[[i]]`
#' @inheritSection SIPd State Variables
#' @inheritSection SIPd Parameters
#'
#' @inheritParams ramp.xds::setup_XH_obj
#' @return an **`xds`** model object
#' @keywords internal
#' @export
setup_XH_obj.SIPd = function(XHname, residence, xds_obj, options=list(), i=1){
  xds_obj$XHname = "SIPd"
  xds_obj$XH_obj[[i]] = make_XH_obj_SIPd(xds_obj$nStrata[i], options)
  xds_obj <- setup_XH_inits(xds_obj, options=options, i=i)
  xds_obj <- setup_timespent("setup", xds_obj, list(residence=residence), i)
  xds_obj <- setup_demog_obj(xds_obj, i)
  xds_obj <- setup_blood_search_weights("setup", xds_obj, i=i)
  xds_obj <- setup_time_away("setup", xds_obj, i=i)
  xds_obj <- setup_travel_eir("no_travel", xds_obj, i=i)
  return(xds_obj)
}


#' @title Get parameters for `SIPd` (**XH**)
#' @seealso [make_XH_obj_SIPd], [change_XH_pars.SIPd]
#' @description Return the parameters for the `SIPd` module as a named list.
#' @note Parameter values for the \eqn{i^{th}} host are
#' stored as `xds_obj$XH_obj[[i]]`.
#' @inheritParams ramp.xds::get_XH_pars
#' @inheritSection SIPd Parameters
#' @return a named [list] of parameter values
#' @seealso [SIPd], [make_XH_obj_SIPd]
#' @keywords internal
#' @export
get_XH_pars.SIPd <- function(xds_obj, i=1) {
  with(xds_obj$XH_obj[[i]],list(q=q, b=b, c=c, r=r, rho=rho, eta=eta, xi=xi))
}


#' @title Change parameters for `SIPd` (**XH**)
#' @seealso [make_XH_obj_SIPd], [get_XH_pars.SIPd]
#' @description Change parameter values named in `options` for the `SIPd`
#' module, leaving other values unchanged.
#' @inheritParams ramp.xds::change_XH_pars
#' @inheritSection SIPd Parameters
#' @return an **`xds`** model object
#' @keywords internal
#' @seealso [SIPd]
#' @export
change_XH_pars.SIPd <- function(xds_obj, options=list(), i=1) {
  nHabitats <- xds_obj$nHabitats
  with(xds_obj$XH_obj[[i]], with(options,{
    xds_obj$XH_obj[[i]]$b <- b
    xds_obj$XH_obj[[i]]$c <- c
    xds_obj$XH_obj[[i]]$r <- r
    xds_obj$XH_obj[[i]]$rho <- rho
    xds_obj$XH_obj[[i]]$eta <- eta
    xds_obj$XH_obj[[i]]$xi <- xi
    xds_obj$XH_obj[[i]]$q <- q
    return(xds_obj)
  }))}

#' @title Make initial values for `SIPd` (**XH**)
#' @description Create a named list of initial values for the `SIPd` module.
#' Values in `options` override the defaults.
#' @details Initial values must be numeric scalars or have length `nStrata`.
#' Scalars are expanded to that length.
#' Initialize `cfoi` to zero in each stratum. This accessory value is constructed
#' internally rather than independently set through `options`.
#' @inheritSection SIPd State Variables
#' @param nStrata the number of population strata
#' @param options a [list] that could override defaults
#' @param H the initial human population density
#' @param I the initial number of infectious humans
#' @param P the initial number of humans under chemoprotection
#' @return a named [list] of initial values
#' @keywords internal
#' @seealso [SIPd]
#' @export
#' @seealso [setup_XH_inits.SIPd], [change_XH_inits.SIPd]
make_XH_inits_SIPd = function(nStrata,
                              options = list(),
                              H=1000, I=1, P=0){
  with(options,{
    H = checkIt(H, nStrata)
    I = checkIt(I, nStrata)
    P = checkIt(P, nStrata)
    return(list(H=H, I=I, P=P, cfoi=P*0))
})}


#' @title Set up initial values for `SIPd` (**XH**)
#' @description Set up initial values for the `SIPd` module using
#' [make_XH_inits_SIPd]. Values in `options` override the defaults.
#' @inheritParams ramp.xds::setup_XH_inits
#' @inheritSection SIPd State Variables
#' @return an **`xds`** model object
#' @keywords internal
#' @seealso [SIPd]
#' @export
#' @seealso [change_XH_inits.SIPd]
#' @note Initial values are stored at `xds_obj$XH_obj[[i]]$inits`.
#' @seealso [ramp.xds::get_XH_inits]
setup_XH_inits.SIPd = function(xds_obj, options=list(), i=1){
  xds_obj$XH_obj[[i]]$inits = make_XH_inits_SIPd(xds_obj$nStrata[i], options)
  return(xds_obj)
}


#' @title Change initial values for `SIPd` (**XH**)
#' @description Change initial values named in `options` for the `SIPd`
#' module, leaving other values unchanged.
#' Supported values are `H`, `I`, and `P`; `cfoi` is left unchanged.
#' @inheritParams ramp.xds::change_XH_inits
#' @inheritSection SIPd State Variables
#' @return an **`xds`** model object
#' @keywords internal
#' @seealso [SIPd]
#' @export
#' @seealso [make_XH_inits_SIPd], [setup_XH_inits.SIPd]
#' @note Stored initial values for this module are at
#' `xds_obj$XH_obj[[i]]$inits`.
change_XH_inits.SIPd <- function(xds_obj, options=list(), i=1) {
  with(get_XH_inits(xds_obj, i), with(options,{
    xds_obj$XH_obj[[i]]$inits$H = H
    xds_obj$XH_obj[[i]]$inits$I = I
    xds_obj$XH_obj[[i]]$inits$P = P
    return(xds_obj)
  }))}

#' @title Set up variable indices for `SIPd` (**XH**)
#' @description Assign consecutive variable indices for the `SIPd` module,
#' starting after `xds_obj$max_ix`, and update the maximum index.
#' Allocate blocks in the order `S_ix`, `I_ix`, `P_ix`, `treated_ix`.
#' Each block has length `nStrata[i]`.
#' @inheritParams ramp.xds::setup_XH_ix
#' @return an **`xds`** model object
#' @importFrom utils tail
#' @keywords internal
#' @seealso [SIPd]
#' @export
#' @note Indices are stored at `xds_obj$XH_obj[[i]]$ix`.
#' The maximum index is stored at `xds_obj$max_ix`.
setup_XH_ix.SIPd <- function(xds_obj, i) {with(xds_obj,{

  S_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(S_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(I_ix, 1)

  P_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(P_ix, 1)

  treated_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(treated_ix, 1)

  xds_obj$max_ix = max_ix
  xds_obj$XH_obj[[i]]$ix = list(S_ix=S_ix, I_ix=I_ix, P_ix=P_ix, treated_ix=treated_ix)
  return(xds_obj)
})}

#' @title Parse outputs for `SIPd` (**XH**)
#' @description Extract the model state outputs for the `SIPd` module
#' from `outputs` and return them by name in a list.
#' Also compute `H = S+I+P`, `ni` using [ramp.xds::F_ni], and `true_pr`
#' using [ramp.xds::F_prevalence].
#' @note Variable indices are stored at `xds_obj$XH_obj[[i]]$ix`.
#' The output-processing pipeline stores the returned orbits at
#' `xds_obj$outputs$orbits$XH[[i]]`.
#'
#' @inheritParams ramp.xds::parse_XH_orbits
#' @return a [list]
#' @keywords internal
#' @seealso [SIPd]
#' @export
parse_XH_orbits.SIPd <- function(outputs, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
    S = outputs[,S_ix]
    I = outputs[,I_ix]
    P = outputs[,P_ix]
    H = S+I+P
    vars <- list(S=S, I=I, P=P, H=H)
    vars$ni <- F_ni(vars, xds_obj$XH_obj[[i]])
    vars$true_pr <- F_prevalence(vars, xds_obj$XH_obj[[i]])
    return(vars)
  })}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_prevalence] for the `SIPd` model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIPd]
#' @export
F_prevalence.SIPd <- function(vars, XH_obj){
  pr = with(vars, I/H)
  return(pr)
}


#' @title Compute observed prevalence
#' @description Return observed prevalence for the `SIPd` module as
#' `q` times the true prevalence returned by [ramp.xds::F_prevalence].
#' @inheritParams ramp.xds::F_pr_obs
#' @return a [numeric] vector of length `nStrata`
#' @seealso [SIPd]
#' @keywords internal
#' @export
F_pr_obs.SIPd <- function(vars, XH_obj) {
  XH_obj$q * ramp.xds::F_prevalence(vars, XH_obj)
}

#' @title Compute the net infectiousness
#' @description Implements [F_ni] for the `SIPd` model.
#' @inheritParams ramp.xds::F_ni
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIPd]
#' @export
F_ni.SIPd <- function(vars, XH_obj) {
  ni = with(vars, XH_obj$c*I/H)
  return(ni)
}


#' @title Compute the HTC for the `SIPd` model
#' @description Implements [get_HTC] for the `SIPd` model with demography.
#' @inheritParams ramp.xds::get_HTC
#' @return a [numeric] vector
#' @keywords internal
#' @seealso [SIPd]
#' @export
get_HTC.SIPd <- function(xds_obj, i) {
  with(xds_obj$XH_obj[[i]],
       return((1-rho)*b/(r+xi)*xi/(eta+xi))
  )
}

#' @title Compute observed prevalence by light microscopy
#' @description Implements [F_prevalence] for the `SIPd` model.
#' @inheritParams ramp.xds::F_pfpr_by_lm
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIPd]
#' @export
F_pfpr_by_lm.SIPd <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute observed prevalence by RDT
#' @description Implements [F_prevalence] for the `SIPd` model.
#' @inheritParams ramp.xds::F_pfpr_by_rdt
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIPd]
#' @export
F_pfpr_by_rdt.SIPd <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute observed prevalence by PCR
#' @description Implements [F_prevalence] for the `SIPd` model.
#' @inheritParams ramp.xds::F_pfpr_by_pcr
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIPd]
#' @export
F_pfpr_by_pcr.SIPd <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Check `SIPd` (**XH**)
#'
#' @inheritParams ramp.xds::check_XH
#'
#' @return an **`xds`** model object
#' @keywords internal
#' @seealso [SIPd]
#' @export
check_XH.SIPd = function(xds_obj, i){
  return(xds_obj)
}

#' @title Get the skill set for `SIPd` (**XH**)
#'
#' @description Return the capabilities and compatibilities of the `SIPd` module as a list.
#'
#' @inheritParams ramp.xds::get_skillset_XH
#'
#' @return the **XH** module skill set, as a list
#' @keywords internal
#'
#' @seealso [SIPd]
#' @export
get_skillset_XH.SIPd = function(xds_obj,i){
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
    mda         = FALSE,
    msat        = FALSE
  )
  return(skills)
}
#' @title Compute steady states for `SIPd` (**XH**)
#' @description Compute steady states given daily FoI `foi` and fixed host density `H`.
#' @note Delegate to [steady_state_X_SIPd], which dispatches on `xds_obj$xds`.
#' @inheritParams ramp.xds::steady_state_X
#' @return a named [list] for `xde` or a named [numeric] vector for `dts`,
#' each containing `S`, `I`, and `P`
#' @export
steady_state_X.SIPd = function(foi, H, xds_obj, i=1){
  steady_state_X_SIPd(foi, H, xds_obj, i)
}

#' @rdname steady_state_X.SIPd
#' @export
steady_state_X_SIPd = function(foi, H, xds_obj, i=1){
  UseMethod("steady_state_X_SIPd", xds_obj$xds)
}

#' @title Compute steady states for `SIPd` (**XH**)
#' @description Compute steady-state infection densities for the `SIPd` module
#' given daily FoI `foi` and fixed host population density `H` for the
#' delay differential-equation system.
#' @inheritParams steady_state_X_SIPd
#' @return a named [list] with elements `S`, `I`, and `P`
#' @keywords internal
#' @seealso [SIPd]
#' @export
steady_state_X_SIPd.xde = function(foi, H, xds_obj, i=1){
  with(xds_obj$XH_obj[[i]],{
    Ieq = (foi*H*eta*(1-rho))/((foi+r+xi)*(eta+xi) +foi*(r-eta)*rho)
    Peq  = (H*xi*(foi+r+xi) + (foi*H*r*rho))/((foi+r+xi)*(eta+xi) +foi*(r-eta)*rho)
    Seq = H -Ieq - Peq
    return(list(S=as.vector(Seq), I=as.vector(Ieq), P = as.vector(Peq)))
})}

#' @title Compute steady states for `SIPd` (**XH**)
#' @description Compute steady-state infection densities for the `SIPd` module
#' given daily FoI `foi` and fixed host population density `H` for the
#' discrete-time system. Convert `foi` to `ar = exp(-foi)` for this calculation.
#' @inheritParams steady_state_X_SIPd
#' @return a named [numeric] vector containing `S`, `I`, and `P`
#' @keywords internal
#' @seealso [SIPd]
#' @export
steady_state_X_SIPd.dts = function(foi, H, xds_obj, i=1){
  ar = exp(-foi)
  with(xds_obj$XH_obj[[i]],{
  Iteq =(ar*H*eta*(1-rho))/((eta+xi)*(ar*(1-r)+(r+xi)+ar*(eta*(r-1)+r+xi)))
  Pteq = (H*(xi*(ar+r*(1-ar)+xi)+ar*(1+xi)*rho))/((eta+xi)*(ar*(1-r)+(r+xi)+ar*(eta*(r-1)+r+xi)))
  Steq = H-Iteq-Pteq
  return(c(S=Steq, I=Iteq,P =Pteq))
})}

#' Plot the density of infected individuals for the `SIPd` model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @keywords internal
#' @seealso [SIPd]
#' @export
xds_plot_X.SIPd = function(xds_obj, i=1, clrs=c("darkblue", "darkred", "darkgreen"), llty=1, add=FALSE){
  XH = xds_obj$outputs$orbits$XH[[i]]
  time = xds_obj$outputs$time

  if(add==FALSE)
    plot(time, 0*time, type = "n", ylim = c(0, max(XH$H)),
         ylab = "# Infected", xlab = "Time")

  add_lines_X_SIS(time, XH, xds_obj$nStrata[i], clrs, llty)
}


#' Add lines for the density of infected individuals for the `SIPd` model
#'
#' @param time time points for the observations
#' @param XH a list with the outputs of parse_XH_orbits.SIPd
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @keywords internal
#' @seealso [SIPd]
#' @export
add_lines_X_SIPd = function(time, XH, nStrata, clrs=c("darkblue", "darkred", "darkgreen"), llty=1){
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
