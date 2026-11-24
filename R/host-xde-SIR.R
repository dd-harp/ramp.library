
# TODO: Review Update_XHt.SIR: it computes Ht but returns undefined St.
# Reconcile the update equations, state order, and description before documenting them.
# TODO: Review change_XH_inits.SIR repairs: write to inits and update H
# instead of S. Confirm constructor consistency, partial updates, and species isolation.

#' @title The `SIR` module for the **XH** component
#' @description
#' Implements the **XH** component using a Susceptible-Infectious-Recovered
#' (`SIR`) compartmental model of human infection dynamics.
#' Infections clear naturally, and recovered individuals have lasting immunity.
#' This implementation includes human demographic change and support for
#' progressive stratification.
#'
#' @section State Variables:
#' \describe{
#'   \item{`H`}{total human (or host) population density}
#'   \item{`I`}{density of infectious humans}
#'   \item{`R`}{density of recovered humans}
#'   \item{`S=H-I-R`}{density of susceptible humans or hosts}
#' }
#'
#' @inheritSection ramp.xds::SIS Demography
#'
#' @section Dynamical Term:
#' \describe{
#'   \item{\eqn{h} or `foi`}{The force of infection is computed and stored by [Exposure] as `FoI`}
#' }
#' A constant fraction of infectious bites causes an infection: the parameter \eqn{b}.
#'
#' @section Parameters:
#' [make_XH_obj_SIR] supplies defaults for the parameters, which can be changed
#' using [change_XH_pars.SIR].
#' \describe{
#'   \item{`q`}{test positivity (default `0.8`)}
#'   \item{`b`}{transmission probability from mosquito to human}
#'   \item{`c`}{transmission probability from human to mosquito}
#'   \item{`r`}{clearance rate for infections}
#' }
#'
#' @section Infection Dynamics:
#' Susceptible individuals become infectious at rate \eqn{hS}.
#' Infections clear at the per-capita rate \eqn{r}, moving individuals into
#' the recovered class. Recovered individuals do not lose immunity.
#' \deqn{
#' \begin{array}{rl}
#' dI/dt &= hS - rI - \mathbf{D}(I) \\
#' dR/dt &= rI - \mathbf{D}(R) \\
#' \end{array}}
#' where \eqn{h} is the force of infection and \eqn{S = H - I - R}.
#' The mortality function uses the same model for \eqn{I} and \eqn{R} as for
#' \eqn{H}. Susceptible density is redundant and is computed after solving.
#'
#' @section Prevalence:
#' [F_pr_obs.SIR] returns observed prevalence \eqn{x_q=qI/H}.
#'
#' [F_prevalence.SIR] returns \eqn{x=I/H}. The parsed `true_pr` output
#' also uses \eqn{I/H}. This module has no diagnostic sensitivity parameters
#' or module-specific light microscopy, RDT, or PCR prevalence methods.
#'
#' @section Infectiousness:
#' The probability that a mosquito becomes infected after feeding on an
#' infectious host is \eqn{c}. Infectious density is \eqn{cI}, and net
#' infectiousness is \eqn{cI/H}.
#'
#' @section Mass Treatment:
#' This implementation does not initialize or use `mda(t)` or `msat(t)`,
#' and has no treatment parameters `rho` or `sigma`.
#'
#' @seealso
#' + Basic Setup: [make_XH_obj_SIR], [setup_XH_obj.SIR]
#' + Dynamics: [dXHdt.SIR]
#' + Parameters: [change_XH_pars.SIR]
#' + Variables: [make_XH_inits_SIR], [get_XH_inits], [change_XH_inits.SIR], [setup_XH_inits.SIR]
#' + [SIS]
#'
#' @name SIR
#' @rdname SIR
NULL


#' @title Compute the derivatives for `SIR` (**XH**)
#' @description The derivatives for the `SIR` module are:
#' \deqn{\begin{array}{rl}
#' dH/dt &= \mathbf{B}(t,y)-\mathbf{D}(t,H) \\
#' dI/dt &= hS-rI-\mathbf{D}(t,I) \\
#' dR/dt &= rI-\mathbf{D}(t,R)
#' \end{array}}
#' Here \eqn{h} is `foi`, \eqn{S=H-I-R},
#' \eqn{\mathbf{B}(t,y)} is returned by [ramp.xds::Births], and
#' \eqn{\mathbf{D}(t,X)} is returned by [ramp.xds::Deaths] for population \eqn{X}.
#' @inheritParams ramp.xds::dXHdt
#' @return a [numeric] vector
#' @keywords internal
#' @inheritSection SIR Demography
#' @inheritSection SIR Infection Dynamics
#' @seealso [SIR]
#' @export
dXHdt.SIR<- function(t, y, xds_obj, i) {

  foi <- xds_obj$terms$FoI[[i]]

  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {

      dH <- Births(t, y, xds_obj, i) - Deaths(t, H, xds_obj, i)
      dI <- foi*S - r*I - Deaths(t, I, xds_obj, i)
      dR <- r*I - Deaths(t, R, xds_obj, i)

      return(c(dH, dI, dR))
    })
  })
}

#' @title Update states for `SIR` (**XH**)
#' @description Implements [Update_XHt] for the `SIR` model
#' @inheritParams ramp.xds::Update_XHt
#' @return a [numeric] vector
#' @keywords internal
#' @seealso [SIR]
#' @export
Update_XHt.SIR<- function(t, y, xds_obj, i) {

  ar <- xds_obj$terms$AR[[i]]

  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {

      Ht <- Births(t, y, xds_obj, i) + D_matrix %*% H
      It <- (1-r)*I + ar*S + D_matrix %*% I
      Rt <- R + r*I + D_matrix %*% R

      return(c(S=unname(St), I=unname(It), R = unname(Rt)))
    })
  })
}


#' @title Get the state variables for `SIR` (**XH**)
#' @description Extract the state variables for the `SIR` module
#' from `y` and return them by name in a list for use during model solving.
#' Also compute susceptible density `S = H - I - R`.
#' @inheritParams ramp.xds::get_XH_vars
#' @return a named [list] of state variables and module-specific derived quantities
#' @keywords internal
#' @inheritSection SIR State Variables
#' @seealso [SIR]
#' @export
#' @note Variable indices are stored at `xds_obj$XH_obj[[i]]$ix`.
get_XH_vars.SIR <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
    H = y[H_ix]
    I = y[I_ix]
    R = y[R_ix]
    S = H-I-R
    return(list(S=S, I=I, R=R, H=H))})
}

#' @title Compute infectious density for `SIR` (**XH**)
#' @description Implements [F_I] for the `SIR` model.
#' @inheritParams ramp.xds::F_I
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @inheritSection SIR Infectiousness
#' @seealso [SIR]
#' @export
F_I.SIR <- function(t, y, xds_obj, i) {
  I = y[xds_obj$XH_obj[[i]]$ix$I_ix]
  Y = with(xds_obj$XH_obj[[i]], c*I)
  return(Y)
}


#' @title Compute population density
#' @description Implements [F_H] for the `SIR` model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @inheritSection SIR Demography
#' @seealso [SIR]
#' @export
F_H.SIR <- function(t, y, xds_obj, i){
  with(get_XH_vars(y, xds_obj, i),
    return(H))
}


#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_infectivity] for the `SIR` model.
#' @inheritParams ramp.xds::F_infectivity
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @inheritSection SIR Dynamical Term
#' @seealso [SIR]
#' @export
F_infectivity.SIR <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]], b)
}

#' @title Make a module object for `SIR` (**XH**)
#' @description Create an `SIR` module object and initialize its parameters.
#' Inspect parameters with [get_XH_pars.SIR] and change supported values
#' with [change_XH_pars.SIR]. Values in `options` override parameter arguments.
#' @seealso [get_XH_pars.SIR], [change_XH_pars.SIR]
#' @param nStrata is the number of population strata
#' @param options a [list] that could override defaults
#' @param b the proportion of infective bites that cause an infection
#' @param r the clearance rate for infections
#' @param c the proportion of bites on infected humans that infect a mosquito
#' @return an **XH** model object
#' @keywords internal
#' @inheritSection SIR Parameters
#' @seealso [SIR]
#' @param q test positivity; a numeric scalar or input of length `nStrata`
#' @export
make_XH_obj_SIR = function(nStrata, options=list(),
                          b=0.55, r=1/180, c=0.15, q=0.8){
  with(options,{
    XH_obj = list()
    class(XH_obj) <- c("SIR")

    XH_obj$b = checkIt(b, nStrata)
    XH_obj$c = checkIt(c, nStrata)
    XH_obj$r = checkIt(r, nStrata)


    XH_obj$q = checkIt(q, nStrata)

    return(XH_obj)
})}

#' @title Set up the module object for `SIR` (**XH**)
#' @description
#'
#' Set up an `SIR` module object for the **XH** component and the \eqn{i^{th}} host species.
#'
#' + Set the module name to `SIR`.
#' + Create the module object with [make_XH_obj_SIR], passing `nStrata[i]` and `options`.
#' + Set initial values using [ramp.xds::setup_XH_inits], passing `options` to override defaults.
#' + Initialize the time-spent matrix using [ramp.xds::setup_timespent.setup].
#' + Initialize birth and death rates using [ramp.xds::setup_demog_obj].
#' + Initialize blood search weights using [ramp.xds::setup_blood_search_weights.setup].
#' + Initialize time away using [ramp.xds::setup_time_away.setup].
#' + Initialize travel EIR using [ramp.xds::setup_travel_eir.no_travel].
#' + Return the updated **`xds`** model object.
#'
#' @seealso [SIR], [make_XH_obj_SIR]
#'
#' @note The module object is stored at `xds_obj$XH_obj[[i]]`
#'
#' @inheritParams ramp.xds::setup_XH_obj
#' @return an **`xds`** model object
#' @keywords internal
#' @inheritSection SIR State Variables
#' @inheritSection SIR Parameters
#' @export
setup_XH_obj.SIR = function(XHname, residence, xds_obj, options=list(), i=1){
  xds_obj$XHname = "SIR"
  xds_obj$XH_obj[[i]] = make_XH_obj_SIR(xds_obj$nStrata[i], options)
  xds_obj <- setup_XH_inits(xds_obj, options=options, i=i)
  xds_obj <- setup_timespent("setup", xds_obj, list(residence=residence), i)
  xds_obj <- setup_demog_obj(xds_obj, i)
  xds_obj <- setup_blood_search_weights("setup", xds_obj, i=i)
  xds_obj <- setup_time_away("setup", xds_obj, i=i)
  xds_obj <- setup_travel_eir("no_travel", xds_obj, i=i)
  return(xds_obj)
}


#' @title Get parameters for `SIR` (**XH**)
#' @description Return the parameters for the `SIR` module as a named list.
#' @inheritParams ramp.xds::get_XH_pars
#' @inheritSection SIR Parameters
#' @seealso [make_XH_obj_SIR], [change_XH_pars.SIR]
#' @return a named [list] of parameter values
#' @keywords internal
#' @export
get_XH_pars.SIR <- function(xds_obj, i=1) {
  with(xds_obj$XH_obj[[i]], list(q=q, b=b, c=c, r=r))
}

#' @title Change parameters for `SIR` (**XH**)
#' @seealso [make_XH_obj_SIR], [get_XH_pars.SIR]
#' @description Change parameter values named in `options` for the `SIR`
#' module, leaving other values unchanged.
#' @inheritParams ramp.xds::change_XH_pars
#' @return an **`xds`** model object
#' @keywords internal
#' @inheritSection SIR Parameters
#' @seealso [SIR]
#' @export
change_XH_pars.SIR <- function(xds_obj, options=list(), i=1) {
  nHabitats <- xds_obj$nHabitats
  with(xds_obj$XH_obj[[i]], with(options,{
    xds_obj$XH_obj[[i]]$b <- b
    xds_obj$XH_obj[[i]]$c <- c
    xds_obj$XH_obj[[i]]$r <- r
    xds_obj$XH_obj[[i]]$q <- q
    return(xds_obj)
  }))}


#' @title Make initial values for `SIR` (**XH**)
#' @description Create a named list of initial values for the `SIR` module.
#' Values in `options` override the defaults.
#' @details Initial values must be numeric scalars or have length `nStrata`.
#' Scalars are expanded to that length.
#' @param nStrata the number of strata in the model
#' @param options a [list] to override defaults
#' @param H the initial value for H
#' @param I the initial value for I
#' @param R the initial values for R
#' @return a named [list] of initial values
#' @keywords internal
#' @inheritSection SIR State Variables
#' @seealso [SIR]
#' @export
#' @seealso [setup_XH_inits.SIR], [change_XH_inits.SIR]
make_XH_inits_SIR = function(nStrata, options = list(), H=1000, I=1, R=0){with(options,{
  H = checkIt(H, nStrata)
  I = checkIt(I, nStrata)
  R = checkIt(R, nStrata)
  return(list(H=H, I=I, R =R))
})}

#' @title Set up initial values for `SIR` (**XH**)
#' @description Set up initial values for the `SIR` module using
#' [make_XH_inits_SIR]. Values in `options` override the defaults.
#' @inheritParams ramp.xds::setup_XH_inits
#' @return an **`xds`** model object
#' @keywords internal
#' @inheritSection SIR State Variables
#' @seealso [SIR]
#' @export
#' @seealso [change_XH_inits.SIR]
#' @note Initial values are stored at `xds_obj$XH_obj[[i]]$inits`.
#' @seealso [ramp.xds::get_XH_inits]
setup_XH_inits.SIR = function(xds_obj, options=list(), i=1){
  xds_obj$XH_obj[[i]]$inits = make_XH_inits_SIR(xds_obj$nStrata[i], options)
  return(xds_obj)
}





#' @title Change initial values for `SIR` (**XH**)
#' @description Change initial values named in `options` for the `SIR`
#' module, leaving other values unchanged.
#' @inheritParams ramp.xds::change_XH_inits
#' @return an **`xds`** model object
#' @keywords internal
#' @inheritSection SIR State Variables
#' @seealso [SIR]
#' @export
#' @seealso [make_XH_inits_SIR], [setup_XH_inits.SIR]
#' @note Stored initial values for this module are at
#' `xds_obj$XH_obj[[i]]$inits`.
change_XH_inits.SIR <- function(xds_obj, options=list(), i=1) {
  with(xds_obj$XH_obj[[i]]$inits,
    with(options,{
      xds_obj$XH_obj[[i]]$inits$H = H
      xds_obj$XH_obj[[i]]$inits$I = I
      xds_obj$XH_obj[[i]]$inits$R = R
      return(xds_obj)
}))}



#' @title Set up variable indices for `SIR` (**XH**)
#' @description Assign consecutive variable indices for the `SIR` module,
#' starting after `xds_obj$max_ix`, and update the maximum index.
#' Allocate blocks in the order `H_ix`, `I_ix`, `R_ix`.
#' Each block has length `nStrata[i]`.
#' @inheritParams ramp.xds::setup_XH_ix
#' @return an **`xds`** model object
#' @importFrom utils tail
#' @keywords internal
#' @seealso [SIR]
#' @export
#' @note Indices are stored at `xds_obj$XH_obj[[i]]$ix`.
#' The maximum index is stored at `xds_obj$max_ix`.
setup_XH_ix.SIR <- function(xds_obj, i) {with(xds_obj,{

  H_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(H_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(I_ix, 1)

  R_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(R_ix, 1)

  xds_obj$max_ix = max_ix
  xds_obj$XH_obj[[i]]$ix = list(H_ix=H_ix, I_ix=I_ix, R_ix=R_ix)
  return(xds_obj)
})}

#' @title Parse outputs for `SIR` (**XH**)
#' @description Extract the model state outputs for the `SIR` module
#' from `outputs` and return them by name in a list.
#' Also compute `S = H-I-R`, `ni = c*I/H`, and `true_pr = I/H`.
#' @note Variable indices are stored at `xds_obj$XH_obj[[i]]$ix`.
#' @inheritParams ramp.xds::parse_XH_orbits
#' @return a [list]
#' @keywords internal
#' @inheritSection SIR State Variables
#' @seealso [SIR]
#' @export
parse_XH_orbits.SIR <- function(outputs, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
    H <- outputs[,H_ix]
    I <- outputs[,I_ix]
    R <- outputs[,R_ix]
    S <- H-I-R
    ni <- xds_obj$XH_obj[[i]]$c*I/H
    true_pr <- I/H
    vars <- list(S=S, I=I, R=R, H=H, ni=ni, true_pr=true_pr)
    return(vars)
})}


#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_prevalence] for the `SIR` model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @inheritSection SIR Prevalence
#' @seealso [SIR]
#' @export
F_prevalence.SIR <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute observed prevalence
#' @description Return observed prevalence for the `SIR` module as
#' `q` times the true prevalence returned by [ramp.xds::F_prevalence].
#' @inheritParams ramp.xds::F_pr_obs
#' @return a [numeric] vector of length `nStrata`
#' @seealso [SIR]
#' @keywords internal
#' @export
F_pr_obs.SIR <- function(vars, XH_obj) {
  XH_obj$q * ramp.xds::F_prevalence(vars, XH_obj)
}

#' @title Compute the net infectiousness
#' @description Implements [F_ni] for the `SIR` model.
#' @inheritParams ramp.xds::F_ni
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @inheritSection SIR Infectiousness
#' @seealso [SIR]
#' @export
F_ni.SIR <- function(vars, XH_obj) {
  with(vars, with(XH_obj, c*I/H))
}


#' @title Compute the HTC for the `SIR` model
#' @description Implements [get_HTC] for the `SIR` model with demography.
#' @inheritParams ramp.xds::get_HTC
#' @return a [numeric] vector
#' @keywords internal
#' @seealso [SIR]
#' @export
get_HTC.SIR <- function(xds_obj, i) {
  with(xds_obj$XH_obj[[i]],
       HTC <- c/r,
       return(HTC)
  )
}

#' @title Check `SIR` (**XH**)
#'
#' @inheritParams ramp.xds::check_XH
#'
#' @return an **`xds`** model object
#' @keywords internal
#' @seealso [SIR]
#' @export
check_XH.SIR = function(xds_obj, i){
  return(xds_obj)
}

#' @title Get the skill set for `SIR` (**XH**)
#'
#' @description Return the capabilities and compatibilities of the `SIR` module as a list.
#'
#' @inheritParams ramp.xds::get_skillset_XH
#'
#' @return the **XH** module skill set, as a list
#' @keywords internal
#'
#' @seealso [SIR]
#' @export
get_skillset_XH.SIR = function(xds_obj,i){
  skills =   list(
    dts         = FALSE,
    xde         = TRUE,
    delay       = FALSE,
    demography  = TRUE,
    malaria     = FALSE,
    pr_obs      = TRUE,
    mda         = FALSE,
    msat        = FALSE
  )
  return(skills)
}
#' @title Compute steady states for `SIR` (**XH**)
#' @description Compute steady states given daily FoI `foi` and fixed host density `H`.
#' @note Delegate to [steady_state_X_SIR], which dispatches on `xds_obj$xds`.
#' @inheritParams ramp.xds::steady_state_X
#' @return a named [numeric] vector containing `S`, `I`, and `R`
#' @export
steady_state_X.SIR = function(foi, H, xds_obj, i=1){
  steady_state_X_SIR(foi, H, xds_obj, i)
}

#' @rdname steady_state_X.SIR
#' @export
steady_state_X_SIR = function(foi, H, xds_obj, i=1){
  UseMethod("steady_state_X_SIR", xds_obj$xds)
}

#' @title Compute steady states for `SIR` (**XH**)
#' @description Compute steady-state infection densities for the `SIR` module
#' given daily FoI `foi` and fixed host population density `H` for the
#' differential-equation system.
#' The current calculation sets `S` and `I` to zero and `R` to `H`.
#' @inheritParams steady_state_X_SIR
#' @return a named [numeric] vector containing `S`, `I`, and `R`
#' @keywords internal
#' @seealso [SIR]
#' @export
steady_state_X_SIR.xde = function(foi, H, xds_obj, i=1){
  with(xds_obj$XH_obj[[i]],{
    Ieq = 0
    Req = H
    Seq = 0
    return(c(S=as.vector(Seq), I=as.vector(Ieq), R=as.vector(Req)))
})}


#' @title Compute steady states for `SIR` (**XH**)
#' @description Compute steady-state infection densities for the `SIR` module
#' given daily FoI `foi` and fixed host population density `H` for the
#' discrete-time system.
#' The current calculation sets `S` and `I` to zero and `R` to `H`.
#' @inheritParams steady_state_X_SIR
#' @return a named [numeric] vector containing `S`, `I`, and `R`
#' @keywords internal
#' @seealso [SIR]
#' @export
steady_state_X_SIR.dts = function(foi, H, xds_obj, i=1){
  ar = exp(-foi)
  with(xds_obj$XH_obj[[i]],{
    Iteq = 0
    Steq = 0
    Rteq = H-Iteq-Steq
    return(c(S=Steq, I=Iteq, R=Rteq))
  })}


#' Plot the density of infected individuals for the `SIR` model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @keywords internal
#' @seealso [SIR]
#' @export
xds_plot_X.SIR = function(xds_obj, i=1, clrs=c("darkblue","darkred","darkgreen"), llty=1,  add=FALSE){
  XH = xds_obj$outputs$orbits$XH[[i]]
  time = xds_obj$outputs$time

  if(add==FALSE)
         plot(time, 0*time, type = "n", ylim = c(0, max(XH$H)),
              ylab = "No of. Infected", xlab = "Time")

  xds_lines_X_SIR(time, XH, xds_obj$nStrata[i], clrs, llty)
}


#' Add lines for the density of infected individuals for the `SIR` model
#'
#' @param time time points for the observations
#' @param XH a list with the outputs of parse_outputs_X_SIR
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#' @keywords internal
#' @seealso [SIR]
#' @export
xds_lines_X_SIR = function(time, XH, nStrata, clrs=c("darkblue","darkred","darkgreen"), llty=1){
  if (length(llty)< nStrata) llty = rep(llty, nStrata)
  with(XH,{
    if(nStrata==1) {
      lines(time, S, col=clrs[1], lty = llty[1])
      lines(time, I, col=clrs[2], lty = llty[1])
      lines(time, R, col=clrs[3], lty = llty[1])
    } else {
      for(i in 1:nStrata)
      lines(time, S[,i], col=clrs[1], lty = llty[i])
      lines(time, I[,i], col=clrs[2], lty = llty[i])
      lines(time, R[,i], col=clrs[3], lty = llty[i])
}})}
