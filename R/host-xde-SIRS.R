# specialized methods for the human SIRS model

#' @title The `SIRS` module for the **XH** component
#' @description
#' Implements the **XH** component using a Susceptible-Infectious-Recovered-Susceptible
#' (`SIRS`) compartmental model of human infection dynamics, extending `SIR` with
#' waning immunity so that recovered individuals return to the susceptible class.
#' This implementation includes human demographic change and support for
#' progressive stratification.
#'
#' @section State Variables:
#' \describe{
#'   \item{`H`}{total human (or host) population density}
#'   \item{`I`}{density of infectious humans}
#'   \item{`R`}{density of recovered (temporarily immune) humans}
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
#' [make_XH_obj_SIRS] supplies defaults for the parameters, which can be changed
#' using [change_XH_pars.SIRS].
#' \describe{
#'   \item{`b`}{transmission probability from mosquito to human}
#'   \item{`c`}{transmission probability from human to mosquito}
#'   \item{`r`}{clearance rate for infections}
#'   \item{`gamma`}{rate of waning immunity (loss of \eqn{R})}
#' }
#'
#' @section Infection Dynamics:
#' Susceptible individuals become infectious at rate \eqn{hS}.
#' Infections clear at the per-capita rate \eqn{r}, moving individuals into
#' the recovered class. Immunity wanes at rate \eqn{\gamma}, returning
#' recovered individuals to susceptibility.
#' \deqn{
#' \begin{array}{rl}
#' dI/dt &= hS - rI - \mathbf{D}(I) \\
#' dR/dt &= rI - \gamma R - \mathbf{D}(R) \\
#' \end{array}}
#' where \eqn{h} is the force of infection and \eqn{S = H - I - R}.
#' The mortality function uses the same model for \eqn{I} and \eqn{R} as for
#' \eqn{H}. Susceptible density is redundant and is computed after solving.
#'
#' @section Prevalence:
#' [F_prevalence.SIRS] returns \eqn{x=I/H}. This module has no diagnostic
#' sensitivity parameters or module-specific light microscopy, RDT, or PCR
#' prevalence methods.
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
#' + Basic Setup: [make_XH_obj_SIRS], [setup_XH_obj.SIRS]
#' + Dynamics: [dXHdt.SIRS]
#' + Parameters: [change_XH_pars.SIRS]
#' + Variables: [make_XH_inits_SIRS], [get_XH_inits], [change_XH_inits.SIRS], [setup_XH_inits.SIRS]
#' + [SIR], [SIS]
#'
#' @name SIRS
#' @rdname SIRS
NULL

#' @title The **XH** module skill set
#'
#' @description The **XH** skill set is a list of
#' a module's capabilities.
#'
#' @note This method dispatches on `class(xds_obj$XH_obj)`
#'
#' @inheritParams ramp.xds::setup_skillset_XH
#'
#' @return an **`xds`** model object
#' @keywords internal
#'
#' @seealso [SIRS]
#' @export
setup_skillset_XH.SIRS = function(xds_obj,i){
  skills =   list(
    demography  = TRUE,
    prevalence  = TRUE,
    malaria     = TRUE,
    diagnostics = "linear",
    mda         = FALSE,
    msat        = FALSE
  )
  xds_obj$XH_obj[[i]]$skill_set = skills
  return(xds_obj)
}

#' Run checks before solving (**XH**)
#'
#' @inheritParams ramp.xds::check_XH
#'
#' @return an **`xds`** model object
#' @keywords internal
#' @seealso [SIRS]
#' @export
check_XH.SIRS = function(xds_obj, i){
  return(xds_obj)
}


#' @title Compute derivatives for `SIRS` (**XH**)
#' @description Implements [dXHdt] for the `SIRS` model
#' @inheritParams ramp.xds::dXHdt
#' @return a [numeric] vector
#' @keywords internal
#' @inheritSection SIRS Demography
#' @inheritSection SIRS Infection Dynamics
#' @seealso [SIRS]
#' @export
dXHdt.SIRS <- function(t, y, xds_obj, i) {

  foi <- xds_obj$terms$FoI[[i]]

  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {

      dH <- Births(t, y, xds_obj, i) - Deaths(t, H, xds_obj, i)
      dI <- foi*S - r*I - Deaths(t, I, xds_obj, i)
      dR <- r*I - gamma*R - Deaths(t, R, xds_obj, i)

      return(c(dH, dI, dR))
    })
  })
}

#' @title DTS updating for `SIRS` (**XH**)
#' @description Implements [Update_XHt] for the `SIRS` model
#' @inheritParams ramp.xds::Update_XHt
#' @return a [numeric] vector
#' @keywords internal
#' @seealso [SIRS]
#' @export
Update_XHt.SIRS<- function(t, y, xds_obj, i) {

  ar <- xds_obj$terms$AR[[i]]

  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {

      Ht <- Births(t, y, xds_obj, i) + D_matrix %*% H
      It <- (1-r)*I + ar*S + D_matrix %*% I
      Rt <- (1-gamma)*R + r*I + D_matrix %*% R

      return(c(S=unname(St), I=unname(It), R = unname(Rt)))
    })
  })
}

#' @title Compute the steady states for the dts `SIRS` model as a function of the daily EIR
#' @description Compute the steady state of the  dts `SIRS` model as a function of the daily eir.
#' @inheritParams ramp.xds::steady_state_X
#' @return the steady states as a named vector
#' @keywords internal
#' @seealso [SIRS]
#' @export
steady_state_X.SIRS_dts = function(foi, H,  xds_obj, i=1){
  ar = exp(-foi)
  with(xds_obj$XH_obj[[i]],{
    Iteq = (ar*H*gamma)/(r*gamma + ar*(r+gamma))
    Rteq = (ar*H*r)/(r*gamma + ar*(r+gamma))
    Steq = H -Iteq-Rteq
  return(c(S=Steq, I=Iteq, R=Rteq))
})}


#' @title Change parameters for `SIRS` (**XH**)
#' @description For the `SIRS` module (**XH**) for the \eqn{i^{th}} host species,
#' change parameter values named in `options` leaving other values unchanged.
#' @inheritParams ramp.xds::change_XH_pars
#' @return an **`xds`** model object
#' @keywords internal
#' @inheritSection SIRS Parameters
#' @seealso [SIRS]
#' @export
change_XH_pars.SIRS <- function(xds_obj, options=list(), i=1) {
  nHabitats <- xds_obj$nHabitats
  with(xds_obj$XH_obj[[i]], with(options,{
    xds_obj$XH_obj[[i]]$b <- b
    xds_obj$XH_obj[[i]]$c <- c
    xds_obj$XH_obj[[i]]$r <- r
    xds_obj$XH_obj[[i]]$gamma <- gamma
    return(xds_obj)
  }))}


#' @title Change initial values for `SIRS` (**XH**)
#' @description For the `SIRS` module (**XH**) for the \eqn{i^{th}} host species,
#' change initial values named in `options` leaving other values unchanged.
#' @inheritParams ramp.xds::change_XH_inits
#' @return an **`xds`** model object
#' @keywords internal
#' @inheritSection SIRS State Variables
#' @seealso [SIRS]
#' @export
change_XH_inits.SIRS <- function(xds_obj, options=list(), i=1) {
  with(xds_obj$XH_obj[[i]], with(options,{
    xds_obj$XH_obj[[i]]$inits$S = S
    xds_obj$XH_obj[[i]]$inits$I = I
    xds_obj$XH_obj[[i]]$inits$R = R
    return(xds_obj)
  }))}


#' @title Compute the steady states for the `SIRS` model as a function of the daily EIR
#' @description Compute the steady state of the `SIRS` model as a function of the daily eir.
#' @inheritParams  ramp.xds::steady_state_X
#' @return the steady states as a named vector
#' @keywords internal
#' @seealso [SIRS]
#' @export
steady_state_X.SIRS = function(foi, H, xds_obj, i=1){
  with(xds_obj$XH_obj[[i]],{
    Ieq = foi*H*gamma/((r*gamma) + foi*(gamma+r))
    Req = foi*r*H/((r*gamma) +foi*(gamma+r))
    Seq = H-Ieq-Req
    return(c(S=as.vector(Seq), I=as.vector(Ieq), R = as.vector(Req)))
})}

#' @title Make initial values for `SIRS` (**XH**)
#' @param nStrata the number of strata in the model
#' @param options a [list] to override defaults
#' @param H the initial value for H
#' @param I the initial value for I
#' @param R the initial values for R
#' @return a [list]
#' @keywords internal
#' @inheritSection SIRS State Variables
#' @seealso [SIRS]
#' @export
make_XH_inits_SIRS = function(nStrata,  options = list(), H=1000, I=1, R=0){with(options,{
  H = checkIt(H, nStrata)
  I = checkIt(I, nStrata)
  R = checkIt(R, nStrata)
  return(list(H=H, I=I, R=R))
})}



#' @title Setup initial values for `SIRS` (**XH**)
#' @description Set initial values for all the variables in the `SIRS` module
#' (**XH**) for the \eqn{i^{th}} host species.
#' + [make_XH_inits_SIRS] supplies default values
#' + The values in `options` override the defaults
#' + To inspect variable names, use [get_XH_inits] or
#' see the inline module documentation.
#' @inheritParams ramp.xds::setup_XH_inits
#' @return an **`xds`** model object
#' @keywords internal
#' @inheritSection SIRS State Variables
#' @seealso [SIRS]
#' @export
setup_XH_inits.SIRS = function(xds_obj, options=list(), i=1){
  xds_obj$XH_obj[[i]]$inits = make_XH_inits_SIRS(xds_obj$nStrata[i], options)
  return(xds_obj)
}


#' @title Add indices for human population to parameter list
#' @description Implements [setup_XH_ix] for the `SIRS` model.
#' @inheritParams ramp.xds::setup_XH_ix
#' @return an **`xds`** model object
#' @importFrom utils tail
#' @keywords internal
#' @seealso [SIRS]
#' @export
setup_XH_ix.SIRS <- function(xds_obj, i) {with(xds_obj,{

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


#' @title List variables for `SIRS` (**XH**)
#' @description This method dispatches on the type of `xds_obj$XH_obj`
#' @inheritParams ramp.xds::get_XH_vars
#' @return a [list]
#' @keywords internal
#' @inheritSection SIRS State Variables
#' @seealso [SIRS]
#' @export
get_XH_vars.SIRS <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
    H = y[H_ix]
    I = y[I_ix]
    R = y[R_ix]
    S = H-I-R
    return(list(S=S, I=I, R=R, H=H))})
}

#' @title Make a `SIRS` **XH** model object
#' @param nStrata is the number of population strata
#' @param options a [list] that could override defaults
#' @param b the proportion of infective bites that cause an infection
#' @param r the clearance rate for infections
#' @param c the proportion of bites on infected humans that infect a mosquito
#' @param gamma the rate of loss of immunity
#' @return an **XH** model object
#' @keywords internal
#' @inheritSection SIRS Parameters
#' @seealso [SIRS]
#' @export
make_XH_obj_SIRS = function(nStrata, options=list(),
                          b=0.55, r=1/180, c=0.15,gamma=0.5){
  with(options,{
    XH_obj = list()
    class(XH_obj) <- c("SIRS")

    XH_obj$b = checkIt(b, nStrata)
    XH_obj$c = checkIt(c, nStrata)
    XH_obj$r = checkIt(r, nStrata)
    XH_obj$gamma = checkIt(gamma, nStrata)


    return(XH_obj)
})}



#' @title Set up `SIRS` (**XH**)
#' @description
#'
#' This sets up an `SIRS` model object for the **XH** module for the \eqn{i^{th}} host species
#' + the object is created by [make_XH_obj_SIRS]
#' + the initial conditions are set by [setup_XH_inits]
#' + `options` is passed to both to override the defaults
#'
#' @seealso [SIRS], [make_XH_obj_SIRS]
#'
#' @note The module object is stored at `xds_obj$XH_obj[[i]]`
#'
#' @inheritParams ramp.xds::setup_XH_obj
#' @return an **`xds`** model object
#' @keywords internal
#' @inheritSection SIRS State Variables
#' @inheritSection SIRS Parameters
#' @export
setup_XH_obj.SIRS = function(XHname, residence, xds_obj, options=list(), i=1){
  xds_obj$XHname = "SIRS"
  xds_obj$XH_obj[[i]] = make_XH_obj_SIRS(xds_obj$nStrata[1], options)
  xds_obj <- setup_XH_inits(xds_obj, options=options, i=i)
  xds_obj <- setup_skillset_XH(xds_obj, i)
  xds_obj <- setup_timespent("setup", xds_obj, list(residence=residence), i)
  xds_obj <- setup_demog_obj(xds_obj, i)
  xds_obj <- setup_blood_search_weights("setup", xds_obj, i=i)
  xds_obj <- setup_time_away("setup", xds_obj, i=i)
  xds_obj <- setup_travel_eir("no_travel", xds_obj, i=i)
  return(xds_obj)
}





#' @title Compute infectious density for `SIRS` (**XH**)
#' @description Implements [F_I] for the `SIRS` model.
#' @inheritParams ramp.xds::F_I
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @inheritSection SIRS Infectiousness
#' @seealso [SIRS]
#' @export
F_I.SIRS <- function(t, y, xds_obj, i) {
  I = y[xds_obj$XH_obj[[i]]$ix$I_ix]
  Y = with(xds_obj$XH_obj[[i]], c*I)
  return(Y)
}





#' @title Compute population density
#' @description Implements [F_H] for the `SIRS` model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @inheritSection SIRS Demography
#' @seealso [SIRS]
#' @export
F_H.SIRS <- function(t, y, xds_obj, i){
  with(get_XH_vars(y, xds_obj, i), return(H))
}





#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_infectivity] for the `SIRS` model.
#' @inheritParams ramp.xds::F_infectivity
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @inheritSection SIRS Dynamical Term
#' @seealso [SIRS]
#' @export
F_infectivity.SIRS <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]],return(b))
}

#' @title Compute the net infectiousness
#' @description Implements [F_ni] for the `SIRS` model.
#' @inheritParams ramp.xds::F_ni
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @inheritSection SIRS Infectiousness
#' @seealso [SIRS]
#' @export
F_ni.SIRS <- function(vars, XH_obj) {
  with(vars, with(XH_obj, c*I/H))
}

#' @title Parse outputs for `SIRS` (**XH**)
#' @description Implements [parse_XH_orbits] for the `SIRS` model
#' @inheritParams ramp.xds::parse_XH_orbits
#' @return a named list of module variables
#' @keywords internal
#' @inheritSection SIRS State Variables
#' @seealso [SIRS]
#' @export
parse_XH_orbits.SIRS <- function(outputs, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
    H = outputs[,H_ix]
    I = outputs[,I_ix]
    R = outputs[,R_ix]
    S = H-I-R
    return(list(S=S, I=I, R=R, H=H))
})}


#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_prevalence] for the `SIRS` model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @inheritSection SIRS Prevalence
#' @seealso [SIRS]
#' @export
F_prevalence.SIRS <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}





#' @title Compute the HTC for the `SIRS` model
#' @description Implements [get_HTC] for the `SIRS` model with demography.
#' @inheritParams ramp.xds::get_HTC
#' @return a [numeric] vector
#' @keywords internal
#' @seealso [SIRS]
#' @export
get_HTC.SIRS <- function(xds_obj, i) {
  with(xds_obj$XH_obj[[i]],
       HTC <- c/r,
       return(HTC)
  )
}



#' Add lines for the density of infected individuals for the `SIRS` model
#'
#' @param times time points for the observations
#' @param XH a list with the outputs of parse_outputs_X_SIRS
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#' @keywords internal
#' @seealso [SIRS]
#' @export
xds_lines_X_SIRS = function(times, XH, nStrata, clrs=c("darkblue","darkred","darkgreen"), llty=1){
  if (length(llty)< nStrata) llty = rep(llty, nStrata)
  with(XH,{
    if(nStrata==1) {
      lines(times, S, col=clrs[1], lty = llty[1])
      lines(times, I, col=clrs[2], lty = llty[1])
      lines(times, R, col=clrs[3], lty = llty[1])
    } else {
      for(i in 1:nStrata)
      lines(times, S[,i], col=clrs[1], lty = llty[i])
      lines(times, I[,i], col=clrs[2], lty = llty[i])
      lines(times, R[,i], col=clrs[3], lty = llty[i])
    }})}





#' Plot the density of infected individuals for the `SIRS` model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @keywords internal
#' @seealso [SIRS]
#' @export
xds_plot_X.SIRS = function(xds_obj, i=1, clrs=c("darkblue","darkred","darkgreen"), llty=1, add=FALSE){
  XH = xds_obj$outputs$orbits$XH[[i]]
  times = xds_obj$outputs$time

  if(add==FALSE)
    plot(times, 0*times, type = "n", ylim = c(0, max(XH$H)),
         ylab = "No of. Infected", xlab = "Time")
  xds_lines_X_SIRS(times, XH, xds_obj$nStrata[i], clrs, llty)
}
