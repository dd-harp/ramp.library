# specialized methods for the human SEIR model

#' @title The `SEIR` module for the **XH** component
#' @description
#' Implements the **XH** component using a Susceptible-Exposed-Infectious-Recovered
#' (`SEIR`) compartmental model of human infection dynamics.
#' Exposed individuals become infectious, then recover with lasting immunity.
#' This implementation includes human demographic change and support for
#' progressive stratification.
#'
#' @section State Variables:
#' \describe{
#'   \item{`H`}{total human (or host) population density}
#'   \item{`E`}{density of exposed (infected but not yet infectious) humans}
#'   \item{`I`}{density of infectious humans}
#'   \item{`R`}{density of recovered humans}
#'   \item{`S=H-E-I-R`}{density of susceptible humans or hosts}
#' }
#'
#' @inheritSection ramp.xds::SIS Demography
#'
#' @note The current `SEIR` implementation adds the result of [`Deaths`]
#' rather than subtracting it as described in the shared Demography section:
#' \eqn{dH/dt = \mathbf{B}(t,H) + \mathbf{D}(H)}.
#'
#' @section Dynamical Term:
#' \describe{
#'   \item{\eqn{h} or `foi`}{The force of infection is computed and stored by [Exposure] as `FoI`}
#' }
#' A constant fraction of infectious bites causes an infection: the parameter \eqn{b}.
#'
#' @section Parameters:
#' [make_XH_obj_SEIR] supplies defaults for the parameters, which can be changed
#' using [change_XH_pars.SEIR].
#' \describe{
#'   \item{`b`}{transmission probability from mosquito to human}
#'   \item{`c`}{transmission probability from human to mosquito}
#'   \item{`tau`}{rate of progression from exposed to infectious}
#'   \item{`r`}{clearance rate for infections}
#' }
#'
#' @section Infection Dynamics:
#' Exposed individuals become infectious at the per-capita rate \eqn{\tau}.
#' Infections clear at the per-capita rate \eqn{r}, moving individuals into
#' the recovered class with lasting immunity.
#' The current implementation uses the following equations, adding [`Deaths`]:
#' \deqn{
#' \begin{array}{rl}
#' dE/dt &= hS - \tau E + \mathbf{D}(E) \\
#' dI/dt &= \tau E - rI + \mathbf{D}(I) \\
#' dR/dt &= rI + \mathbf{D}(R) \\
#' \end{array}}
#' where \eqn{h} is the force of infection and \eqn{S = H - E - I - R}.
#' Susceptible density is redundant and is computed after solving.
#'
#' @section Prevalence:
#' [F_prevalence.SEIR] returns \eqn{x=I/H}. In contrast,
#' [parse_XH_orbits.SEIR] stores `true_pr` as \eqn{(E+I)/H}.
#' This module has no diagnostic sensitivity parameters or module-specific
#' light microscopy, RDT, or PCR prevalence methods.
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
#' + Basic Setup: [make_XH_obj_SEIR], [setup_XH_obj.SEIR]
#' + Dynamics: [dXHdt.SEIR]
#' + Parameters: [change_XH_pars.SEIR]
#' + Variables: [make_XH_inits_SEIR], [get_XH_inits], [change_XH_inits.SEIR], [setup_XH_inits.SEIR]
#' + [SIR], [SEIS]
#'
#' @name SEIR
#' @rdname SEIR
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
#' @seealso [SEIR]
#' @export
setup_skillset_XH.SEIR = function(xds_obj,i){
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
#' @seealso [SEIR]
#' @export
check_XH.SEIR = function(xds_obj, i){
  return(xds_obj)
}


#' @title Compute derivatives for `SEIR` (**XH**)
#' @description Implements [dXHdt] for the `SEIR` model
#' @inheritParams ramp.xds::dXHdt
#' @return a [numeric] vector
#' @keywords internal
#' @inheritSection SEIR Demography
#' @inheritSection SEIR Infection Dynamics
#' @seealso [SEIR]
#' @export
dXHdt.SEIR<- function(t, y, xds_obj, i) {
  # do not change this
  foi <- xds_obj$terms$FoI[[i]]

  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {

      dH <- Births(t, y, xds_obj, i) + Deaths(t, H, xds_obj, i)
      dE <- foi*S - tau*E + Deaths(t, E, xds_obj, i)
      dI <- tau*E - r*I + Deaths(t, I, xds_obj, i)
      dR <- r*I + Deaths(t, R, xds_obj, i)

       # concatenate the derivatives
      derivs = c(dH, dE, dI, dR)
      # return the derivatives
      return(derivs)
    })
  })
}


#' @title DTS updating for `SEIR` (**XH**)
#' @description Implements [Update_XHt] for the `SEIR` model
#' @inheritParams ramp.xds::Update_XHt
#' @return a [numeric] vector
#' @keywords internal
#' @seealso [SEIR]
#' @export
Update_XHt.SEIR<- function(t, y, xds_obj, i) {

  ar <- xds_obj$AR[[i]]
  Hpar <- xds_obj$Hpar[[i]]
  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {

      St <- (1-ar)*S  + dHdt(t, S, Hpar) + Births(t, y, xds_obj, i)
      Et <- a*S +(1-tau)*E + D_matrix %*% E
      It <- (1-r)*I + tau*E + D_matrix %*% I
      Rt <- R + r*I + D_matrix %*% R

      return(c(S=unname(St), E = unname(It), I=unname(It), R = unname(Rt)))
    })
  })
}

#' @title Compute the steady states for the dts `SEIR` model as a function of the daily EIR
#' @description Compute the steady state of the  dts `SEIR` model as a function of the daily eir.
#' @inheritParams ramp.xds::steady_state_X
#' @return the steady states as a named vector
#' @keywords internal
#' @seealso [SEIR]
#' @export
steady_state_X.SEIR_dts = function(foi, H, xds_obj, i=1){
  ar = exp(-foi)
  with(xds_obj$XH_obj[[i]],{
    Steq = 0
    Eteq = 0
    Iteq = 0
    Rteq = H -Steq - Eteq-Iteq
    return(c(S=Steq,  E= Eteq, I=Iteq, R=Rteq))
})}


#' @title Compute the steady states for the `SEIR` model as a function of the daily EIR
#' @description Compute the steady state of the `SEIR` model as a function of the daily eir.
#' @inheritParams  ramp.xds::steady_state_X
#' @return the steady states as a named vector
#' @keywords internal
#' @seealso [SEIR]
#' @export
steady_state_X.SEIR = function(foi, H,  xds_obj, i=1){
  with(xds_obj$XH_obj[[i]],{
    Eeq = 0
    Ieq = 0
    Req = H
    Seq = H-Ieq-Req
  return(c(S=as.vector(Seq),E = as.vector(Eeq), I=as.vector(Ieq), R = as.vector(Req)))
})}

#' @title Setup initial values for `SEIR` (**XH**)
#' @description Set initial values for all the variables in the `SEIR` module
#' (**XH**) for the \eqn{i^{th}} host species.
#' + [make_XH_inits_SEIR] supplies default values
#' + The values in `options` override the defaults
#' + To inspect variable names, use [get_XH_inits] or
#' see the inline module documentation.
#' @inheritParams ramp.xds::setup_XH_inits
#' @return an **`xds`** model object
#' @keywords internal
#' @inheritSection SEIR State Variables
#' @seealso [SEIR]
#' @export
setup_XH_inits.SEIR = function(xds_obj, options=list(), i=1){
  xds_obj$XH_obj[[i]]$inits = make_XH_inits_SEIR(xds_obj$nStrata[i], options)
  return(xds_obj)
}

#' @title Make initial values for `SEIR` (**XH**)
#' @param nStrata the number of strata in the model
#' @param options a [list] to override defaults
#' @param H the initial value for H
#' @param E the initial value for E
#' @param I the initial value for I
#' @param R the initial values for R
#' @return a [list]
#' @keywords internal
#' @inheritSection SEIR State Variables
#' @seealso [SEIR]
#' @export
make_XH_inits_SEIR = function(nStrata, options=list(), H=1000, I=1, E=0, R=1){
  with(options,{
    H = checkIt(H, nStrata)
    E = checkIt(E, nStrata)
    I = checkIt(I, nStrata)
    R = checkIt(R, nStrata)
    return(list(H=H, E=E, I=I, R=R))
})}



#' @title Change parameters for `SEIR` (**XH**)
#' @description For the `SEIR` module (**XH**) for the \eqn{i^{th}} host species,
#' change parameter values named in `options` leaving other values unchanged.
#' @inheritParams ramp.xds::change_XH_pars
#' @return an **`xds`** model object
#' @keywords internal
#' @inheritSection SEIR Parameters
#' @seealso [SEIR]
#' @export
change_XH_pars.SEIR <- function(xds_obj, options=list(), i=1) {
  nHabitats <- xds_obj$nHabitats
  with(xds_obj$XH_obj[[i]], with(options,{
    xds_obj$XH_obj[[i]]$b <- b
    xds_obj$XH_obj[[i]]$c <- c
    xds_obj$XH_obj[[i]]$r <- r
    xds_obj$XH_obj[[i]]$tau <- tau
    return(xds_obj)
  }))}


#' @title Change initial values for `SEIR` (**XH**)
#' @description For the `SEIR` module (**XH**) for the \eqn{i^{th}} host species,
#' change initial values named in `options` leaving other values unchanged.
#' @inheritParams ramp.xds::change_XH_inits
#' @return an **`xds`** model object
#' @keywords internal
#' @inheritSection SEIR State Variables
#' @seealso [SEIR]
#' @export
change_XH_inits.SEIR <- function(xds_obj, options=list(), i=1) {
  with(xds_obj$XH_obj[[i]], with(options,{
    xds_obj$XH_inits[[i]]$inits$H = H
    xds_obj$XH_inits[[i]]$inits$E = E
    xds_obj$XH_inits[[i]]$inits$I = I
    xds_obj$XH_inits[[i]]$inits$R = R
    return(xds_obj)
}))}





#' @title Add indices for human population to parameter list
#' @description Implements [setup_XH_ix] for the `SEIR` model.
#' @inheritParams ramp.xds::setup_XH_ix
#' @return an **`xds`** model object
#' @importFrom utils tail
#' @keywords internal
#' @seealso [SEIR]
#' @export
setup_XH_ix.SEIR <- function(xds_obj, i) {with(xds_obj,{

  H_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(H_ix, 1)

  E_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(E_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(I_ix, 1)

  R_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(R_ix, 1)

  xds_obj$max_ix = max_ix
  xds_obj$XH_obj[[i]]$ix = list(H_ix=H_ix,  E_ix=E_ix, I_ix=I_ix, R_ix=R_ix)
  return(xds_obj)
})}


#' @title List variables for `SEIR` (**XH**)
#' @description This method dispatches on the type of `xds_obj$XH_obj`
#' @inheritParams ramp.xds::get_XH_vars
#' @return a [list]
#' @keywords internal
#' @inheritSection SEIR State Variables
#' @seealso [SEIR]
#' @export
get_XH_vars.SEIR <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
       H = y[H_ix]
       E = y[E_ix]
       I = y[I_ix]
       R = y[R_ix]
       S = H - E - I - R
       return(list(S=S,E=E,I=I,R=R,H=H))})
}



#' @title Make a `SEIR` **XH** model object
#' @param nStrata is the number of population strata
#' @param options a [list] that could override defaults
#' @param tau  incubation rate
#' @param b the proportion of infective bites that cause an infection
#' @param r the clearance rate for infections
#' @param c the proportion of bites on infected humans that infect a mosquito
#' @return an **XH** model object
#' @keywords internal
#' @inheritSection SEIR Parameters
#' @seealso [SEIR]
#' @export
make_XH_obj_SEIR = function(nStrata, options=list(),
                          b=0.55, r=1/180, c=0.15,tau= 0.5){
  with(options,{
    XH_obj = list()
    class(XH_obj) <- c("SEIR")

    XH_obj$b = checkIt(b, nStrata)
    XH_obj$tau = checkIt(tau, nStrata)
    XH_obj$c = checkIt(c, nStrata)
    XH_obj$r = checkIt(r, nStrata)


    return(XH_obj)
  })}





#' @title Set up `SEIR` (**XH**)
#' @description
#'
#' This sets up an `SEIR` model object for the **XH** module for the \eqn{i^{th}} host species
#' + the object is created by [make_XH_obj_SEIR]
#' + the initial conditions are set by [setup_XH_inits]
#' + `options` is passed to both to override the defaults
#'
#' @seealso [SEIR], [make_XH_obj_SEIR]
#'
#' @note The module object is stored at `xds_obj$XH_obj[[i]]`
#'
#' @inheritParams ramp.xds::setup_XH_obj
#' @return an **`xds`** model object
#' @keywords internal
#' @inheritSection SEIR State Variables
#' @inheritSection SEIR Parameters
#' @export
setup_XH_obj.SEIR = function(XHname, residence, xds_obj, options=list(), i=1){
  xds_obj$XHname = "SEIR"
  xds_obj$XH_obj[[i]] = make_XH_obj_SEIR(xds_obj$nStrata[1], options)
  xds_obj <- setup_XH_inits(xds_obj, options=options, i=i)
  xds_obj <- setup_skillset_XH(xds_obj, i)
  xds_obj <- setup_timespent("setup", xds_obj, list(residence=residence), i)
  xds_obj <- setup_demog_obj(xds_obj, i)
  xds_obj <- setup_blood_search_weights("setup", xds_obj, i=i)
  xds_obj <- setup_time_away("setup", xds_obj, i=i)
  xds_obj <- setup_travel_eir("no_travel", xds_obj, i=i)
  return(xds_obj)
}


#' @title Compute infectious density for `SEIR` (**XH**)
#' @description Implements [F_I] for the `SEIR` model.
#' @inheritParams ramp.xds::F_I
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @inheritSection SEIR Infectiousness
#' @seealso [SEIR]
#' @export
F_I.SEIR <- function(t,y, xds_obj, i) {
  I = y[xds_obj$XH_obj[[i]]$ix$I_ix]
  Y = with(xds_obj$XH_obj[[i]], c*I)
  return(Y)
}


#' @title Compute population density
#' @description Implements [F_H] for the `SEIR` model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @inheritSection SEIR Demography
#' @seealso [SEIR]
#' @export
F_H.SEIR <- function(t, y, xds_obj, i){
  with(get_XH_vars(y, xds_obj, i), return(H))
}


#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_infectivity] for the `SEIR` model.
#' @inheritParams ramp.xds::F_infectivity
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @inheritSection SEIR Dynamical Term
#' @seealso [SEIR]
#' @export
F_infectivity.SEIR <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]],return(b))
}

#' @title Compute the net infectiousness
#' @description Implements [F_ni] for the `SEIR` model.
#' @inheritParams ramp.xds::F_ni
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @inheritSection SEIR Infectiousness
#' @seealso [SEIR]
#' @export
F_ni.SEIR <- function(vars, XH_obj) {
  with(vars, with(XH_obj, c*I/H))
}

#' @title Parse outputs for `SEIR` (**XH**)
#' @description Implements [parse_XH_orbits] for the `SEIR` model
#' @inheritParams ramp.xds::parse_XH_orbits
#' @return a named list of module variables
#' @keywords internal
#' @inheritSection SEIR State Variables
#' @seealso [SEIR]
#' @export
parse_XH_orbits.SEIR <- function(outputs, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
    H = outputs[,H_ix]
    E = outputs[,E_ix]
    I = outputs[,I_ix]
    R = outputs[,R_ix]
    S = H-E-I-R
    ni <- xds_obj$XH_obj[[i]]$c*I/H
    true_pr <- (I+E)/H
    return(list(S=S, E=E,I=I, R=R, H=H,ni=ni, true_pr= true_pr))
})}


#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_prevalence] for the `SEIR` model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @inheritSection SEIR Prevalence
#' @seealso [SEIR]
#' @export
F_prevalence.SEIR <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}


#' @title Compute the HTC for the `SEIR` model
#' @description Implements [get_HTC] for the `SEIR` model with demography.
#' @inheritParams ramp.xds::get_HTC
#' @return a [numeric] vector
#' @keywords internal
#' @seealso [SEIR]
#' @export
get_HTC.SEIR <- function(xds_obj, i) {
  with(xds_obj$XH_obj[[i]],
       HTC <- c/r,
       return(HTC)
  )
}


#' Add lines for the density of infected individuals for the `SEIR` model
#' @param time time points for the observations
#' @param XH a list with the outputs of parse_outputs_X_SEIR
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @keywords internal
#' @seealso [SEIR]
#' @export
xds_lines_X_SEIR = function(time, XH, nStrata, clrs=c("black","darkblue","darkred","darkgreen"), llty=1){
  if (length(llty)< nStrata) llty = rep(llty, nStrata)
  with(XH,{
    if(nStrata==1) {
      lines(time, S, col=clrs[1], lty = llty[1])
      lines(time, E, col=clrs[2], lty = llty[1])
      lines(time, I, col=clrs[2], lty = llty[1])
      lines(time, R, col=clrs[3], lty = llty[1])
    } else {
      for(i in 1:nStrata)
        lines(time, S[,i], col=clrs[1], lty = llty[i])
      lines(time, E[,i], col=clrs[2], lty = llty[i])
      lines(time, I[,i], col=clrs[2], lty = llty[i])
      lines(time, R[,i], col=clrs[3], lty = llty[i])
}})}


#' Plot the density of infected individuals for the `SEIR` model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @keywords internal
#' @seealso [SEIR]
#' @export
xds_plot_X.SEIR = function(xds_obj, i=1, clrs=c("black","darkblue","darkred","darkgreen"), llty=1,  add=FALSE){
  XH = xds_obj$outputs$orbits$XH[[i]]
  time = xds_obj$outputs$time

  if(add==FALSE)
    plot(time, 0*time, type = "n", ylim = c(0, max(XH$H)),
         ylab = "No of. Infected", xlab = "Time")
  xds_lines_X_SEIR(time, XH, xds_obj$nStrata[i], clrs, llty)
}
