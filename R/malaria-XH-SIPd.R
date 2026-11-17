# specialized methods for the human SIPd model

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
#' @seealso [SIPd]
#' @export
setup_skillset_XH.SIPd = function(xds_obj,i){
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
#' @seealso [SIPd]
#' @export
check_XH.SIPd = function(xds_obj, i){
  return(xds_obj)
}

#' @title Compute derivatives for `SIPd` (**XH**)
#' @description Computes derivatives for the `SIPd` module (**XH**).
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
        treated_eta = lagderiv(t=t-eta, nr=xds_obj$ix$X[[i]]$treated_ix)
      }

      dH <- Births(t, y, xds_obj, i) - Deaths(t, H, xds_obj, i)
      dI <- (1-rho)*foi*S - (r+xi)*I - Deaths(t, I, xds_obj, i)
      dP <- rho*foi*S + xi*(S+I) - treated_eta - Deaths(t, P, xds_obj, i)
      dtreated <- (rho*foi+xi)*S + xi*I

      return(c(dH, dI, dP, dtreated))
    })
  })
}

#' @title Set up `SIPd` (**XH**)
#' @description
#'
#' This sets up an `SIPd` model object for the **XH** module for the \eqn{i^{th}} host species
#' + the object is created by [make_XH_obj_SIPd]
#' + the initial conditions are set by [setup_XH_inits]
#' + `options` is passed to both to override the defaults
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
  xds_obj = ramp.xds::ode_to_dde(xds_obj)
  xds_obj$XHname = "SIPd"
  xds_obj$XH_obj[[i]] = make_XH_obj_SIPd(xds_obj$nStrata[1], options)
  xds_obj <- setup_XH_inits(xds_obj, options=options, i=i)
  xds_obj <- setup_skillset_XH(xds_obj, i)
  xds_obj <- setup_timespent("setup", xds_obj, list(residence=residence), i)
  xds_obj <- setup_demog_obj(xds_obj, i)
  xds_obj <- setup_blood_search_weights("setup", xds_obj, i=i)
  xds_obj <- setup_time_away("setup", xds_obj, i=i)
  xds_obj <- setup_travel_eir("no_travel", xds_obj, i=i)
  return(xds_obj)
}


#' @title Make a `SIPd` **XH** model object
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
#' @export
make_XH_obj_SIPd = function(nStrata, options=list(),
                            b=0.55, r=1/180, c=0.15,
                            rho=.1, eta=25, xi=1/365){
  with(options,{
    XH_obj = list()
    class(XH_obj) <- c("SIPd")

    XH_obj$b = checkIt(b, nStrata)
    XH_obj$c = checkIt(c, nStrata)
    XH_obj$r = checkIt(r, nStrata)
    XH_obj$rho = checkIt(rho, nStrata)
    XH_obj$eta = checkIt(eta, nStrata)
    XH_obj$xi = checkIt(xi, nStrata)



    return(XH_obj)
  })}


#' @title Change parameters for `SIPd` (**XH**)
#' @description For the `SIPd` module (**XH**) for the \eqn{i^{th}} host species,
#' change parameter values named in `options` leaving other values unchanged.
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
    return(xds_obj)
  }))}


#' @title Change initial values for `SIPd` (**XH**)
#' @description For the `SIPd` module (**XH**) for the \eqn{i^{th}} host species,
#' change initial values named in `options` leaving other values unchanged.
#' @inheritParams ramp.xds::change_XH_inits
#' @inheritSection SIPd State Variables
#' @return an **`xds`** model object
#' @keywords internal
#' @seealso [SIPd]
#' @export
change_XH_inits.SIPd <- function(xds_obj, options=list(), i=1) {
  with(get_XH_inits(xds_obj, i), with(options,{
    xds_obj$XH_obj[[i]]$inits[[i]]$S = get_H(xds_obj,i)-I-P
    xds_obj$XH_obj[[i]]$inits[[i]]$I = I
    xds_obj$XH_obj[[i]]$inits[[i]]$P = P
    return(xds_obj)
  }))}


#' @title Get parameters for `SIPd` (**XH**)
#' @description Parameter values for the \eqn{i^{th}} host are
#' stored as `xds_obj$XH_obj[[i]]`. This returns the stored parameter
#' values as a list.
#' @inheritParams ramp.xds::get_XH_pars
#' @inheritSection SIPd Parameters
#' @return a [list]
#' @seealso [SIPd], [make_XH_obj_SIPd]
#' @keywords internal
#' @export
get_XH_pars.SIPd <- function(xds_obj, i=1) {
  with(xds_obj$XH_obj[[i]],list(b=b, c=c, r=r, rho=rho, eta=eta, xi=xi))
}

#' @title Derivatives for human population
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

#' @title Compute the steady states for the  dts `SIPd` model as a function of the daily EIR
#' @description Compute the steady state of the  dts `SIPd` model as a function of the daily eir.
#' @inheritParams ramp.xds::steady_state_X
#' @return the steady states as a named vector
#' @keywords internal
#' @seealso [SIPd]
#' @export
steady_state_X.SIPd = function(ar, H, XH_obj){with(XH_obj,{
  Iteq =(ar*H*eta*(1-rho))/((eta+xi)*(ar*(1-r)+(r+xi)+ar*(eta*(r-1)+r+xi)))
  Pteq = (H*(xi*(ar+r*(1-ar)+xi)+ar*(1+xi)*rho))/((eta+xi)*(ar*(1-r)+(r+xi)+ar*(eta*(r-1)+r+xi)))
  Steq = H-Iteq-Pteq
  return(c(S=Steq, I=Iteq,P =Pteq))
})}

#' @title Size of effective infectious human population
#' @description Implements [F_I] for the `SIPd` model.
#' @inheritParams ramp.xds::F_I
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIPd]
#' @export
F_I.SIPd <- function(t, y, xds_obj, i) {
  I = y[xds_obj$ix$X[[i]]$I_ix]
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



#' @title List variables for `SIPd` (**XH**)
#' @description This method dispatches on the type of `xds_obj$XH_obj`
#' @inheritParams ramp.xds::get_XH_vars
#' @return a [list]
#' @keywords internal
#' @seealso [SIPd]
#' @export
get_XH_vars.SIPd <- function(y, xds_obj, i) {
  with(xds_obj$ix$X[[i]],{
    S = y[S_ix]
    I = y[I_ix]
    P = y[P_ix]
    H = S+I+P
    return(list(S=S,I=I,P=P,H=H))})
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


#' @title Setup initial values for `SIPd` (**XH**)
#' @description Set initial values for all the variables in the `SIPd` module
#' (**XH**) for the \eqn{i^{th}} host species.
#' + [make_XH_inits_SIPd] supplies default values
#' + The values in `options` override the defaults
#' + To inspect variable names, use [get_XH_inits] or
#' see the inline module documentation.
#' @inheritParams ramp.xds::setup_XH_inits
#' @inheritSection SIPd State Variables
#' @return an **`xds`** model object
#' @keywords internal
#' @seealso [SIPd]
#' @export
setup_XH_inits.SIPd = function(xds_obj, options=list(), i=1){
  xds_obj$XH_obj[[i]]$inits = make_XH_inits_SIPd(xds_obj$nStrata[i], options)
  return(xds_obj)
}

#' @title Make initial values for `SIPd` (**XH**)
#' @inheritSection SIPd State Variables
#' @param nStrata the number of population strata
#' @param options a [list] that could override defaults
#' @param H the initial human population density
#' @param I the initial number of infectious humans
#' @param P the initial number of humans under chemoprotection
#' @return a [list]
#' @keywords internal
#' @seealso [SIPd]
#' @export
make_XH_inits_SIPd = function(nStrata,
                              options = list(),
                              H=1000, I=1, P=0){
  with(options,{
    H = checkIt(H, nStrata)
    I = checkIt(I, nStrata)
    P = checkIt(P, nStrata)
    return(list(H=H, I=I, P=P, cfoi=P*0))
})}

#' @title Parse outputs for `SIPd` (**XH**)

#' @description Implements [parse_XH_orbits] for the `SIPd` model
#'
#' @inheritParams ramp.xds::parse_XH_orbits
#' @return a named list of module variables and prevalence summaries
#' @keywords internal
#' @seealso [SIPd]
#' @export
parse_XH_orbits.SIPd <- function(outputs, xds_obj, i) {
  with(xds_obj$X_obj[[i]]$ix,{
    S = outputs[,S_ix]
    I = outputs[,I_ix]
    P = outputs[,P_ix]
    H = S+I+P
    vars <- list(S=S, I=I, P=P, H=H)
    vars$ni <- F_ni(vars, xds_obj$XH_obj[[i]])
    vars$true_pr <- F_prevalence(vars, xds_obj$XH_obj[[i]])
    return(vars)
  })}

#' @title Add indices for human population to parameter list
#' @description Implements [setup_XH_ix] for the `SIPd` model.
#' @inheritParams ramp.xds::setup_XH_ix
#' @return an **`xds`** model object
#' @importFrom utils tail
#' @keywords internal
#' @seealso [SIPd]
#' @export
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
  xds_obj$ix$X[[i]] = list(S_ix=S_ix, I_ix=I_ix, P_ix=P_ix, treated_ix=treated_ix)
  return(xds_obj)
})}


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

#' @title Compute the prevalence of infection by light microscopy
#' @description Implements [F_prevalence] for the `SIPd` model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIPd]
#' @export
F_pfpr_by_lm.SIPd <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute the prevalence of infection by RDT
#' @description Implements [F_prevalence] for the `SIPd` model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIPd]
#' @export
F_pfpr_by_rdt.SIPd <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute the prevalence of infection by PCR
#' @description Implements [F_prevalence] for the `SIPd` model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SIPd]
#' @export
F_pfpr_by_pcr.SIPd <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}

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

#' @title Compute the steady states for the `SIPd` model as a function of the daily foi
#' @description Compute the steady state of the `SIPd` model as a function of the daily eir.
#' @inheritParams ramp.xds::steady_state_X
#' @return the steady states as a named list
#' @keywords internal
#' @seealso [SIPd]
#' @export
steady_state_X.SIPd = function(foi, H, xds_obj, i=1){
  with(xds_obj$XH_obj[[i]],{
    Ieq = (foi*H*eta*(1-rho))/((foi+r+xi)*(eta+xi) +foi*(r-eta)*rho)
    Peq  = (H*xi*(foi+r+xi) + (foi*H*r*rho))/((foi+r+xi)*(eta+xi) +foi*(r-eta)*rho)
    Seq = H -Ieq - Peq
    return(list(S=as.vector(Seq), I=as.vector(Ieq), P = as.vector(Peq)))
})}
