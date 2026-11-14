# specialized methods for the human SEIS model

#' @title The `SEIS` module for the **XH** component
#' @description
#' `SEIS` implements the S-E-I-S family of compartmental models for malaria as a
#' system of differential equations: individuals are susceptible \eqn{(S)},
#' exposed (infected but not yet infectious) \eqn{(E)}, or infectious \eqn{(I)}.
#' After clearing infection, individuals return directly to susceptibility.
#'
#' This implementation includes:
#' + parameters to model treating and curing infections;
#' + functions to model mass treatment;
#' + functions to model human demographic change;
#' + support for progressive stratification;
#' + parameters to model detection by diagnostic methods.
#'
#' @section State Variables:
#' \describe{
#'   \item{`H`}{total human (or host) population density}
#'   \item{`E`}{density of exposed (infected but not yet infectious) humans}
#'   \item{`I`}{density of infectious humans}
#'   \item{`S=H-E-I`}{density of susceptible humans or hosts}
#' }
#'
#' @section Demography:
#' Demographic changes (**H**) are given by:
#' \deqn{dH/dt = \mathbf{B}(t,H) - \mathbf{D}(H)}
#' This model has no parameter or terms to describe disease-induced mortality.
#' \describe{
#'   \item{[Births] or \eqn{\mathbf{B}}}{the births function}
#'   \item{[Deaths] or \eqn{\mathbf{D}}}{the deaths function}
#' }
#' **Setup:**
#' + The defaults are no births and no deaths
#' + Use [setup_births] to configure [Births]
#' + Use [setup_deaths] to configure [Deaths]
#' Also, see [xds_info_demographic_change].
#'
#' @section Dynamical Term:
#' \describe{
#'   \item{\eqn{h} or `foi`}{The force of infection is computed and stored by [Exposure] as `FoI`}
#' }
#' In the `SEIS` module, there is no effect of pre-erythrocytic immunity, so
#' a constant fraction of infectious bites causes an infection: the parameter \eqn{b}.
#'
#' @section Parameters:
#' The parameters returned by [get_XH_pars.SEIS] and modifiable by
#' [change_XH_pars.SEIS] describe the dynamics or detection by diagnostic methods.
#' \describe{
#'   \item{`b`}{transmission probability from mosquito to human}
#'   \item{`c`}{transmission probability from human to mosquito}
#'   \item{`nu`}{rate of progression from exposed to infectious}
#'   \item{`r`}{clearance rate for infections}
#'   \item{`rho` or \eqn{\rho}}{the fraction of incident infections that are treated and cured}
#'   \item{`sigma` or \eqn{\sigma}}{the rate that prevalent infections are treated and cured}
#'   \item{`q_lm`}{test positivity by light microscopy}
#'   \item{`q_rdt`}{test positivity by RDT}
#'   \item{`q_pcr`}{test positivity by PCR}
#' }
#'
#' @section Infection Dynamics:
#' Exposed individuals become infectious at the per-capita rate \eqn{\nu}.
#' Infections clear naturally at the per-capita rate \eqn{r}, after which
#' individuals are susceptible to infection.
#' Incident infections are progression from exposed to infectious, \eqn{\nu E},
#' not new exposures \eqn{hS}. A fraction \eqn{\rho} of incident infections
#' are treated and cured, returning directly to susceptibility.
#' Prevalent infections are treated and cured at rate \eqn{\sigma}.
#' Mass treatment adds clearance at rate \eqn{\Xi(t) = mda(t) + msat(t)}.
#' Treatment does not change the exposed-state equation and does not add a
#' \eqn{\rho h} clearance term to prevalent infections.
#' \deqn{
#' \begin{array}{rl}
#' dE/dt &= hS - \nu E - \mathbf{D}(E) \\
#' dI/dt &= (1-\rho)\nu E - (r+\sigma+\Xi(t))I - \mathbf{D}(I) \\
#' \end{array}}
#' where \eqn{h} is the force of infection and \eqn{S = H - E - I}.
#' The mortality function uses the same model for \eqn{E} and \eqn{I} as for
#' \eqn{H}. Susceptible density is redundant and is computed after solving.
#'
#' @section Prevalence:
#' \describe{
#'   \item{\eqn{x}}{True prevalence is \eqn{x = I/H}}
#'   \item{\eqn{x_{lm}}}{Expected prevalence by light microscopy is \eqn{q_{lm}\; x}}
#'   \item{\eqn{x_{rdt}}}{Expected prevalence by RDT is \eqn{q_{rdt}\; x}}
#'   \item{\eqn{x_{pcr}}}{Expected prevalence by PCR is \eqn{q_{pcr}\; x}}
#' }
#' These prevalence functions count infectious individuals, not exposed individuals.
#'
#' @section Infectiousness:
#' The probability that a mosquito becomes infected after feeding on an
#' infectious host is \eqn{c}. Infectious density is \eqn{cI}, and net
#' infectiousness is \eqn{cI/H}.
#'
#' @section Mass Treatment:
#' Mass treatment clears infectious individuals, returning them to susceptibility;
#' exposed individuals are unaffected. At basic setup, \eqn{\Xi(t)=0}.
#' To configure mass treatment, see [setup_mass_treatment].
#' \describe{
#'   \item{`mda(t)`}{mass drug administration}
#'   \item{`msat(t)`}{mass screen and treat}
#' }
#'
#' @seealso
#' + Basic Setup: [make_XH_obj_SEIS], [setup_XH_obj.SEIS]
#' + Dynamics: [dXHdt.SEIS]
#' + Parameters: [get_XH_pars.SEIS], [change_XH_pars.SEIS]
#' + Variables: [make_XH_inits_SEIS], [get_XH_inits], [change_XH_inits.SEIS], [setup_XH_inits.SEIS]
#' + [SIS]
#'
#' @name SEIS
#' @rdname SEIS
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
#' @return the **`xds`** object
#' @keywords internal
#'
#' @seealso [SEIS]
#' @export
setup_skillset_XH.SEIS = function(xds_obj,i){
  skills =   list(
    demography  = TRUE,
    prevalence  = TRUE,
    malaria     = TRUE,
    diagnostics = "linear",
    mda         = TRUE,
    msat        = TRUE
  )
  xds_obj$XH_obj[[i]]$skill_set = skills
  return(xds_obj)
}

#' Run checks before solving (**XH**)
#'
#' @inheritParams ramp.xds::check_XH
#'
#' @return an **`xds`** object
#' @keywords internal
#' @seealso [SEIS]
#' @export
check_XH.SEIS = function(xds_obj, i){
  return(xds_obj)
}

#' @title Compute derivatives for `SEIS` (**XH**)
#' @description Computes derivatives for the `SEIS` module (**XH**).
#' @inheritSection SEIS Demography
#' @inheritSection SEIS Infection Dynamics
#' @inheritParams ramp.xds::dXHdt
#' @return a [numeric] vector
#' @keywords internal
#' @seealso [SEIS]
#' @export
dXHdt.SEIS <- function(t, y, xds_obj, i) {

  foi <- xds_obj$terms$FoI[[i]]

  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {
      dH <- Births(t, xds_obj, i) - Deaths(t, H, xds_obj, i)
      dE <- foi*S - nu*E - Deaths(t, E, xds_obj, i)
      dI <- (1-rho)*nu*E - (r+sigma)*I - Deaths(t, I, xds_obj, i)
      dI <- dI - mda(t)*I - msat(t)*I
      return(c(dH, dE, dI))
    })
  })
}

#' @title DTS updating for the `SEIS` model for human / vertebrate host infections
#' @description Implements [Update_XHt] for the `SEIS` model
#' @inheritParams ramp.xds::Update_XHt
#' @return a [numeric] vector
#' @keywords internal
#' @seealso [SEIS]
#' @export
Update_XHt.SEIS <- function(t, y, xds_obj, i) {

  ar <- xds_obj$terms$AR[[i]]

  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {

      Ht <- Births(t, xds_obj, i) +
      Et <- ar*S + (1-nr)*ar*I + nu*E + dHdt(t, E, Hpar)
      It <- nr*I + (1-nu)*E + dHdt(t, I, Hpar)

      return(c(H=unname(Ht), E=unname(Et), I=unname(It)))
    })
  })
}

#' @title Compute the steady states for the  dts `SEIS` model as a function of the daily EIR
#' @description Compute the steady state of the dts `SEIS` model as a function of the daily eir.
#' @inheritParams ramp.xds::steady_state_X
#' @return the steady states as a named vector
#' @keywords internal
#' @seealso [SEIS]
#' @export
steady_state_X.SEIS_dts = function(foi, H, xds_obj, i=1){
  ar = exp(-foi)
  with(xds_obj$XH_obj[[i]],{
    Iteq = (ar*H*(1-nu))/((1-nr)+ar-nu+ (nr*nu*(1-ar)))
    Eteq = (ar*H*(1-nr))/((1-nr)+ar-nu+ (nr*nu*(1-ar)))
    Steq = H-Iteq-Eteq

    return(c(S=Steq, E=Eteq, I=Iteq))
})}

#' @title Set up `SEIS` (**XH**)
#' @description
#'
#' This sets up an `SEIS` model object for the **XH** module for the \eqn{i^{th}} host species
#' + the object is created by [make_XH_obj_SEIS]
#' + the initial conditions are set by [setup_XH_inits]
#' + `options` is passed to both to override the defaults
#'
#' @seealso [SEIS], [make_XH_obj_SEIS]
#'
#' @note The module object is stored at `xds_obj$XH_obj[[i]]`
#' @inheritSection SEIS State Variables
#' @inheritSection SEIS Parameters
#'
#' @inheritParams ramp.xds::setup_XH_obj
#' @return an **`xds`** object
#' @keywords internal
#' @export
setup_XH_obj.SEIS = function(XHname, residence, xds_obj, options=list(), i=1){
  xds_obj$XHname = "SEIS"
  xds_obj$XH_obj[[i]] = make_XH_obj_SEIS(xds_obj$nStrata[1], options)
  xds_obj <- setup_XH_inits(xds_obj, options=options, i=i)
  xds_obj <- setup_skillset_XH(xds_obj, i)
  xds_obj <- setup_mass_treatment(xds_obj, i=i)
  xds_obj <- setup_timespent("setup", xds_obj, list(residence=residence), i)
  xds_obj <- setup_demog_obj(xds_obj, i)
  xds_obj <- setup_blood_search_weights("setup", xds_obj, i=i)
  xds_obj <- setup_time_away("setup", xds_obj, i=i)
  xds_obj <- setup_travel_eir("no_travel", xds_obj, i=i)
  return(xds_obj)
}

#' @title Make a `SEIS` **XH** model object
#' @inheritSection SEIS Parameters
#' @param nStrata is the number of population strata
#' @param options a [list] that could override defaults
#' @param b transmission probability (efficiency) from mosquito to human
#' @param r recovery rate
#' @param nu 1/latent period
#' @param c transmission probability (efficiency) from human to mosquito
#' @param q_lm test positivity by light microscopy
#' @param q_rdt test positivity by RDT
#' @param q_pcr test positivity by PCR
#' @param rho the fraction of incident infections that are treated and cured
#' @param sigma the rate prevalent infections are treated and cured
#' @return an **XH** model object
#' @keywords internal
#' @seealso [SEIS]
#' @export
make_XH_obj_SEIS = function(nStrata, options=list(),
                            b=0.55, r=1/180, nu=1/14, c=0.15,
                            q_lm=0.8, q_rdt=0.8, q_pcr=0.9,
                            rho=0, sigma=0){
  with(options,{
    XH_obj = list()
    class(XH_obj) <- "SEIS"

    XH_obj$b = checkIt(b, nStrata)
    XH_obj$c = checkIt(c, nStrata)
    XH_obj$r = checkIt(r, nStrata)
    XH_obj$nu = checkIt(nu, nStrata)
    XH_obj$rho = checkIt(rho, nStrata)
    XH_obj$sigma = checkIt(sigma, nStrata)
    XH_obj$q_lm = checkIt(q_lm, nStrata)
    XH_obj$q_rdt = checkIt(q_rdt, nStrata)
    XH_obj$q_pcr = checkIt(q_pcr, nStrata)

    return(XH_obj)
  })}


#' @title Change parameters for `SEIS` (**XH**)
#' @description For the `SEIS` module (**XH**) for the \eqn{i^{th}} host species,
#' change parameter values named in `options` leaving other values unchanged.
#' @inheritParams ramp.xds::change_XH_pars
#' @inheritSection SEIS Parameters
#' @return an **`xds`** object
#' @keywords internal
#' @seealso [SEIS]
#' @export
change_XH_pars.SEIS <- function(xds_obj, options=list(), i=1) {
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
    return(xds_obj)
  }))}

#' @title Get parameters for `SEIS` (**XH**)
#' @description Parameter values for the \eqn{i^{th}} host are
#' stored as `xds_obj$XH_obj[[i]]`. This returns the stored parameter
#' values as a list.
#' @inheritParams ramp.xds::get_XH_pars
#' @inheritSection SEIS Parameters
#' @return a [list]
#' @seealso [SEIS], [make_XH_obj_SEIS]
#' @keywords internal
#' @export
get_XH_pars.SEIS <- function(xds_obj, i=1) {
  with(xds_obj$XH_obj[[i]],list(b=b, c=c, r=r, nu=nu,
                                q_lm=q_lm, q_rdt=q_rdt, q_pcr=q_pcr,
                                rho=rho, sigma=sigma))
}

# specialized methods for the human SEIS model

#' @title Size of effective infectious human population
#' @description Implements [F_I] for the `SEIS` model.
#' @inheritParams ramp.xds::F_I
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SEIS]
#' @export
F_I.SEIS <- function(t, y, xds_obj, i) {
  I = y[xds_obj$XH_obj[[i]]$ix$I_ix]
  X = with(xds_obj$XH_obj[[i]], c*I)
  return(X)
}

#' @title Compute population density
#' @description Implements [F_H] for the `SEIS` model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SEIS]
#' @export
F_H.SEIS <- function(t, y, xds_obj, i){
  with(get_XH_vars(y, xds_obj, i), return(H))
}


#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_infectivity] for the `SEIS` model.
#' @inheritParams ramp.xds::F_infectivity
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SEIS]
#' @export
F_infectivity.SEIS <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]], b)
}

#' @title Make initial values for `SEIS` (**XH**)
#' @inheritSection SEIS State Variables
#' @param nStrata the number of strata in the model
#' @param options a [list] to override defaults
#' @param H the initial human population density
#' @param E the initial number of exposed humans
#' @param I the initial number of infectious humans
#' @return a [list]
#' @keywords internal
#' @seealso [SEIS]
#' @export
make_XH_inits_SEIS = function(nStrata, options = list(), H=1000, E=0, I=1){with(options,{
  H = checkIt(H, nStrata)
  E = checkIt(E, nStrata)
  I = checkIt(I, nStrata)
  return(list(H=H, E=E, I=I))
})}

#' @title Change initial values for `SEIS` (**XH**)
#' @description For the `SEIS` module (**XH**) for the \eqn{i^{th}} host species,
#' change initial values named in `options` leaving other values unchanged.
#' @inheritParams ramp.xds::change_XH_inits
#' @inheritSection SEIS State Variables
#' @return an **`xds`** object
#' @keywords internal
#' @seealso [SEIS]
#' @export
change_XH_inits.SEIS <- function(xds_obj, options=list(), i=1) {
  with(xds_obj$XH_obj[[i]]$inits,
       with(options,{
         xds_obj$XH_obj[[i]]$inits[[i]]$H = H
         xds_obj$XH_obj[[i]]$inits[[i]]$E = E
         xds_obj$XH_obj[[i]]$inits[[i]]$I = I
         xds_obj$XH_obj[[i]]$inits[[i]]$S = H-E-I
         return(xds_obj)
       }))}



#' @title List variables for `SEIS` (**XH**)
#' @description This method dispatches on the type of `xds_obj$XH_obj`
#' @inheritParams ramp.xds::get_XH_vars
#' @return a [list]
#' @keywords internal
#' @seealso [SEIS]
#' @export
get_XH_vars.SEIS <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
    H = y[H_ix]
    E = y[E_ix]
    I = y[I_ix]
    S = H-E-I
    return(list(S=S,E=E,I=I,H=H))})
}

#' @title Setup initial values for `SEIS` (**XH**)
#' @description Set initial values for all the variables in the `SEIS` module
#' (**XH**) for the \eqn{i^{th}} host species.
#' + [make_XH_inits_SEIS] supplies default values
#' + The values in `options` override the defaults
#' + To inspect variable names, use [get_XH_inits] or
#' see the inline module documentation.
#' @inheritParams ramp.xds::setup_XH_inits
#' @inheritSection SEIS State Variables
#' @return an **`xds`** object
#' @keywords internal
#' @seealso [SEIS]
#' @export
setup_XH_inits.SEIS = function(xds_obj, options=list(), i=1){
  xds_obj$XH_obj[[i]]$inits = make_XH_inits_SEIS(xds_obj$nStrata[i], options)
  return(xds_obj)
}

#' @title Add indices for human population to parameter list
#' @description Implements [setup_XH_ix] for the `SEIS` model.
#' @inheritParams ramp.xds::setup_XH_ix
#' @return none
#' @importFrom utils tail
#' @keywords internal
#' @seealso [SEIS]
#' @export
setup_XH_ix.SEIS <- function(xds_obj, i) {with(xds_obj,{

  H_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(H_ix, 1)

  E_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(E_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(I_ix, 1)

  xds_obj$max_ix = max_ix
  xds_obj$XH_obj[[i]]$ix = list(H_ix=H_ix, E_ix=E_ix, I_ix=I_ix)
  return(xds_obj)
})}

#' @title Parse outputs for `SEIS` (**XH**)
#' @description Implements [parse_XH_orbits] for the `SEIS` model
#' @inheritParams ramp.xds::parse_XH_orbits
#' @return none
#' @keywords internal
#' @seealso [SEIS]
#' @export
parse_XH_orbits.SEIS <- function(outputs, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
    H = outputs[,H_ix]
    E = outputs[,E_ix]
    I = outputs[,I_ix]
    S = H-E-I
    vars <- list(S=S, E=E, I=I, H=H)
    return(vars)
})}

#' @title Compute the HTC for the `SEIS` model
#' @description Implements [get_HTC] for the `SEIS` model with demography.
#' @inheritParams ramp.xds::get_HTC
#' @return a [numeric] vector
#' @keywords internal
#' @seealso [SEIS]
#' @export
get_HTC.SEIS <- function(xds_obj, i) {
  with(xds_obj$XH_obj[[i]],
       return(c/r)
  )
}

#' @title Compute the NI
#' @description Implements [F_ni] for the `SEIS` model.
#' @inheritParams ramp.xds::F_ni
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SEIS]
#' @export
F_ni.SEIS <- function(vars, XH_obj) {
  ni = with(vars, XH_obj$c*I/H)
  return(ni)
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_prevalence] for the `SEIS` model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SEIS]
#' @export
F_prevalence.SEIS <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute the prevalence of infection by light microscopy
#' @description Implements [F_prevalence] for the `SEIS` model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SEIS]
#' @export
F_pfpr_by_lm.SEIS <- function(vars, XH_obj) {
  pr = with(XH_obj, with(vars, q_lm*I/H))
  return(pr)
}

#' @title Compute the prevalence of infection by RDT
#' @description Implements [F_prevalence] for the `SEIS` model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SEIS]
#' @export
F_pfpr_by_rdt.SEIS <- function(vars, XH_obj) {
  pr = with(XH_obj, with(vars, q_rdt*I/H))
  return(pr)
}

#' @title Compute the prevalence of infection by PCR
#' @description Implements [F_prevalence] for the `SEIS` model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @seealso [SEIS]
#' @export
F_pfpr_by_pcr.SEIS <- function(vars, XH_obj) {
  pr = with(XH_obj, with(vars, q_pcr*I/H))
  return(pr)
}


#' Plot the density of infected individuals for the `SEIS` model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @keywords internal
#' @seealso [SEIS]
#' @export
xds_plot_X.SEIS = function(xds_obj, i=1, clrs=c("darkblue","darkred"), llty=1, add=FALSE){
  XH = xds_obj$outputs$orbits$XH[[i]]
  time = xds_obj$outputs$time

  if(add==FALSE)
    plot(time, 0*time, type = "n", ylim = c(0, max(XH$H)),
         ylab = "# Infected", xlab = "Time")

  add_lines_X_SIS(time, XH, xds_obj$nStrata[i], clrs, llty)
}

#' Add lines for the density of infected individuals for the `SEIS` model
#'
#' @param time time points for the observations
#' @param XH a list with the outputs of parse_XH_orbits_SEIS
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @keywords internal
#' @seealso [SEIS]
#' @export
add_lines_X_SEIS = function(time, XH, nStrata, clrs=c("darkblue","darkred"), llty=1){
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

#' @title Compute the steady states for the `SEIS` model as a function of the daily EIR
#' @description Compute the steady state of the `SEIS` model as a function of the daily eir.
#' @note Assumes no demographic change or mass treatment. Includes `rho` and `sigma`.
#' @inheritParams ramp.xds::steady_state_X
#' @return the steady states as a named vector
#' @keywords internal
#' @seealso [SEIS]
#' @export
steady_state_X.SEIS = function(foi, H, xds_obj, i=1){
  with(xds_obj$XH_obj[[i]],{
    Seq = H/(1 + foi/nu + (1-rho)*foi/(r+sigma))
    Eeq = foi*Seq/nu
    Ieq = (1-rho)*foi*Seq/(r+sigma)
    return(c(S=Seq, E=Eeq, I=Ieq, H=H))
})}
