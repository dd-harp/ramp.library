# The SEIS


#' @title The **XH** Module Skill Set
#'
#' @description The **XH** skill set is a list of
#' an module's capabilities.
#'
#' @note This method dispatches on `class(xds_obj$XH_obj)`
#'
#' @inheritParams ramp.xds::skill_set_XH
#'
#' @return the skill set, as a list
#'
#' @export
skill_set_XH.SEIS = function(Xname = "SIP"){
  return(list(
    demography  = TRUE,
    prevalence  = TRUE,
    malaria     = TRUE,
    diagnostics = FALSE,
    mda = FALSE,
    msat = FALSE
  ))
}

#' Check / update before solving
#'
#' @inheritParams  ramp.xds::check_XH
#'
#' @returns an **`xds`** model object
#' @export
check_XH.SEIS = function(xds_obj, i){
  return(xds_obj)
}

#' @title \eqn{\cal X} Component Derivatives for the `SEIS` Model
#' @description Implements [dXHdt] for the SEIS model
#' @inheritParams ramp.xds::dXHdt
#' @return a [numeric] vector
#' @export
dXHdt.SEIS <- function(t, y, xds_obj, i) {

  foi <- xds_obj$terms$FoI[[i]]

  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {
      dH <- Births(t, H, births) + D_matrix %*% H
      dE <- foi*S - nu*E +  D_matrix %*% E
      dI <- nu*E - r*I + D_matrix %*% I
      return(c(dH, dE, dI))
    })
  })
}

#' @title DTS updating for the SEIS model for human / vertebrate host infections
#' @description Implements [Update_XHt] for the SEIS model
#' @inheritParams ramp.xds::Update_XHt
#' @return a [numeric] vector
#' @export
Update_XHt.SEIS <- function(t, y, xds_obj, i) {

  ar <- xds_obj$terms$AR[[i]]

  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {

      Ht <- Births(t, H, births) +
      Et <- ar*S + (1-nr)*ar*I + nu*E + dHdt(t, E, Hpar)
      It <- nr*I + (1-nu)*E + dHdt(t, I, Hpar)

      return(c(H=unname(Ht), E=unname(Et), I=unname(It)))
    })
  })
}

#' @title Compute the steady states for the  dts SEIS model as a function of the daily EIR
#' @description Compute the steady state of the  dts SIS model as a function of the daily eir.
#' @inheritParams ramp.xds::steady_state_X
#' @return the steady states as a named vector
#' @export
steady_state_X.SEIS_dts = function(foi, H, xds_obj, i=1){
  ar = exp(-foi)
  with(xds_obj$XH_obj[[i]],{
    Iteq = (ar*H*(1-nu))/((1-nr)+ar-nu+ (nr*nu*(1-ar)))
    Eteq = (ar*H*(1-nr))/((1-nr)+ar-nu+ (nr*nu*(1-ar)))
    Steq = H-Iteq-Eteq

    return(c(S=Steq, E=Eteq, I=Iteq))
})}

#' @title Setup XH_obj.SEIS
#' @description Implements [setup_XH_obj] for the SEIS model
#' @inheritParams ramp.xds::setup_XH_obj
#' @return a [list] vector
#' @export
setup_XH_obj.SEIS = function(Xname, xds_obj, i, options=list()){
  XH_obj <- make_XH_obj_SEIS(xds_obj$nStrata[i], options)
  class(XH_obj) <- c("SEIS", paste("SEIS_", xds_obj$xds, sep=""))
  xds_obj$XH_obj[[i]] = XH_obj
  return(xds_obj)
}

#' @title Make parameters for SEIS xde human model, with defaults
#' @param nStrata is the number of population strata
#' @param options a [list] that could overwrite defaults
#' @param b transmission probability (efficiency) from mosquito to human
#' @param r recovery rate
#' @param nu 1/latent period
#' @param c transmission probability (efficiency) from human to mosquito
#' @return a [list]
#' @export
make_XH_obj_SEIS = function(nStrata, options=list(),
                            b=0.55, r=1/180, nu=1/14, c=0.15){
  with(options,{
    XH_obj = list()
    class(XH_obj) <- "SEIS"

    XH_obj$b = checkIt(b, nStrata)
    XH_obj$c = checkIt(c, nStrata)
    XH_obj$r = checkIt(r, nStrata)
    XH_obj$nu = checkIt(nu, nStrata)

    # Ports for demographic models
    XH_obj$D_matrix = diag(0, nStrata)
    births = "zero"
    class(births) = births
    XH_obj$births = births

    return(XH_obj)
  })}


#' @title Return the parameters as a list
#' @description This method dispatches on the type of `xds_obj$XH_obj[[i]]`.
#' @inheritParams ramp.xds::change_XH_pars
#' @return an **`xds`** object
#' @export
change_XH_pars.SEIS <- function(xds_obj, i=1, options=list()) {
  with(xds_obj$XH_obj[[i]], with(options,{
    xds_obj$XH_obj[[i]]$b <- b
    xds_obj$XH_obj[[i]]$c <- c
    xds_obj$XH_obj[[i]]$r <- r
    xds_obj$XH_obj[[i]]$nu <- nu
    return(xds_obj)
  }))}

#' @title Return the parameters as a list
#' @description Parameter values for the \eqn{i^{th}} host are
#' stored as `xds_obj$XH_obj[[i]]`. This returns the stored parameter
#' values as a list.
#' @inheritParams ramp.xds::get_XH_pars
#' @return a [list]
#' @seealso [make_XH_obj_SEIS]
#' @export
get_XH_pars.SEIS <- function(xds_obj, i=1) {
  with(xds_obj$XH_obj[[i]],list(b=b, c=c, r=r, nu=nu))
}

# specialized methods for the human SEIS model

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SEIS model.
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SEIS <- function(t, y, xds_obj, i) {
  I = y[xds_obj$XH_obj[[i]]$ix$I_ix]
  X = with(xds_obj$XH_obj[[i]], c*I)
  return(X)
}

#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SEIS model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SEIS <- function(t, y, xds_obj, i){
  with(get_XH_vars(y, xds_obj, i), return(H))
}


#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_infectivity] for the SEIS model.
#' @inheritParams ramp.xds::F_infectivity
#' @return a [numeric] vector of length `nStrata`
#' @export
F_infectivity.SEIS <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]], b)
}

#' @title Make initial values for the SEIS xde human model, with defaults
#' @param nStrata the number of strata in the model
#' @param H the initial human population density
#' @param options a [list] to overwrite defaults
#' @param E the initial values of the parameter E
#' @param I the initial values of the parameter I
#' @return a [list]
#' @export
make_XH_inits_SEIS = function(nStrata, H, options = list(), E=0, I=1){with(options,{
  H = checkIt(H, nStrata)
  E = checkIt(E, nStrata)
  I = checkIt(I, nStrata)
  return(list(H=H, E=E, I=I))
})}

#' @title Return the parameters as a list
#' @description This method dispatches on the type of `xds_obj$XH_obj[[i]]`.
#' @inheritParams ramp.xds::change_XH_inits
#' @return an **`xds`** object
#' @export
change_XH_inits.SEIS <- function(xds_obj, i=1, options=list()) {
  with(xds_obj$XH_obj[[i]]$inits,
       with(options,{
         xds_obj$XH_obj[[i]]$inits[[i]]$H = H
         xds_obj$XH_obj[[i]]$inits[[i]]$E = E
         xds_obj$XH_obj[[i]]$inits[[i]]$I = I
         xds_obj$XH_obj[[i]]$inits[[i]]$S = H-E-I
         return(xds_obj)
       }))}



#' @title Return the variables as a list
#' @description This method dispatches on the type of `xds_obj$XH_obj`
#' @inheritParams ramp.xds::get_XH_vars
#' @return a [list]
#' @export
get_XH_vars.SEIS <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
    H = y[H_ix]
    E = y[E_ix]
    I = y[I_ix]
    S = H-E-I
    return(list(S=S,E=E,I=I,H=H))})
}

#' @title Setup Xinits.SEIS
#' @description Implements [setup_XH_inits] for the SEIS model
#' @inheritParams ramp.xds::setup_XH_inits
#' @return a [list] vector
#' @export
setup_XH_inits.SEIS = function(xds_obj, H, i, options=list()){
  xds_obj$XH_obj[[i]]$inits = make_XH_inits_SEIS(xds_obj$nStrata[i], H, options)
  return(xds_obj)
}

#' @title Add indices for human population to parameter list
#' @description Implements [setup_XH_ix] for the SEIS model.
#' @inheritParams ramp.xds::setup_XH_ix
#' @return none
#' @importFrom utils tail
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

#' @title parse the output of deSolve and return variables for the SEIS model
#' @description Implements [parse_XH_orbits] for the SEIS model
#' @inheritParams ramp.xds::parse_XH_orbits
#' @return none
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

#' @title Compute the HTC for the SEIS model
#' @description Implements [HTC] for the SEIS model with demography.
#' @inheritParams ramp.xds::HTC
#' @return a [numeric] vector
#' @export
HTC.SEIS <- function(xds_obj, i) {
  with(xds_obj$XH_obj[[i]],
       return(c/r)
  )
}

#' @title Compute the NI
#' @description Implements [F_ni] for the SEIS model.
#' @inheritParams ramp.xds::F_ni
#' @return a [numeric] vector of length `nStrata`
#' @export
F_ni.SEIS <- function(vars, XH_obj) {
  ni = with(vars, XH_obj$c*I/H)
  return(ni)
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_prevalence] for the SEIS model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @export
F_prevalence.SEIS <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute the prevalence of infection by light microscopy
#' @description Implements [F_prevalence] for the SEIS model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pfpr_by_lm.SEIS <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute the prevalence of infection by RDT
#' @description Implements [F_prevalence] for the SEIS model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pfpr_by_rdt.SEIS <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute the prevalence of infection by pcr
#' @description Implements [F_prevalence] for the SEIS model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pfpr_by_pcr.SEIS <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}


#' Plot the density of infected individuals for the SEIS model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @export
xds_plot_X.SEIS = function(xds_obj, i=1, clrs=c("darkblue","darkred"), llty=1, add=FALSE){
  XH = xds_obj$outputs$orbits$XH[[i]]
  time = xds_obj$outputs$time

  if(add==FALSE)
    plot(time, 0*time, type = "n", ylim = c(0, max(XH$H)),
         ylab = "# Infected", xlab = "Time")

  add_lines_X_SIS(time, XH, xds_obj$nStrata[i], clrs, llty)
}

#' Add lines for the density of infected individuals for the SEIS model
#'
#' @param time time points for the observations
#' @param XH a list with the outputs of parse_XH_orbits_SEIS
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
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

#' @title Compute the steady states for the SEIS model as a function of the daily EIR
#' @description Compute the steady state of the SEIS model as a function of the daily eir.
#' @inheritParams ramp.xds::steady_state_X
#' @return the steady states as a named vector
#' @export
steady_state_X.SEIS_ode = function(foi, H, xds_obj, i=1){
  with(xds_obj$XH_obj[[i]],{
    Ieq = (foi*H*nu)/(foi*(r+nu) +r*nu)
    Eeq = (foi*H*r)/(foi*(r+nu) +r*nu)
    Seq = H- Ieq-Eeq
    return(c(S=Seq, E=Eeq, I=Ieq, H=H))
})}

