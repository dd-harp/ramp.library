# specialized methods for the human SIP model

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
skill_set_XH.SIP = function(Xname = "SIP"){
  return(list(
    demography  = TRUE,
    prevalence  = TRUE,
    malaria     = TRUE,
    diagnostics = FALSE
  ))
}

#' Check / update before solving
#'
#' @inheritParams ramp.xds::check_XH
#'
#' @returns an **`xds`** model object
#' @export
check_XH.SIP = function(xds_obj, i){
  return(xds_obj)
}

#' @title Derivatives for the `SIP` Module (**XH** Component)
#' @description Compute the derivatives for SIP compartmental model:
#' \deqn{
#' \begin{array}{rl}
#' dH/dt =& B(t,H) + D \cdot H \\
#' dI/dt =& (1-\rho) h S  - (r+\xi+\sigma ) I + D \cdot I \\
#' dP/dt =& \rho h S + \xi (S+I)  + \sigma I - \eta P + D \cdot P
#' \end{array}
#' }
#' where \eqn{S = H-I-P}; \eqn{B(t, H)} is the
#' time-dependent birth rate; and \eqn{D} the demographic matrix.
#'
#' The model sets up ports for mass treatment:
#'
#' + the parameter `xi_t` is computed as \eqn{\xi'(t) = \xi + \mbox{mda}(t)}
#'
#' + the parameter `sigma_t` is computed as \eqn{\sigma'(t) = \sigma + \mbox{msat}(t)}
#'
#' The functions `mda` and `msat` are set to zero. To configure mass treatment, see `ramp.control`
#'
#' @seealso The parameters are defined in [make_XH_obj_SIP]
#' @inheritParams ramp.xds::dXHdt
#' @return a [numeric] vector
#' @export
dXHdt.SIP <- function(t, y, xds_obj, i){

  foi <- xds_obj$terms$FoI[[i]]

  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {

      xi_t = xi + mda(t)
      sigma_t = sigma + msat(t)

      dH <- Births(t, H, births) + D_matrix %*% H
      dI <- (1-rho)*foi*S - (r+xi_t+sigma_t)*I + D_matrix %*% H
      dP <- rho*foi*S + xi_t*(S+I) + sigma_t*I - eta*P + D_matrix %*% P

      return(c(dH, dI, dP))
    })
  })
}

#' @title Setup `XH_obj` for an `SIP`
#' @description Set up the `SIP` model with parameters
#' @inheritParams ramp.xds::setup_XH_obj
#' @return a [list] vector
#' @export
setup_XH_obj.SIP = function(Xname, xds_obj, i, options=list()){
  XH_obj <- make_XH_obj_SIP(xds_obj$nStrata[i], options)
  class(XH_obj) <- c("SIP", paste("SIP_", xds_obj$xds, sep=""))
  xds_obj$XH_obj[[i]] = XH_obj
  return(xds_obj)
}


#' @title Make parameters for SIP human model, with defaults
#' @param nStrata the number of population strata
#' @param options a [list] that could overwrite defaults
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @param rho probability of successful treatment upon infection
#' @param eta prophylaxis waning rate
#' @param xi background treatment rate
#' @param sigma increased treatment rate while infected
#' @param F_mass_treat mass treatment rates as a function of time
#' @importFrom ramp.xds F_zero
#' @return a [list]
#' @export
make_XH_obj_SIP = function(nStrata, options=list(),
                           b=0.55, r=1/180, c=0.15,
                           rho =.1, eta=1/25,
                           xi =1/365, sigma = 0,
                           F_mass_treat = F_zero){
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

    # Ports for demographic models
    XH_obj$D_matrix = diag(0, nStrata)
    births = "zero"
    class(births) = births
    XH_obj$births = births

    XH_obj$mda = F_zero
    XH_obj$msat = F_zero

    return(XH_obj)
  })}



#' @title Return the parameters as a list
#' @description This method dispatches on the type of `xds_obj$XH_obj[[i]]`.
#' @inheritParams ramp.xds::change_XH_pars
#' @return an **`xds`** object
#' @export
change_XH_pars.SIP <- function(xds_obj, i=1, options=list()) {
  nHabitats <- xds_obj$nHabitats
  with(xds_obj$XH_obj[[i]], with(options,{
    xds_obj$XH_obj[[i]]$b <- b
    xds_obj$XH_obj[[i]]$c <- c
    xds_obj$XH_obj[[i]]$r <- r
    xds_obj$XH_obj[[i]]$rho <- rho
    xds_obj$XH_obj[[i]]$eta <- eta
    xds_obj$XH_obj[[i]]$xi <- xi
    xds_obj$XH_obj[[i]]$mda <- mda
    xds_obj$XH_obj[[i]]$msat <- msat
    return(xds_obj)
  }))}


#' @title Return the parameters as a list
#' @description This method dispatches on the type of `xds_obj$XH_obj[[i]]`.
#' @inheritParams ramp.xds::change_XH_inits
#' @return an **`xds`** object
#' @export
change_XH_inits.SIP <- function(xds_obj, i=1, options=list()) {
  with(get_XH_inits(xds_obj, i), with(options,{
    xds_obj$XH_obj[[i]]$inits$I = I
    xds_obj$XH_obj[[i]]$inits$P = P
    return(xds_obj)
  }))}


#' @title Return the parameters as a list
#' @description Parameter values for the \eqn{i^{th}} host are
#' stored as `xds_obj$XH_obj[[i]]`. This returns the stored parameter
#' values as a list.
#' @inheritParams ramp.xds::get_XH_pars
#' @return a [list]
#' @seealso [make_XH_obj_SIP]
#' @export
get_XH_pars.SIP <- function(xds_obj, i=1) {
  with(xds_obj$XH_obj[[i]],list(b=b, c=c, r=r, rho=rho, eta=eta, xi=xi))
}

#' @title Derivatives for human population
#' @description Implements [Update_XHt] for the SIP model.
#' @inheritParams ramp.xds::Update_XHt
#' @return a [numeric] vector
#' @export
Update_XHt.SIP <- function(t, y, xds_obj, i){

  attack <- xds_obj$AR[[i]]

  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {

      St <- (1-attack)*(S+r*I) + eta*P - xi*S
      It <- (1-r)*I + attack*(1-rho)*(S+r*I) - xi*I
      Pt <- xi*(S+I) + attack*rho*(S+r*I) + (1-eta)*P

      St <- dHdt(t, St, xds_obj$Hpar[[i]]) + Births(t, H, xds_obj$Hpar[[i]])
      It <- dHdt(t, It, xds_obj$Hpar[[i]])
      Pt <- dHdt(t, Pt, xds_obj$Hpar[[i]])

      return(c(St, It, Pt))
    })
  })
}

#' @title Compute the steady states for the  dts SIP model as a function of the daily EIR
#' @description Compute the steady state of the  dts SIP model as a function of the daily eir.
#' @inheritParams ramp.xds::steady_state_X
#' @return the steady states as a named vector
#' @export
steady_state_X.SIP_dts = function(foi, H, xds_obj, i=1){
  ar = exp(-foi)
  with(xds_obj$XH_obj[[i]],{
    Iteq =(ar*H*eta*(1-rho))/((eta+xi)*(ar*(1-r)+(r+xi)+ar*(eta*(r-1)+r+xi)))
    Pteq = (H*(xi*(ar+r*(1-ar)+xi)+ar*(1+xi)*rho))/((eta+xi)*(ar*(1-r)+(r+xi)+ar*(eta*(r-1)+r+xi)))
    Steq = H-Iteq-Pteq
    return(c(S=Steq, I=Iteq,P =Pteq))
})}

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIP model.
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SIP <- function(t, y, xds_obj, i) {
  I = y[xds_obj$XH_obj[[i]]$ix$I_ix]
  X = with(xds_obj$XH_obj[[i]], c*I)
  return(X)
}

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIP model.
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SIP <- function(t, y, xds_obj, i){
  with(get_XH_vars(y, xds_obj, i), return(H))
}


#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_infectivity] for the SIP model.
#' @inheritParams ramp.xds::F_infectivity
#' @return a [numeric] vector of length `nStrata`
#' @export
F_infectivity.SIP <- function(y, xds_obj,i) {
  with(xds_obj$XH_obj[[i]], b)
}


#' @title Return the variables as a list
#' @description This method dispatches on the type of `xds_obj$XH_obj`
#' @inheritParams ramp.xds::get_XH_vars
#' @return a [list]
#' @export
get_XH_vars.SIP <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
    H = y[H_ix]
    I = y[I_ix]
    P = y[P_ix]
    S = H-I-P
    return(list(S=S,I=I,P=P,H=H))})
}


#' @title Compute the HTC for the SIP model
#' @description Implements [HTC] for the SIP model with demography.
#' @inheritParams ramp.xds::HTC
#' @return a [numeric] vector
#' @export
HTC.SIP <- function(xds_obj, i) {
  with(xds_obj$XH_obj[[i]],
       return((1-rho)*b/(r+xi)*xi/(eta+xi))
  )
}

#' @title Setup Xinits.SIP
#' @description Implements [setup_XH_inits] for the SIP model
#' @inheritParams ramp.xds::setup_XH_inits
#' @return a [list] vector
#' @export
setup_XH_inits.SIP = function(xds_obj, H, i, options=list()){
  xds_obj$XH_obj[[i]]$inits = make_XH_inits_SIP(xds_obj$nStrata[i], H, options)
  return(xds_obj)
}

#' @title Make initial values for the SIP human model, with defaults
#' @param nStrata the number of population strata
#' @param H the initial human population density
#' @param options a [list] that could overwrite defaults
#' @param I the initial values of the parameter I
#' @param P the initial values of the parameter P
#' @return a [list]
#' @export
make_XH_inits_SIP = function(nStrata, H, options = list(),
                           I=1, P=0){with(options,{
  H = H
  I = checkIt(I, nStrata)
  P = checkIt(P, nStrata)
  return(list(H=H, I=I, P=P))
})}

#' @title Parse Solutions Matrix for the SIP model
#'
#' @title In [xds_solve], after solving a system of equations,
#' the orbits for \eqn{H}, \eqn{I} and \eqn{P} are
#' pulled from the solutions matrix, and the model
#' also computes \eqn{S=H-I-P}. The variables are returned
#' as a named list.
#'
#' @inheritParams ramp.xds::parse_XH_orbits
#'
#' @return SIP module variables as a named list
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

#' @title Add indices for human population to parameter list
#' @description Implements [setup_XH_ix] for the SIP model.
#' @inheritParams ramp.xds::setup_XH_ix
#' @return none
#' @importFrom utils tail
#' @export
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


#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_ni] for the SIP model.
#' @inheritParams ramp.xds::F_ni
#' @return a [numeric] vector of length `nStrata`
#' @export
F_ni.SIP <- function(vars, XH_obj) {
  ni = with(vars, XH_obj$c*I/H)
  return(ni)
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_prevalence] for the SIP model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @export
F_prevalence.SIP <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute the prevalence of infection by light microscopy
#' @description Implements [F_prevalence] for the SIP model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pfpr_by_lm.SIP <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute the prevalence of infection by RDT
#' @description Implements [F_prevalence] for the SIP model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pfpr_by_rdt.SIP <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute the prevalence of infection by pcr
#' @description Implements [F_prevalence] for the SIP model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pfpr_by_pcr.SIP <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}

#' Plot the density of infected individuals for the SIP model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @export
xds_plot_X.SIP = function(xds_obj, i=1, clrs=c("darkblue", "darkred", "darkgreen"), llty=1, add=FALSE){
  XH = xds_obj$outputs$orbits$XH[[i]]
  time = xds_obj$outputs$time

  if(add==FALSE)
    plot(time, 0*time, type = "n", ylim = c(0, max(XH$H)),
         ylab = "# Infected", xlab = "Time")

  add_lines_X_SIP(time, XH, xds_obj$nStrata[i], clrs, llty)
}


#' Add lines for the density of infected individuals for the SIP model
#'
#' @param time time points for the observations
#' @param XH a list with the outputs of parse_XH_orbits.SIP
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
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

#' @title Compute the steady states for the SIP model as a function of the daily foi
#' @description Compute the steady state of the SIP model as a function of the daily eir.
#' @inheritParams ramp.xds::steady_state_X
#' @return the steady states as a named vector
#' @export
steady_state_X.SIP_ode = function(foi, H, xds_obj, i=1){
  with(xds_obj$XH_obj[[i]],{
    stopifnot(sigma == 0)
    Ieq = (foi*H*eta*(1-rho))/((foi+r+xi+sigma)*(eta+xi) +foi*((r-eta)*rho+sigma)
    Peq  = (H*xi*(foi+r+xi+sigma) + (foi*H*(r*rho + sigma)))/((foi+r+xi+sigma)*(eta+xi) +foi*((r-eta)*rho+sigma)
    return(list(H=H, I=as.vector(Ieq), P = as.vector(Peq)))
})}

