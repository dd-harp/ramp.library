
#' @title The `SIR` module for the XH component
#' @description
#' Implements the **XH** component using a Susceptible-Infectious-Recovered
#' (SIR) compartmental model of human infection dynamics.
#'
#' @section State Variables:
#' \describe{
#'   \item{`H`}{total human (or host) population density}
#'   \item{`I`}{density of infectious humans}
#'   \item{`R`}{density of recovered humans}
#' }
#' Note: susceptible density \eqn{S = H - I - R}.
#'
#' @section Parameters:
#' \describe{
#'   \item{`b`}{transmission probability from mosquito to human}
#'   \item{`c`}{transmission probability from human to mosquito}
#'   \item{`r`}{clearance rate for infections}
#'   \item{`B`}{time-dependent birth rate function \eqn{B(t, H)}}
#'   \item{`D`}{linear operator (matrix) for mortality, migration, aging, and transfers}
#' }
#'
#' @section Dynamics:
#' \deqn{
#' \begin{array}{rl}
#' dH/dt &= B(t,H) + D \cdot H \\
#' dI/dt &= hS - rI + D \cdot I \\
#' dR/dt &= rI + D \cdot R \\
#' \end{array}}
#' where \eqn{h} is the force of infection and \eqn{S = H - I - R}.
#'
#' @name SIR
#' @rdname SIR
NULL

#' @title The **XH** module skill set for `SIR`
#'
#' @description The **XH** skill set is a list of
#' a module's capabilities.
#'
#' @note This method dispatches on `class(xds_obj$XH_obj)`
#'
#' @inheritParams ramp.xds::skill_set_XH
#'
#' @return the skill set, as a list
#'
#' @keywords internal
#' @export
skill_set_XH.SIR = function(Xname = "SIR"){
  return(list(
    demography  = TRUE,
    prevalence  = TRUE,
    malaria     = TRUE,
    diagnostics = FALSE
  ))
}

#' Run checks before solving (**XH**)
#'
#' @inheritParams ramp.xds::check_XH
#'
#' @return an **`xds`** object
#' @keywords internal
#' @export
check_XH.SIR = function(xds_obj, i){
  return(xds_obj)
}


#' @title Compute the derivatives for parasite infection dynamics in human population strata
#' @description Implements [dXHdt] for the SIR model
#' @inheritParams ramp.xds::dXHdt
#' @return a [numeric] vector
#' @keywords internal
#' @export
dXHdt.SIR<- function(t, y, xds_obj, i) {

  foi <- xds_obj$terms$FoI[[i]]

  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {

      dH <- Births(t, H, births) + D_matrix %*% H
      dI <- foi*S - r*I + D_matrix %*% I
      dR <- r*I + D_matrix %*% R

      return(c(dH, dI, dR))
    })
  })
}

#' @title DTS updating for the SIS model for human / vertebrate host infections
#' @description Implements [Update_XHt] for the SIS model
#' @inheritParams ramp.xds::Update_XHt
#' @return a [numeric] vector
#' @keywords internal
#' @export
Update_XHt.SIR<- function(t, y, xds_obj, i) {

  ar <- xds_obj$terms$AR[[i]]

  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {

      Ht <- Births(t, H, births) + D_matrix %*% H
      It <- (1-r)*I + ar*S + D_matrix %*% I
      Rt <- R + r*I + D_matrix %*% R

      return(c(S=unname(St), I=unname(It), R = unname(Rt)))
    })
  })
}


#' @title Make initial values for the SIR human model, with defaults
#' @param nStrata the number of strata in the model
#' @param options a [list] to overwrite defaults
#' @param H the initial value for H
#' @param I the initial value for I
#' @param R the initial values for R
#' @return a [list]
#' @keywords internal
#' @export
make_XH_inits_SIR = function(nStrata, H, options = list(), I=1, R=0){with(options,{
  I = checkIt(I, nStrata)
  R = checkIt(R, nStrata)
  return(list(H=H, I=I, R =R))
})}





#' @title Return the parameters as a list
#' @description This method dispatches on the type of `xds_obj$XH_obj[[i]]`.
#' @inheritParams ramp.xds::change_XH_inits
#' @return an **`xds`** object
#' @keywords internal
#' @export
change_XH_inits.SIR <- function(xds_obj, i=1, options=list()) {
  with(xds_obj$XH_obj[[i]]$inits,
    with(options,{
      xds_obj$XH_obj[[i]]$S = S
      xds_obj$XH_obj[[i]]$I = I
      xds_obj$XH_obj[[i]]$R = R
      return(xds_obj)
}))}

#' @title Setup XH_obj.SIR
#' @description Implements [setup_XH_inits] for the SIR model
#' @inheritParams ramp.xds::setup_XH_inits
#' @return a [list] vector
#' @keywords internal
#' @export
setup_XH_inits.SIR = function(xds_obj, H, i, options=list()){
  xds_obj$XH_obj[[i]]$inits = make_XH_inits_SIR(xds_obj$nStrata[i], H, options)
  return(xds_obj)
}



#' @title Add indices for human population to parameter list
#' @description Implements [setup_XH_ix] for the SIR model.
#' @inheritParams ramp.xds::setup_XH_ix
#' @return none
#' @importFrom utils tail
#' @keywords internal
#' @export
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


#' @title Return the variables as a list
#' @description This method dispatches on the type of `xds_obj$XH_obj`
#' @inheritParams ramp.xds::get_XH_vars
#' @return a [list]
#' @keywords internal
#' @export
get_XH_vars.SIR <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
    H = y[H_ix]
    I = y[I_ix]
    R = y[R_ix]
    S = H-I-R
    return(list(S=S, I=I, R=R, H=H))})
}

#' @title Make parameters for SIR human model, with defaults
#' @param nStrata is the number of population strata
#' @param options a [list] that could overwrite defaults
#' @param b the proportion of infective bites that cause an infection
#' @param r the the duration of an infection
#' @param c the proportion of bites on infected humans that infect a mosquito
#' @return a [list]
#' @keywords internal
#' @export
make_XH_obj_SIR = function(nStrata, options=list(),
                          b=0.55, r=1/180, c=0.15){
  with(options,{
    XH_obj = list()
    class(XH_obj) <- c("SIR")

    XH_obj$b = checkIt(b, nStrata)
    XH_obj$c = checkIt(c, nStrata)
    XH_obj$r = checkIt(r, nStrata)

    # Ports for demographic models
    XH_obj$D_matrix = diag(0, nStrata)
    births = "zero"
    class(births) = births
    XH_obj$births = births

    return(XH_obj)
})}


#' @title Setup XH_obj.SIR
#' @description Implements [setup_XH_obj] for the SIR model
#' @inheritParams ramp.xds::setup_XH_obj
#' @return a [list] vector
#' @keywords internal
#' @export
setup_XH_obj.SIR = function(Xname, xds_obj, i, options=list()){
  XH_obj <- make_XH_obj_SIR(xds_obj$nStrata[i], options)
  class(XH_obj) <- c("SIR", paste("SIR_", xds_obj$xds, sep=""))
  xds_obj$XH_obj[[i]] = XH_obj
  return(xds_obj)
}

#' @title Return the parameters as a list
#' @description This method dispatches on the type of `xds_obj$XH_obj[[i]]`.
#' @inheritParams ramp.xds::change_XH_pars
#' @return an **`xds`** object
#' @keywords internal
#' @export
change_XH_pars.SIR <- function(xds_obj, i=1, options=list()) {
  nHabitats <- xds_obj$nHabitats
  with(xds_obj$XH_obj[[i]], with(options,{
    xds_obj$XH_obj[[i]]$b <- b
    xds_obj$XH_obj[[i]]$c <- c
    xds_obj$XH_obj[[i]]$r <- r
    return(xds_obj)
  }))}

#' @title Size of effective infectious human population
#' @description Implements [F_I] for the SIS model.
#' @inheritParams ramp.xds::F_I
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @export
F_I.SIR <- function(t, y, xds_obj, i) {
  I = y[xds_obj$XH_obj[[i]]$ix$I_ix]
  Y = with(xds_obj$XH_obj[[i]], c*I)
  return(Y)
}


#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SIR model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @export
F_H.SIR <- function(t, y, xds_obj, i){
  with(get_XH_vars(y, xds_obj, i),
    return(H))
}


#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_infectivity] for the SIR model.
#' @inheritParams ramp.xds::F_infectivity
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @export
F_infectivity.SIR <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]], b)
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_ni] for the SIR model.
#' @inheritParams ramp.xds::F_ni
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @export
F_ni.SIR <- function(vars, XH_obj) {
  with(vars, with(XH_obj, c*I/H))
}

#' @title parse the output of deSolve and return variables for the SIR model
#' @description Implements [parse_XH_orbits] for the SIR model
#' @inheritParams ramp.xds::parse_XH_orbits
#' @return none
#' @keywords internal
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
#' @description Implements [F_prevalence] for the SIR model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @keywords internal
#' @export
F_prevalence.SIR <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}


#' @title Compute the HTC for the SIR model
#' @description Implements [get_HTC] for the SIR model with demography.
#' @inheritParams ramp.xds::get_HTC
#' @return a [numeric] vector
#' @keywords internal
#' @export
get_HTC.SIR <- function(xds_obj, i) {
  with(xds_obj$XH_obj[[i]],
       HTC <- c/r,
       return(HTC)
  )
}


#' Add lines for the density of infected individuals for the SIR model
#'
#' @param time time points for the observations
#' @param XH a list with the outputs of parse_outputs_X_SIR
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#' @keywords internal
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


#' Plot the density of infected individuals for the SIR model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @keywords internal
#' @export
xds_plot_X.SIR = function(xds_obj, i=1, clrs=c("darkblue","darkred","darkgreen"), llty=1,  add=FALSE){
  XH = xds_obj$outputs$orbits$XH[[i]]
  time = xds_obj$outputs$time

  if(add==FALSE)
         plot(time, 0*time, type = "n", ylim = c(0, max(XH$H)),
              ylab = "No of. Infected", xlab = "Time")

  xds_lines_X_SIR(time, XH, xds_obj$nStrata[i], clrs, llty)
}


#' @title Compute the steady states for the  dts SEIS model as a function of the daily EIR
#' @description Compute the steady state of the  dts SIS model as a function of the daily eir.
#' @inheritParams ramp.xds::steady_state_X
#' @return the steady states as a named vector
#' @keywords internal
#' @export
steady_state_X.SIR_dts = function(foi, H, xds_obj, i=1){
  ar = exp(-foi)
  with(xds_obj$XH_obj[[i]],{
    Iteq = 0
    Steq = 0
    Rteq = H-Iteq-Steq
    return(c(S=Steq, I=Iteq, R=Rteq))
  })}

#' @title Compute the steady states for the SIR model as a function of the daily EIR
#' @description Compute the steady state of the SIR model as a function of the daily eir.
#' @inheritParams ramp.xds::steady_state_X
#' @return the steady states as a named vector
#' @keywords internal
#' @export
steady_state_X.SIR_ode = function(foi, H, xds_obj, i=1){
  with(xds_obj$XH_obj[[i]],{
    Ieq = 0
    Req = H
    Seq = 0
    return(c(S=as.vector(Seq), I=as.vector(Ieq), R=as.vector(Req)))
})}
