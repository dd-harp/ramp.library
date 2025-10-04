#' @title The **XH** Module Skill Set
#'
#' @description The **XH** skill set is a list of
#' an module's capabilities.
#'
#' @note This method dispatches on `class(xds_obj$XH_obj)`
#'
#' @inheritParams ramp.xds::skill_set_XH
#'
#' @return the `SEIR` *XH* module skill set, a list
#'
#' @export
skill_set_XH.SEIR = function(Xname = "SEIR"){
  return(list(
    H_dynamics = TRUE,
    mda        = FALSE,
    msat       = FALSE,
    malaria    = FALSE,
    pr_obs     = TRUE,
    pf_lm      = FALSE,
    pf_rdt     = FALSE,
    pf_pcr     = FALSE
  ))
}

#' Check / update before solving
#'
#' @inheritParams ramp.xds::check_XH
#'
#' @returns an **`xds`** model object
#' @export
check_XH.SEIR = function(xds_obj, i){
  return(xds_obj)
}


#' @title Compute the derivatives for parasite infection dynamics in human population strata
#' @description Implements [dXHdt] for the SEIR model
#' @inheritParams ramp.xds::dXHdt
#' @return a [numeric] vector
#' @export
dXHdt.SEIR<- function(t, y, xds_obj, i) {
  # do not change this
  foi <- xds_obj$terms$FoI[[i]]

  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {

      dH <- Births(t, H, births) + D_matrix %*% H
      dE <- foi*S - tau*E + D_matrix %*% E
      dI <- tau*E - r*I + D_matrix %*% I
      dR <- r*I + D_matrix %*% R

       # concatenate the derivatives
      derivs = c(dH, dE, dI, dR)
      # return the derivatives
      return(derivs)
    })
  })
}


#' @title DTS updating for the SIS model for human / vertebrate host infections
#' @description Implements [Update_XHt] for the SIS model
#' @inheritParams ramp.xds::Update_XHt
#' @return a [numeric] vector
#' @export
Update_XHt.SEIR<- function(t, y, xds_obj, i) {

  ar <- xds_obj$AR[[i]]
  Hpar <- xds_obj$Hpar[[i]]
  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]], {

      St <- (1-ar)*S  + dHdt(t, S, Hpar) + Births(t, H, Hpar)
      Et <- a*S +(1-tau)*E + D_matrix %*% E
      It <- (1-r)*I + tau*E + D_matrix %*% I
      Rt <- R + r*I + D_matrix %*% R

      return(c(S=unname(St), E = unname(It), I=unname(It), R = unname(Rt)))
    })
  })
}

#' @title Compute the steady states for the  dts SEIS model as a function of the daily EIR
#' @description Compute the steady state of the  dts SIS model as a function of the daily eir.
#' @inheritParams ramp.xds::steady_state_X
#' @return the steady states as a named vector
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


#' @title Compute the steady states for the SIRS model as a function of the daily EIR
#' @description Compute the steady state of the SIRS model as a function of the daily eir.
#' @inheritParams  ramp.xds::steady_state_X
#' @return the steady states as a named vector
#' @export
steady_state_X.SEIR_ode = function(foi, H,  xds_obj, i=1){
  with(xds_obj$XH_obj[[i]],{
    Eeq = 0
    Ieq = 0
    Req = H
    Seq = H-Ieq-Req
  return(c(S=as.vector(Seq),E = as.vector(Eeq), I=as.vector(Ieq), R = as.vector(Req)))
})}

#' @title Setup XH_obj.SEIR
#' @description Implements [setup_XH_inits] for the SEIR model
#' @inheritParams ramp.xds::setup_XH_inits
#' @return a [list] vector
#' @export
setup_XH_inits.SEIR = function(xds_obj, H, i, options=list()){
  xds_obj$XH_obj[[i]]$inits = make_XH_inits_SEIR(xds_obj$nStrata[i], H, options)
  return(xds_obj)
}

#' @title Make initial values for the SEIR human model, with defaults
#' @param nStrata the number of strata in the model
#' @param H the initial value for H
#' @param options a [list] to overwrite defaults
#' @param E the initial value for E
#' @param I the initial value for I
#' @param R the initial values for R
#' @return a [list]
#' @export
make_XH_inits_SEIR = function(nStrata, H, options=list(), I=1, E=0, R=1){
  with(options,{
    E = checkIt(E, nStrata)
    I = checkIt(I, nStrata)
    R = checkIt(R, nStrata)
    return(list(H=H, E=E, I=I, R=R))
})}



#' @title Return the parameters as a list
#' @description This method dispatches on the type of `xds_obj$XH_obj[[i]]`.
#' @inheritParams ramp.xds::change_XH_pars
#' @return an **`xds`** object
#' @export
change_XH_pars.SEIR <- function(xds_obj, i=1, options=list()) {
  nHabitats <- xds_obj$nHabitats
  with(xds_obj$XH_obj[[i]], with(options,{
    xds_obj$XH_obj[[i]]$b <- b
    xds_obj$XH_obj[[i]]$c <- c
    xds_obj$XH_obj[[i]]$r <- r
    xds_obj$XH_obj[[i]]$tau <- tau
    return(xds_obj)
  }))}


#' @title Return the parameters as a list
#' @description This method dispatches on the type of `xds_obj$XH_obj[[i]]`.
#' @inheritParams ramp.xds::change_XH_inits
#' @return an **`xds`** object
#' @export
change_XH_inits.SEIR <- function(xds_obj, i=1, options=list()) {
  with(xds_obj$XH_obj[[i]], with(options,{
    xds_obj$XH_inits[[i]]$inits$H = H
    xds_obj$XH_inits[[i]]$inits$E = E
    xds_obj$XH_inits[[i]]$inits$I = I
    xds_obj$XH_inits[[i]]$inits$R = R
    return(xds_obj)
}))}





#' @title Add indices for human population to parameter list
#' @description Implements [setup_XH_ix] for the SEIR model.
#' @inheritParams ramp.xds::setup_XH_ix
#' @return none
#' @importFrom utils tail
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


#' @title Return the variables as a list
#' @description This method dispatches on the type of `xds_obj$XH_obj`
#' @inheritParams ramp.xds::get_XH_vars
#' @return a [list]
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



#' @title Make parameters for SEIR human model, with defaults
#' @param nStrata is the number of population strata
#' @param options a [list] that could overwrite defaults
#' @param tau  incubation rate
#' @param b the proportion of infective bites that cause an infection
#' @param r the the duration of an infection
#' @param c the proportion of bites on infected humans that infect a mosquito
#' @return a [list]
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

    # Ports for demographic models
    XH_obj$D_matrix = diag(0, nStrata)
    births = "zero"
    class(births) = births
    XH_obj$births = births

    return(XH_obj)
  })}





#' @title Setup XH_obj.SEIR
#' @description Implements [setup_XH_obj] for the SEIR model
#' @inheritParams ramp.xds::setup_XH_obj
#' @return a [list] vector
#' @export
setup_XH_obj.SEIR = function(Xname, xds_obj, i, options=list()){
  XH_obj <- make_XH_obj_SEIR(xds_obj$nStrata[i], options)
  class(XH_obj) <- c("SEIR", paste("SEIR_", xds_obj$xds, sep=""))
  xds_obj$XH_obj[[i]] = XH_obj
  return(xds_obj)
}


#' @title Size of effective infectious human population
#' @description Implements [F_I] for the SIS model.
#' @inheritParams ramp.xds::F_I
#' @return a [numeric] vector of length `nStrata`
#' @export
F_I.SEIR <- function(t,y, xds_obj, i) {
  I = y[xds_obj$XH_obj[[i]]$ix$I_ix]
  Y = with(xds_obj$XH_obj[[i]], c*I)
  return(Y)
}


#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SEIR model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SEIR <- function(t, y, xds_obj, i){
  with(get_XH_vars(y, xds_obj, i), return(H))
}


#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_infectivity] for the SEIR model.
#' @inheritParams ramp.xds::F_infectivity
#' @return a [numeric] vector of length `nStrata`
#' @export
F_infectivity.SEIR <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]],return(b))
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_ni] for the SIR model.
#' @inheritParams ramp.xds::F_ni
#' @return a [numeric] vector of length `nStrata`
#' @export
F_ni.SEIR <- function(vars, XH_obj) {
  with(vars, with(XH_obj, c*I/H))
}

#' @title parse the output of deSolve and return variables for the SEIR model
#' @description Implements [parse_XH_orbits] for the SEIR model
#' @inheritParams ramp.xds::parse_XH_orbits
#' @return none
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
#' @description Implements [F_prevalence] for the SEIR model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @export
F_prevalence.SEIR <- function(vars, XH_obj) {
  pr = with(vars, I/H)
  return(pr)
}


#' @title Compute the HTC for the SEIR model
#' @description Implements [HTC] for the SEIR model with demography.
#' @inheritParams ramp.xds::HTC
#' @return a [numeric] vector
#' @export
HTC.SEIR <- function(xds_obj, i) {
  with(xds_obj$XH_obj[[i]],
       HTC <- c/r,
       return(HTC)
  )
}


#' Add lines for the density of infected individuals for the SEIR model
#' @param time time points for the observations
#' @param XH a list with the outputs of parse_outputs_X_SEIR
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
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


#' Plot the density of infected individuals for the SEIR model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @export
xds_plot_X.SEIR = function(xds_obj, i=1, clrs=c("black","darkblue","darkred","darkgreen"), llty=1,  add=FALSE){
  XH = xds_obj$outputs$orbits$XH[[i]]
  time = xds_obj$outputs$time

  if(add==FALSE)
    plot(time, 0*time, type = "n", ylim = c(0, max(XH$H)),
         ylab = "No of. Infected", xlab = "Time")
  xds_lines_X_SEIR(time, XH, xds_obj$nStrata[i], clrs, llty)
}

