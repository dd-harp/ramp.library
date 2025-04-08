#' @title Compute the derivatives for parasite infection dynamics in human population strata
#' @description Implements [dXdt] for the SIRS model
#' @inheritParams ramp.xds::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SIRS <- function(t, y, pars, i) {

  foi <- pars$FoI[[i]]
  Hpar <- pars$Hpar[[i]]

  with(list_Xvars(y, pars, i),{
    with(pars$Xpar[[i]], {

      dH <- Births(t, H, Hpar) + dHdt(t, H, Hpar)
      dI <- foi*S - r*I + dHdt(t, I, Hpar)
      dR <- r*I - gamma*R + dHdt(t, R, Hpar)

      return(c(dH, dI, dR))
    })
  })
}

#' @title DTS updating for the SIS model for human / vertebrate host infections
#' @description Implements [Update_Xt] for the SIS model
#' @inheritParams ramp.xds::Update_Xt
#' @return a [numeric] vector
#' @export
Update_Xt.SIRS<- function(t, y, pars, i) {

  ar <- pars$AR[[i]]
  Hpar <- pars$Hpar[[i]]
  with(list_Xvars(y, pars, i),{
    with(pars$Xpar[[i]], {

      St <- (1-ar)*S  + gamma*R + dHdt(t, S, Hpar) + Births(t, H, Hpar)
      It <- (1-r)*I + ar*S + dHdt(t, I, Hpar)
      Rt <- (1-gamma)*R + r*I + dHdt(t, R, Hpar)

      return(c(S=unname(St), I=unname(It), R = unname(Rt)))
    })
  })
}

#' @title Compute the steady states for the  dts SEIS model as a function of the daily EIR
#' @description Compute the steady state of the  dts SIS model as a function of the daily eir.
#' @inheritParams ramp.xds::dts_steady_state_X
#' @return the steady states as a named vector
#' @export
dts_steady_state_X.SIRS= function(ar, H, Xpar){with(Xpar,{
  Iteq = (ar*H*gamma)/(r*gamma + ar*(r+gamma))
  Rteq = (ar*H*r)/(r*gamma + ar*(r+gamma))
  Steq = H -Iteq-Rteq

  return(c(S=Steq, I=Iteq, R=Rteq))
})}


#' @title Return the parameters as a list
#' @description This method dispatches on the type of `pars$Xpar[[i]]`.
#' @inheritParams ramp.xds::set_Xpars
#' @return an **`xds`** object
#' @export
set_Xpars.SIRS <- function(pars, i=1, Xopts=list()) {
  nHabitats <- pars$nHabitats
  with(pars$Xpar[[i]], with(Xopts,{
    pars$Xpar[[i]]$b <- b
    pars$Xpar[[i]]$c <- c
    pars$Xpar[[i]]$r <- r
    pars$Xpar[[i]]$gamma <- gamma
    return(pars)
  }))}


#' @title Return the parameters as a list
#' @description This method dispatches on the type of `pars$Xpar[[i]]`.
#' @inheritParams ramp.xds::set_Xinits
#' @return an **`xds`** object
#' @export
set_Xinits.SIRS <- function(pars, i=1, Xopts=list()) {
  with(pars$Xpar[[i]], with(Xopts,{
    pars$Xinits[[i]]$S = S
    pars$Xinits[[i]]$I = I
    pars$Xinits[[i]]$R = R
    return(pars)
  }))}


#' @title Compute the steady states for the SIRS model as a function of the daily EIR
#' @description Compute the steady state of the SIRS model as a function of the daily eir.
#' @inheritParams  ramp.xds::xde_steady_state_X
#' @return the steady states as a named vector
#' @export
xde_steady_state_X.SIRS = function(foi, H, Xpar){with(Xpar,{
  Ieq = foi*H*gamma/((r*gamma) + foi*(gamma+r))
  Req = foi*r*H/((r*gamma) +foi*(gamma+r))
  Seq = H-Ieq-Req
  return(c(S=as.vector(Seq), I=as.vector(Ieq), R = as.vector(Req)))
})}

#' @title Make initial values for the SIRS human model, with defaults
#' @param nStrata the number of strata in the model
#' @param Xopts a [list] to overwrite defaults
#' @param H the initial value for H
#' @param I the initial value for I
#' @param R the initial values for R
#' @return a [list]
#' @export
make_Xinits_SIRS = function(nStrata, H,  Xopts = list(), I=1, R=0){with(Xopts,{
  H = H
  I = checkIt(I, nStrata)
  R = checkIt(R, nStrata)
  return(list(H=H, I=I, R=R))
})}



#' @title Setup Xinits.SIRS
#' @description Implements [setup_Xinits] for the SIRS model
#' @inheritParams ramp.xds::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SIRS = function(pars, H, i, Xopts=list()){
  pars$Xinits[[i]] = make_Xinits_SIRS(pars$nStrata[i], H, Xopts)
  return(pars)
}


#' @title Add indices for human population to parameter list
#' @description Implements [setup_Xix] for the SIRS model.
#' @inheritParams ramp.xds::setup_Xix
#' @return none
#' @importFrom utils tail
#' @export
setup_Xix.SIRS <- function(pars, i) {with(pars,{

  H_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(H_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(I_ix, 1)

  R_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(R_ix, 1)

  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(H_ix=H_ix, I_ix=I_ix, R_ix=R_ix)
  return(pars)
})}


#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xds::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SIRS <- function(y, pars, i) {
  with(pars$ix$X[[i]],{
    H = y[H_ix]
    I = y[I_ix]
    R = y[R_ix]
    S = H-I-R
    return(list(S=S, I=I, R=R, H=H))})
}


#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xds::get_Xinits
#' @return a [numeric] vector
#' @export
get_Xinits.SIRS <- function(pars, i){
  pars$Xinits[[i]]
}


#' @title Update inits for the SIRS human model from a vector of states
#' @inheritParams ramp.xds::update_Xinits
#' @return none
#' @export
update_Xinits.SIRS <- function(pars, y, i) {
  with(list_Xvars(y, pars, i),{
    pars$Xinits[[i]] = make_Xinits_SIRS(pars$nStrata[i], H, list(), I=I, R=R)
    return(pars)
  })}


#' @title Make parameters for SIRS human model, with defaults
#' @param nStrata is the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param b the proportion of infective bites that cause an infection
#' @param r the the duration of an infection
#' @param c the proportion of bites on infected humans that infect a mosquito
#' @param gamma the rate of loss of immunity
#' @return a [list]
#' @export
make_Xpar_SIRS = function(nStrata, Xopts=list(),
                          b=0.55, r=1/180, c=0.15,gamma=0.5){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SIRS")

    Xpar$b = checkIt(b, nStrata)
    Xpar$c = checkIt(c, nStrata)
    Xpar$r = checkIt(r, nStrata)
    Xpar$gamma = checkIt(gamma, nStrata)

    return(Xpar)
})}



#' @title Setup Xpar.SIRS
#' @description Implements [setup_Xpar] for the SIRS model
#' @inheritParams ramp.xds::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SIRS = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SIRS(pars$nStrata[i], Xopts)
  return(pars)
}





#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SIRS <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  Y = with(pars$Xpar[[i]], c*I)
  return(Y)
}





#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SIRS model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SIRS <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i), return(H))
}





#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SIRS model.
#' @inheritParams ramp.xds::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SIRS <- function(y, pars, i) {
  with(pars$Xpar[[i]],return(b))
}


#' @title Parse the output of deSolve and return variables for the SIRS model
#' @description Implements [parse_Xorbits] for the SIRS model
#' @inheritParams ramp.xds::parse_Xorbits
#' @return none
#' @export
parse_Xorbits.SIRS <- function(outputs, pars, i) {with(pars$ix$X[[i]],{
    H = outputs[,H_ix]
    I = outputs[,I_ix]
    R = outputs[,R_ix]
    S = H-I-R
    ni <- pars$Xpar[[i]]$c*I/H
    true_pr <- I/H
    return(list(S=S, I=I, R=R, H=H,ni=ni, true_pr= true_pr))
  })}





#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SIRS model.
#' @inheritParams ramp.xds::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SIRS <- function(vars, Xpar) {
  pr = with(vars, I/H)
  return(pr)
}





#' @title Compute the HTC for the SIRS model
#' @description Implements [HTC] for the SIRS model with demography.
#' @inheritParams ramp.xds::HTC
#' @return a [numeric] vector
#' @export
HTC.SIRS <- function(pars, i) {
  with(pars$Xpar[[i]],
       HTC <- c/r,
       return(HTC)
  )
}



#' Add lines for the density of infected individuals for the SIRS model
#'
#' @param times time points for the observations
#' @param XH a list with the outputs of parse_outputs_X_SIRS
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
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





#' Plot the density of infected individuals for the SIRS model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @export
xds_plot_X.SIRS = function(pars, i=1, clrs=c("darkblue","darkred","darkgreen"), llty=1, add=FALSE){
  XH = pars$outputs$orbits$XH[[i]]
  times = pars$outputs$time

  if(add==FALSE)
    plot(times, 0*times, type = "n", ylim = c(0, max(XH$H)),
         ylab = "No of. Infected", xlab = "Time")
  xds_lines_X_SIRS(times, XH, pars$nStrata[i], clrs, llty)
}

