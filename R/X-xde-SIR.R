#' @title Compute the derivatives for parasite infection dynamics in human population strata
#' @description Implements [dXdt] for the SIR model
#' @inheritParams ramp.xds::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SIR<- function(t, y, pars, i) {

  foi <- pars$FoI[[i]]

  with(list_Xvars(y, pars, i),{
    Hpar <- pars$Hpar[[i]]
    with(pars$Xpar[[i]], {

      dS <- Births(t, H, Hpar) + dHdt(t, S, Hpar) - foi*S
      dI <- foi*S - r*I + dHdt(t, I, Hpar)
      dR <- r*I + dHdt(t, R, Hpar)

      return(c(dS, dI, dR))
    })
  })
}

#' @title DTS updating for the SIS model for human / vertebrate host infections
#' @description Implements [Update_Xt] for the SIS model
#' @inheritParams ramp.xds::Update_Xt
#' @return a [numeric] vector
#' @export
Update_Xt.SIR<- function(t, y, pars, i) {

  ar <- pars$AR[[i]]
  Hpar <- pars$Hpar[[i]]
  with(list_Xvars(y, pars, i),{
    with(pars$Xpar[[i]], {

      St <- (1-ar)*S  + dHdt(t, S, Hpar) + Births(t, H, Hpar)
      It <- (1-r)*I + ar*S + dHdt(t, I, Hpar)
      Rt <- R + r*I + dHdt(t, R, Hpar)

      return(c(S=unname(St), I=unname(It), R = unname(Rt)))
    })
  })
}

#' @title Compute the steady states for the  dts SEIS model as a function of the daily EIR
#' @description Compute the steady state of the  dts SIS model as a function of the daily eir.
#' @inheritParams ramp.xds::dts_steady_state_X
#' @return the steady states as a named vector
#' @export
dts_steady_state_X.SIR = function(ar, H, Xpar){with(Xpar,{
  Iteq = 0
  Steq = 0
  Rteq = H-Iteq-Steq

  return(c(S=Steq, I=Iteq, R=Rteq))
})}

#' @title Compute the steady states for the SIR model as a function of the daily EIR
#' @description Compute the steady state of the SIR model as a function of the daily eir.
#' @inheritParams ramp.xds::xde_steady_state_X
#' @return the steady states as a named vector
#' @export
xde_steady_state_X.SIR = function(foi, H, Xpar){with(Xpar,{
  Ieq = 0
  Req = H
  Seq = 0
  return(c(S=as.vector(Seq), I=as.vector(Ieq), R=as.vector(Req)))
})}

#' @title Make initial values for the SIR human model, with defaults
#' @param nStrata the number of strata in the model
#' @param Xopts a [list] to overwrite defaults
#' @param H the initial value for H
#' @param I the initial value for I
#' @param R the initial values for R
#' @return a [list]
#' @export
create_Xinits_SIR = function(nStrata,H, Xopts = list(),I=1, R=0){with(Xopts,{
  S = checkIt(H-I-R, nStrata)
  I = checkIt(I, nStrata)
  R = checkIt(R, nStrata)
  return(list(S=S, I=I, R =R))
})}


#' @title Return the parameters as a list
#' @description This method dispatches on the type of `pars$Xpar[[i]]`.
#' @inheritParams ramp.xds::set_Xpars
#' @return an **`xds`** object
#' @export
set_Xpars.SIR <- function(pars, i=1, Xopts=list()) {
  nHabitats <- pars$nHabitats
  with(pars$Xpar[[i]], with(Xopts,{
    pars$Xpar[[i]]$b <- b
    pars$Xpar[[i]]$c <- c
    pars$Xpar[[i]]$r <- r
    return(pars)
  }))}


#' @title Return the parameters as a list
#' @description This method dispatches on the type of `pars$Xpar[[i]]`.
#' @inheritParams ramp.xds::set_Xinits
#' @return an **`xds`** object
#' @export
set_Xinits.SIR <- function(pars, i=1, Xopts=list()) {
  with(pars$Xpar[[i]], with(Xopts,{
    pars$Xinits[[i]]$S = S
    pars$Xinits[[i]]$I = I
    pars$Xinits[[i]]$R = R
    return(pars)
  }))}





#' @title Setup Xinits.SIR
#' @description Implements [make_Xinits] for the SIR model
#' @inheritParams ramp.xds::make_Xinits
#' @return a [list] vector
#' @export
make_Xinits.SIR = function(pars, H, i, Xopts=list()){
  pars$Xinits[[i]] = create_Xinits_SIR(pars$nStrata[i], H, Xopts)
  return(pars)
}



#' @title Add indices for human population to parameter list
#' @description Implements [make_X_indices] for the SIR model.
#' @inheritParams ramp.xds::make_X_indices
#' @return none
#' @importFrom utils tail
#' @export
make_X_indices.SIR <- function(pars, i) {with(pars,{

  S_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(S_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(I_ix, 1)

  R_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(R_ix, 1)



  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(S_ix=S_ix, I_ix=I_ix, R_ix=R_ix)
  return(pars)
})}





#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xds::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SIR <- function(y, pars, i) {
  with(pars$ix$X[[i]],{
    S = y[S_ix]
    I = y[I_ix]
    R = y[R_ix]
    H =S+I+R
    return(list(S=S, I=I, R=R, H=H))})
}



#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xds::get_Xinits
#' @return a [numeric] vector
#' @export
get_Xinits.SIR <- function(pars, i){
  pars$Xinits[[i]]
}






#' @title Update inits for the SIR human model from a vector of states
#' @inheritParams ramp.xds::update_Xinits
#' @return none
#' @export
update_Xinits.SIR <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars$Xinits[[i]] = create_Xinits_SIR(pars$nStrata[i], pars$H0, list(), I=I, R=R)
    return(pars)
  })}



#' @title Make parameters for SIR human model, with defaults
#' @param nStrata is the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param b the proportion of infective bites that cause an infection
#' @param r the the duration of an infection
#' @param c the proportion of bites on infected humans that infect a mosquito
#' @return a [list]
#' @export
create_Xpar_SIR = function(nStrata, Xopts=list(),
                          b=0.55, r=1/180, c=0.15){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SIR")

    Xpar$b = checkIt(b, nStrata)
    Xpar$c = checkIt(c, nStrata)
    Xpar$r = checkIt(r, nStrata)

    return(Xpar)
})}


#' @title Setup Xpar.SIR
#' @description Implements [make_Xpar] for the SIR model
#' @inheritParams ramp.xds::make_Xpar
#' @return a [list] vector
#' @export
make_Xpar.SIR = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = create_Xpar_SIR(pars$nStrata[i], Xopts)
  return(pars)
}





#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SIR <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  Y = with(pars$Xpar[[i]], c*I)
  return(Y)
}





#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SIR model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SIR <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i),
    return(H))
  }





#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SIR model.
#' @inheritParams ramp.xds::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SIR <- function(y, pars, i) {
  with(pars$Xpar[[i]], b)
}

#' @title Parse the output of deSolve and return variables for the SIR model
#' @description Implements [parse_Xorbits] for the SIR model
#' @inheritParams ramp.xds::parse_Xorbits
#' @return none
#' @export
parse_Xorbits.SIR <- function(outputs, pars, i) {with(pars$ix$X[[i]],{
  S <- outputs[,S_ix]
  I <- outputs[,I_ix]
  R <- outputs[,R_ix]
  H <- S+I+R
  ni <- pars$Xpar[[i]]$c*I/H
  true_pr <- I/H
  vars <- list(S=S, I=I,R=R, H=H, ni=ni, true_pr=true_pr)
  return(vars)
})}





#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SIR model.
#' @inheritParams ramp.xds::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SIR <- function(vars, Xpar) {
  pr = with(vars, I/H)
  return(pr)
}





#' @title Compute the HTC for the SIR model
#' @description Implements [HTC] for the SIR model with demography.
#' @inheritParams ramp.xds::HTC
#' @return a [numeric] vector
#' @export
HTC.SIR <- function(pars, i) {
  with(pars$Xpar[[i]],
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
#' @export
xds_plot_X.SIR = function(pars, i=1, clrs=c("darkblue","darkred","darkgreen"), llty=1,  add=FALSE){
  XH = pars$outputs$orbits$XH[[i]]
  time = pars$outputs$time

  if(add==FALSE)
         plot(time, 0*time, type = "n", ylim = c(0, max(XH$H)),
              ylab = "No of. Infected", xlab = "Time")

  xds_lines_X_SIR(time, XH, pars$nStrata[i], clrs, llty)
}

