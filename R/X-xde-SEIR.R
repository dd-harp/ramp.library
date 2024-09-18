#' @title Compute the derivatives for parasite infection dynamics in human population strata
#' @description Implements [dXdt] for the SEIR model
#' @inheritParams ramp.xds::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SEIR<- function(t, y, pars, i) {
  # do not change this
  foi <- pars$FoI[[i]]
  Hpar <- pars$Hpar[[i]]

  # attach the variables by name
  with(list_Xvars(y, pars, i),{
    # compute H (if it isn't one of the variables)

    # expose the parameters (see create_Xpar_SEIR)
    with(pars$Xpar[[i]], {
      # compute the derivatives
      dS <- Births(t, H, Hpar) - foi*S + dHdt(t, S, Hpar)
      dE <- foi*S - tau*E + dHdt(t, E, Hpar)
      dI <- tau*E - r*I + dHdt(t, I, Hpar)
      dR <- r*I + dHdt(t, R, Hpar)

       # concatenate the derivatives
      derivs = c(dS, dE, dI, dR)

      # return the derivatives
      return(derivs)
    })
  })
}


#' @title DTS updating for the SIS model for human / vertebrate host infections
#' @description Implements [Update_Xt] for the SIS model
#' @inheritParams ramp.xds::Update_Xt
#' @return a [numeric] vector
#' @export
Update_Xt.SEIR<- function(t, y, pars, i) {

  ar <- pars$AR[[i]]
  Hpar <- pars$Hpar[[i]]
  with(list_Xvars(y, pars, i),{
    with(pars$Xpar[[i]], {

      St <- (1-ar)*S  + dHdt(t, S, Hpar) + Births(t, H, Hpar)
      Et <- a*S +(1-tau)*E
      It <- (1-r)*I + tau*E + dHdt(t, I, Hpar)
      Rt <- R + r*I + dHdt(t, R, Hpar)

      return(c(S=unname(St), E = unname(It), I=unname(It), R = unname(Rt)))
    })
  })
}

#' @title Compute the steady states for the  dts SEIS model as a function of the daily EIR
#' @description Compute the steady state of the  dts SIS model as a function of the daily eir.
#' @inheritParams ramp.xds::dts_steady_state_X
#' @return the steady states as a named vector
#' @export
dts_steady_state_X.SEIR= function(ar, H, Xpar){with(Xpar,{
  Steq = 0
  Eteq = 0
  Iteq = 0
  Rteq = H -Steq - Eteq-Iteq

  return(c(S=Steq,  E= Eteq, I=Iteq, R=Rteq))
})}


#' @title Compute the steady states for the SIRS model as a function of the daily EIR
#' @description Compute the steady state of the SIRS model as a function of the daily eir.
#' @inheritParams  ramp.xds::xde_steady_state_X
#' @return the steady states as a named vector
#' @export
xde_steady_state_X.SEIR = function(foi, H, Xpar){with(Xpar,{
  Eeq = 0
  Ieq = 0
  Req =H
  Seq = H-Ieq-Req
  return(c(S=as.vector(Seq),E = as.vector(Eeq), I=as.vector(Ieq), R = as.vector(Req)))
})}


#' @title Make initial values for the SEIR human model, with defaults
#' @param nStrata the number of strata in the model
#' @param H the initial value for H
#' @param Xopts a [list] to overwrite defaults
#' @param E the initial value for E
#' @param I the initial value for I
#' @param R the initial values for R
#' @return a [list]
#' @export
create_Xinits_SEIR = function(nStrata, H, Xopts = list(), I=1, E=0, R = 1){with(Xopts,{
  S = checkIt(H-I-E-R, nStrata)
  E = checkIt(E, nStrata)
  I = checkIt(I, nStrata)
  R = checkIt(R, nStrata)
  return(list(S=S, E=E, I=I, R =R))
})}






#' @title Setup Xinits.SEIR
#' @description Implements [make_Xinits] for the SEIR model
#' @inheritParams ramp.xds::make_Xinits
#' @return a [list] vector
#' @export
make_Xinits.SEIR = function(pars, H, i, Xopts=list()){
  pars$Xinits[[i]] = create_Xinits_SEIR(pars$nStrata[i], H, Xopts)
  return(pars)
}





#' @title Add indices for human population to parameter list
#' @description Implements [make_X_indices] for the SEIR model.
#' @inheritParams ramp.xds::make_X_indices
#' @return none
#' @importFrom utils tail
#' @export
make_X_indices.SEIR <- function(pars, i) {with(pars,{

  S_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(S_ix, 1)

  E_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(E_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(I_ix, 1)

  R_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(R_ix, 1)



  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(S_ix=S_ix,  E_ix=E_ix, I_ix=I_ix, R_ix=R_ix)
  return(pars)
})}





#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xds::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SEIR <- function(y, pars, i) {
  with(pars$ix$X[[i]],{
       S = y[S_ix]
       E = y[E_ix]
       I = y[I_ix]
       R = y[R_ix]
       H =S+E+I+R
       return(list(S=S,E=E,I=I,R=R,H=H))})
}





#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xds::get_Xinits
#' @return a [numeric] vector
#' @export
get_Xinits.SEIR <- function(pars, i){
  pars$Xinits[[i]]
}






#' @title Update inits for the SEIR human model from a vector of states
#' @inheritParams ramp.xds::update_Xinits
#' @return none
#' @export
update_Xinits.SEIR <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = create_Xinits_SEIR(pars$nStrata[i], pars$H0, list(), E=E, I=I, R=R)
    return(pars)
  })}



#' @title Make parameters for SEIR human model, with defaults
#' @param nStrata is the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param tau  incubation rate
#' @param b the proportion of infective bites that cause an infection
#' @param r the the duration of an infection
#' @param c the proportion of bites on infected humans that infect a mosquito
#' @return a [list]
#' @export
create_Xpar_SEIR = function(nStrata, Xopts=list(),
                          b=0.55, r=1/180, c=0.15,tau= 0.5){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SEIR")

    Xpar$b = checkIt(b, nStrata)
    Xpar$tau = checkIt(tau, nStrata)
    Xpar$c = checkIt(c, nStrata)
    Xpar$r = checkIt(r, nStrata)

    return(Xpar)
  })}





#' @title Setup Xpar.SEIR
#' @description Implements [make_Xpar] for the SEIR model
#' @inheritParams ramp.xds::make_Xpar
#' @return a [list] vector
#' @export
make_Xpar.SEIR = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = create_Xpar_SEIR(pars$nStrata[i], Xopts)
  return(pars)
}





#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SEIR <- function(t,y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  Y = with(pars$Xpar[[i]], c*I)
  return(Y)
}


#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SEIR model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SEIR <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i), return(H))
}





#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SEIR model.
#' @inheritParams ramp.xds::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SEIR <- function(y, pars, i) {
  with(pars$Xpar[[i]],return(b))
}


#' @title Parse the output of deSolve and return variables for the SEIR model
#' @description Implements [parse_Xorbits] for the SEIR model
#' @inheritParams ramp.xds::parse_Xorbits
#' @return none
#' @export
parse_Xorbits.SEIR <- function(outputs, pars, i) {with(pars$ix$X[[i]],{
    S = outputs[,S_ix]
    E = outputs[,E_ix]
    I = outputs[,I_ix]
    R = outputs[,R_ix]
    H = S+I+R
    ni <- pars$Xpar[[i]]$c*I/H
    true_pr <- (I+E)/H
    return(list(S=S, E=E,I=I, R=R, H=H,ni=ni, true_pr= true_pr))
  })}





#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIR model.
#' @inheritParams ramp.xds::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIR <- function(vars, Xpar) {
  pr = with(vars, I/H)
  return(pr)
}





#' @title Compute the HTC for the SEIR model
#' @description Implements [HTC] for the SEIR model with demography.
#' @inheritParams ramp.xds::HTC
#' @return a [numeric] vector
#' @export
HTC.SEIR <- function(pars, i) {
  with(pars$Xpar[[i]],
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
xds_plot_X.SEIR = function(pars, i=1, clrs=c("black","darkblue","darkred","darkgreen"), llty=1,  add=FALSE){
  XH = pars$outputs$orbits$XH[[i]]
  time = pars$outputs$time

  if(add==FALSE)
    plot(time, 0*time, type = "n", ylim = c(0, max(XH$H)),
         ylab = "No of. Infected", xlab = "Time")
  xds_lines_X_SEIR(time, XH, pars$nStrata[i], clrs, llty)
}

