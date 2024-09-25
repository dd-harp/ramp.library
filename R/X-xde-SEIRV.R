#' @title Compute the derivatives for parasite infection dynamics in human population strata
#' @description Implements [dXdt] for the SEIRV model
#' @inheritParams ramp.xds::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SEIRV<- function(t, y, pars, i) {
  foi <- pars$FoI[[i]]
  Hpar <- pars$Hpar[[i]]

  with(list_Xvars(y, pars, i),{
    with(pars$Xpar[[i]], {

      dS <- (1-alpha)*Births(t, H, Hpar) - foi*S + dHdt(t, S, Hpar)+ gamma*R
      dE <- foi*S - tau*E + dHdt(t, E, Hpar)
      dI <- tau*E - r*I + dHdt(t, I, Hpar)
      dR <- (1-varepsilon)*r*I -gamma*R+ dHdt(t, R, Hpar)
      dV <- alpha*Births(t, H, Hpar) + varepsilon*r*I + dHdt(t, V, Hpar)

      derivs = c(dS, dE, dI, dR, dV)

      return(derivs)
    })
  })
}

#' @title DTS updating for the SIS model for human / vertebrate host infections
#' @description Implements [Update_Xt] for the SIS model
#' @inheritParams ramp.xds::Update_Xt
#' @return a [numeric] vector
#' @export
Update_Xt.SEIRV<- function(t, y, pars, i) {

  ar <- pars$AR[[i]]
  Hpar <- pars$Hpar[[i]]
  with(list_Xvars(y, pars, i),{
    with(pars$Xpar[[i]], {

      St <- (1-ar)*S  + gamma*R + dHdt(t, S, Hpar) + Births(t, H, Hpar)
      Et <- a*S +(1-tau)*E
      It <- (1-r)*I + tau*E + dHdt(t, I, Hpar)
      Rt <- (1-gamma)*R + (1-varepsilon)*r*I + dHdt(t, R, Hpar)
      Vt  <- V + varepsilon*r* I +  dHdt(t, V, Hpar)

      return(c(S=unname(St), E = unname(It), I=unname(It), R = unname(Rt), V = unname(Vt)))
    })
  })
}

#' @title Return the parameters as a list
#' @description This method dispatches on the type of `pars$Xpar[[i]]`.
#' @inheritParams ramp.xds::set_Xpars
#' @return an **`xds`** object
#' @export
set_Xpars.SEIRV <- function(pars, i=1, Xopts=list()) {
  nHabitats <- pars$nHabitats
  with(pars$Xpar[[i]], with(Xopts,{
    pars$Xpar[[i]]$b <- b
    pars$Xpar[[i]]$c <- c
    pars$Xpar[[i]]$r <- r
    pars$Xpar[[i]]$tau <- tau
    pars$Xpar[[i]]$gamma <- gamma
    return(pars)
  }))}


#' @title Return the parameters as a list
#' @description This method dispatches on the type of `pars$Xpar[[i]]`.
#' @inheritParams ramp.xds::set_Xinits
#' @return an **`xds`** object
#' @export
set_Xinits.SEIRV <- function(pars, i=1, Xopts=list()) {
  with(pars$Xpar[[i]], with(Xopts,{
    pars$Xinits[[i]]$S = S
    pars$Xinits[[i]]$E = E
    pars$Xinits[[i]]$I = I
    pars$Xinits[[i]]$R = R
    pars$Xinits[[i]]$V = V
    return(pars)
  }))}

#' @title Compute the steady states for the  dts SEIS model as a function of the daily EIR
#' @description Compute the steady state of the  dts SIS model as a function of the daily eir.
#' @inheritParams ramp.xds::dts_steady_state_X
#' @return the steady states as a named vector
#' @export
dts_steady_state_X.SEIRV= function(ar, H, Xpar){with(Xpar,{
  Steq = 0
  Eteq = 0
  Iteq = 0
  Rteq = 0
  Vteq = H -Steq - Eteq-Iteq-Rteq

  return(c(S=Steq,  E= Eteq, I=Iteq, R=Rteq, V= Vteq))
})}



#' @title Make initial values for the SEIRV human model, with defaults
#' @param nStrata the number of strata in the model
#' @param Xopts a [list] to overwrite defaults
#' @param H the initial value for H
#' @param E the initial value for E
#' @param I the initial value for I
#' @param R the initial values for R
#' @param V the initial values for V
#' @return a [list]
#' @export
create_Xinits_SEIRV = function(nStrata,H, Xopts = list(), I=1, E=0,R = 1,V = 1){with(Xopts,{
  S = checkIt(H-E-I-R-V, nStrata)
  E = checkIt(E, nStrata)
  I = checkIt(I, nStrata)
  R = checkIt(R, nStrata)
  V = checkIt(V, nStrata)
  return(list(S=S, E=E, I=I, R =R,V=V))
})}






#' @title Setup Xinits.SEIRV
#' @description Implements [make_Xinits] for the SEIRV model
#' @inheritParams ramp.xds::make_Xinits
#' @return a [list] vector
#' @export
make_Xinits.SEIRV = function(pars, H, i, Xopts=list()){
  pars$Xinits[[i]] = create_Xinits_SEIRV(pars$nStrata[i], H, Xopts)
  return(pars)
}





#' @title Add indices for human population to parameter list
#' @description Implements [make_X_indices] for the SEIRV model.
#' @inheritParams ramp.xds::make_X_indices
#' @return none
#' @importFrom utils tail
#' @export
make_X_indices.SEIRV <- function(pars, i) {with(pars,{

  S_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(S_ix, 1)

  E_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(E_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(I_ix, 1)

  R_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(R_ix, 1)

  V_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(V_ix, 1)



  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(S_ix=S_ix,  E_ix=E_ix, I_ix=I_ix, R_ix=R_ix, V_ix=V_ix)
  return(pars)
})}





#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xds::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SEIRV <- function(y, pars, i) {
    with(pars$ix$X[[i]],{
      S = y[S_ix]
      E = y[E_ix]
      I = y[I_ix]
      R = y[R_ix]
      V = y[V_ix]
      H =S+E+I+R+V
      return(list(S=S,E=E,I=I,R=R,V=V,H=H))})
}




#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xds::get_Xinits
#' @return a [numeric] vector
#' @export
get_Xinits.SEIRV <- function(pars, i){
  pars$Xinits[[i]]
}



#' @title Update inits for the SEIRV human model from a vector of states
#' @inheritParams ramp.xds::update_Xinits
#' @return none
#' @export
update_Xinits.SEIRV <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars$Xinits[[i]] = create_Xinits_SEIRV(pars$nStrata[i], pars$H0,  list(), E= E, I=I, R=R,V =V)
    return(pars)
  })}



#' @title Make parameters for SEIRV human model, with defaults
#' @param nStrata is the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param tau  incubation rate
#' @param b the proportion of infective bites that cause an infection
#' @param r the rate infections clear
#' @param c the proportion of bites on infected humans that infect a mosquito
#' @param alpha  proportion of vaccinated humans
#' @param gamma   loss of immunity rate
#' @param varepsilon proportion of recovered humans that are protected from the pathogen
#' @return a [list]
#' @export
create_Xpar_SEIRV = function(nStrata, Xopts=list(),
                          alpha =0.1,b=0.55, r=1/180, c=0.15,tau= 0.5,varepsilon=0.5,gamma=0.5){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SEIRV")

    Xpar$alpha = checkIt(alpha, nStrata)
    Xpar$b = checkIt(b, nStrata)
    Xpar$c = checkIt(c, nStrata)
    Xpar$r = checkIt(r, nStrata)
    Xpar$tau = checkIt(tau, nStrata)
    Xpar$gamma = checkIt(gamma, nStrata)
    Xpar$varepsilon = checkIt(varepsilon, nStrata)

    return(Xpar)
  })}

#' @title Compute the steady states for the SEIRV model as a function of the daily EIR
#' @description Compute the steady state of the SIS model as a function of the daily eir.
#' @inheritParams ramp.xds::xde_steady_state_X
#' @return the steady states as a named vector
#' @export
xde_steady_state_X.SEIRV = function(foi, H, Xpar){with(Xpar,{
  Ieq = 0
  Seq = 0
  Eeq = 0
  Req = 0
  Veq = H
  return(c(S=as.vector(Seq), E =as.vector(Eeq),I=as.vector(Ieq), R =as.vector(Req), V=as.vector(Veq)))
})}





#' @title Setup Xpar.SEIRV
#' @description Implements [make_Xpar] for the SEIRV model
#' @inheritParams ramp.xds::make_Xpar
#' @return a [list] vector
#' @export
make_Xpar.SEIRV = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = create_Xpar_SEIRV(pars$nStrata[i], Xopts)
  return(pars)
}





#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SEIRV <- function(t,y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  Y = with(pars$Xpar[[i]], c*I)
  return(Y)
}





#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SEIRV model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SEIRV <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i), {
    H <- S +E+ I+R+V
    return(H)
  })}





#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SEIRV model.
#' @inheritParams ramp.xds::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SEIRV <- function(y, pars, i) {
  with(pars$Xpar[[i]],return(b))
}


#' @title Parse the output of deSolve and return variables for the SEIRV model
#' @description Implements [parse_Xorbits] for the SEIRV model
#' @inheritParams ramp.xds::parse_Xorbits
#' @return none
#' @export
parse_Xorbits.SEIRV <- function(outputs, pars, i) {with(pars$ix$X[[i]],{
    S = outputs[,S_ix]
    E = outputs[,E_ix]
    I = outputs[,I_ix]
    R = outputs[,R_ix]
    V = outputs[,V_ix]
    H = S+E+I+R+V
    ni <- pars$Xpar[[i]]$c*I/H
    true_pr <- (E+I)/H
    return(list(S=S, E=E,I=I, R=R, V=V,H=H,ni=ni, true_pr = true_pr))
  })}





#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIRV model.
#' @inheritParams ramp.xds::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIRV <- function(vars, Xpar) {
  pr = with(vars, I/H)
  return(pr)
}





#' @title Compute the HTC for the SEIRV model
#' @description Implements [HTC] for the SEIRV model with demography.
#' @inheritParams ramp.xds::HTC
#' @return a [numeric] vector
#' @export
HTC.SEIRV <- function(pars, i) {
  with(pars$Xpar[[i]],
       HTC <- c/r,
       return(HTC)
  )
}





#' Add lines for the density of infected individuals for the SEIRV model
#'
#' @param time time points for the observations
#' @param XH a list with the outputs of parse_outputs_X_SEIRV
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
xds_lines_X_SEIRV = function(time, XH, nStrata, clrs=c("black","darkblue","darkred","darkgreen", "purple"), llty=1){
  if (length(llty)< nStrata) llty = rep(llty, nStrata)
  with(XH,{
    if(nStrata==1) {
      lines(time, S, col=clrs[1], lty = llty[1])
      lines(time, E, col=clrs[2], lty = llty[1])
      lines(time, I, col=clrs[3], lty = llty[1])
      lines(time, R, col=clrs[4], lty = llty[1])
      lines(time, V, col=clrs[5], lty = llty[1])
    } else {
      for(i in 1:nStrata)
        lines(time, S[,i], col=clrs[1,i], lty = llty[i])
        lines(time, E[,i], col=clrs[2,i], lty = llty[i])
        lines(time, I[,i], col=clrs[3,i], lty = llty[i])
        lines(time, R[,i], col=clrs[4,i], lty = llty[i])
        lines(time, V[,i], col=clrs[5,i], lty = llty[i])
      }
  })}






#' Plot the density of infected individuals for the SEIRV model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @export
xds_plot_X.SEIRV = function(pars, i=1, clrs=c("black","darkblue","darkred","darkgreen","purple"), llty=1, add=FALSE){
  XH = pars$outputs$orbits$XH[[i]]
  time = pars$outputs$time

  if(add==FALSE)
    plot(time, 0*time, type = "n", ylim = c(0, max(XH$H)),
         ylab = "No of. Infected", xlab = "Time")
  xds_lines_X_SEIRV(time, XH, pars$nStrata[i], clrs, llty)
}
