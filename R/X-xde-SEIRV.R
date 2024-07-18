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

      dS <- (1-alpha)*Births(t, H, Hpar) - foi*S + dHdt(t, S, Hpar)+ eps*R
      dE <- foi*S - tau*E + dHdt(t, E, Hpar)
      dI <- tau*E - r*I + dHdt(t, I, Hpar)
      dR <- (1-delta)*r*I -sig*R -eps*R+ dHdt(t, R, Hpar)
      dV <- alpha*Births(t, H, Hpar) + delta*r*I + sig*R + dHdt(t, V, Hpar)

      derivs = c(dS, dE, dI, dR, dV)

      return(derivs)
    })
  })
}



#' @title Make initial values for the SEIRV human model, with defaults
#' @param nStrata the number of strata in the model
#' @param Xopts a [list] to overwrite defaults
#' @param H0 the initial value for H
#' @param S0 the initial value for S
#' @param E0 the initial value for E
#' @param I0 the initial value for I
#' @param R0 the initial values for R
#' @param V0 the initial values for V
#' @return a [list]
#' @export
make_Xinits_SEIRV = function(nStrata, Xopts = list(), H0= NULL, S0=NULL, I0=1, E0=0,R0 = 1,V0 = 1){with(Xopts,{
  if(is.null(S0)) S0 = H0-(E0+I0+R0+V0)
  stopifnot(is.numeric(S0))
  S = checkIt(S0, nStrata)
  E = checkIt(E0, nStrata)
  I = checkIt(I0, nStrata)
  R = checkIt(R0, nStrata)
  V = checkIt(V0, nStrata)
  return(list(S=S, E=E, I=I, R =R,V=V))
})}






#' @title Setup Xinits.SEIRV
#' @description Implements [setup_Xinits] for the SEIRV model
#' @inheritParams ramp.xds::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SEIRV = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_SEIRV(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}





#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SEIRV model.
#' @inheritParams ramp.xds::make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.SEIRV <- function(pars, i) {with(pars,{

  S_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(S_ix, 1)

  E_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(E_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(I_ix, 1)

  R_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(R_ix, 1)

  V_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
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
#' @inheritParams ramp.xds::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.SEIRV <- function(pars, i){
  with(pars$Xinits[[i]], return(c(S,E,I,R,V)))
}






#' @title Update inits for the SEIRV human model from a vector of states
#' @inheritParams ramp.xds::update_inits_X
#' @return none
#' @export
update_inits_X.SEIRV <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_SEIRV(pars, list(), S0=S,E0= E, I0=I, R0=R,V0 =V)
    return(pars)
  })}



#' @title Make parameters for SEIRV human model, with defaults
#' @param nStrata is the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param tau  incubation rate
#' @param b the proportion of infective bites that cause an infection
#' @param r the the duration of an infection
#' @param c the proportion of bites on infected humans that infect a mosquito
#' @param alpha  travallers with immunity from malaria
#' @param eps   loss of immunity rate
#' @param delta   proportion of protected after recovery
#' @param sig  progression rate of recovered individual to protected class
#' @return a [list]
#' @export
make_Xpar_SEIRV = function(nStrata, Xopts=list(),
                          alpha =0.1,b=0.55, r=1/180, c=0.15,tau= 0.5,eps,sig=0.3,delta =0.5){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SEIRV")

    Xpar$alpha = checkIt(alpha, nStrata)
    Xpar$b = checkIt(b, nStrata)
    Xpar$c = checkIt(c, nStrata)
    Xpar$r = checkIt(r, nStrata)
    Xpar$tau = checkIt(tau, nStrata)
    Xpar$eps = checkIt(eps, nStrata)
    Xpar$sig = checkIt(sig, nStrata)
    Xpar$delta = checkIt(delta, nStrata)

    return(Xpar)
  })}





#' @title Setup Xpar.SEIRV
#' @description Implements [xde_setup_Xpar] for the SEIRV model
#' @inheritParams ramp.xds::xde_setup_Xpar
#' @return a [list] vector
#' @export
xde_setup_Xpar.SEIRV = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SEIRV(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}





#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SEIRV <- function(y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  Y = with(pars$Xpar[[i]], c*I)
  return(Y)
}





#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SEIRV model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SEIRV <- function(y, pars, i){
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
#' @description Implements [parse_outputs_X] for the SEIRV model
#' @inheritParams ramp.xds::parse_outputs_X
#' @return none
#' @export
parse_outputs_X.SEIRV <- function(outputs, pars, i) {
  time = outputs[,1]
  with(pars$ix$X[[i]],{
    S = outputs[,S_ix+1]
    E = outputs[,E_ix+1]
    I = outputs[,I_ix+1]
    R = outputs[,R_ix+1]
    V = outputs[,V_ix+1]
    H = S+E+I+R+V
    return(list(time=time, S=S, E=E,I=I, R=R, V=V,H=H))
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
#' @param XH a list with the outputs of parse_outputs_X_SEIRV
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
xde_lines_X_SEIRV = function(XH, nStrata, clrs=c("black","darkblue","darkred","darkgreen", "purple"), llty=1){
  with(XH,{
    if(nStrata==1) {
      lines(time, S, col=clrs[1], lty = llty[1])
      lines(time, E, col=clrs[2], lty = llty[1])
      lines(time, I, col=clrs[3], lty = llty[1])
      lines(time, R, col=clrs[4], lty = llty[1])
      lines(time, V, col=clrs[5], lty = llty[1])
    }
    if(nStrata>1){
      if (length(clrs)==5) clrs=matrix(clrs, 5, nStrata)
      if (length(llty)==1) llty=rep(llty, nStrata)

      for(i in 1:nStrata){
        lines(time, S[,i], col=clrs[1,i], lty = llty[i])
        lines(time, E[,i], col=clrs[2,i], lty = llty[i])
        lines(time, I[,i], col=clrs[3,i], lty = llty[i])
        lines(time, R[,i], col=clrs[4,i], lty = llty[i])
        lines(time, V[,i], col=clrs[5,i], lty = llty[i])
      }
    }
  })}






#' Plot the density of infected individuals for the SEIRV model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @export
xds_plot_X.SEIRV = function(pars, i=1, clrs=c("black","darkblue","darkred","darkgreen","purple"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "No of. Infected", xlab = "Time"))


  xde_lines_X_SEIRV(vars$XH[[i]], pars$Hpar[[i]]$nStrata, clrs, llty)
}
