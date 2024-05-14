#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [dXdt] for the SEIRnew model
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SEIRnew<- function(t, y, pars, i) {
  # do not change this
  foi <- pars$FoI[[i]]
  
  # attach the variables by name
  with(list_Xvars(y, pars, i),{
    # compute H (if it isn't one of the variables)
    H <- F_H(t, y, pars, i)
    
    # expose the parameters (see make_Xpar_SEIRnew)
    with(pars$Xpar[[i]], {
      # compute the derivatives
      dS <- Births(t, H, pars,i) - foi*S + dHdt(t, S, pars,i)
      dE <- foi*S - tau*E + dHdt(t, E, pars,i)
      dI <- tau*E - r*I + dHdt(t, I, pars,i)
      dR <- r*I + dHdt(t, R, pars,i)
      
       # concatenate the derivatives
      derivs = c(dS, dE, dI, dR)
      
      # return the derivatives
      return(derivs)
    })
  })
}



#' @title Make initial values for the SEIRnew human model, with defaults
#' @param nStrata the number of strata in the model
#' @param Xopts a [list] to overwrite defaults
#' @param H0 the initial value for H
#' @param S0 the initial value for S
#' @param E0 the initial value for E
#' @param I0 the initial value for I
#' @param R0 the initial values for R
#' @return a [list]
#' @export
make_Xinits_SEIRnew = function(nStrata, Xopts = list(), H0= NULL, S0=NULL, I0=1, E0=0,R0 = 1){with(Xopts,{
  if(is.null(S0)) S0 = H0-(E0+I0+R0)
  stopifnot(is.numeric(S0))
  S = checkIt(S0, nStrata)
  E = checkIt(E0, nStrata)
  I = checkIt(I0, nStrata)
  R = checkIt(R0, nStrata)
  return(list(S=S, E=E, I=I, R =R))
})}






#' @title Setup Xinits.SEIRnew
#' @description Implements [setup_Xinits] for the SEIRnew model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SEIRnew = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_SEIRnew(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}





#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SEIRnew model.
#' @inheritParams ramp.xde::make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.SEIRnew <- function(pars, i) {with(pars,{
  
  S_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(S_ix, 1)
  
  E_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(E_ix, 1)
  
  I_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(I_ix, 1)
  
  R_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(R_ix, 1)
  
  
  
  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(S_ix=S_ix,  E_ix=E_ix, I_ix=I_ix, R_ix=R_ix)
  return(pars)
})}





#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xde::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SEIRnew <- function(y, pars, i) {
  with(pars$ix$X[[i]],
       return(list(
         S = y[S_ix],
         E = y[E_ix],
         I = y[I_ix],
         R = y[R_ix]
       )
       ))
}





#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.SEIRnew <- function(pars, i){
  with(pars$Xinits[[i]], return(c(S,E,I,R)))
}






#' @title Update inits for the SEIRnew human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.SEIRnew <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_SEIRnew(pars, list(), S0=S,E0= E, I0=I, R0=R)
    return(pars)
  })}



#' @title Make parameters for SEIRnew human model, with defaults
#' @param nStrata is the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param tau  incubation rate
#' @param b the proportion of infective bites that cause an infection
#' @param r the the duration of an infection
#' @param c the proportion of bites on infected humans that infect a mosquito
#' @return a [list]
#' @export
make_Xpar_SEIRnew = function(nStrata, Xopts=list(),
                          b=0.55, r=1/180, c=0.15,tau= 0.5){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SEIRnew")
    
    Xpar$b = checkIt(b, nStrata)
    Xpar$tau = checkIt(tau, nStrata)
    Xpar$c = checkIt(c, nStrata)
    Xpar$r = checkIt(r, nStrata)
    
    return(Xpar)
  })}





#' @title Setup Xpar.SEIRnew
#' @description Implements [setup_Xpar] for the SEIRnew model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SEIRnew = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SEIRnew(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}





#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SEIRnew <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  Y = with(pars$Xpar[[i]], c*I)
  return(Y)
}





#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SEIRnew model.
#' @inheritParams ramp.xde::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SEIRnew <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i), {
    H <- S +E+ I+R
    return(H)
  })}





#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SEIRnew model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SEIRnew <- function(y, pars, i) {
  with(pars$Xpar[[i]],{ 
    ########################
    # retrieve or compute it 
    ########################
    b = pars$Xpar[[i]]$b
    ########################
    # return it 
    ########################
    return(b)
  })
}





#' @title Parse the output of deSolve and return variables for the SEIRnew model
#' @description Implements [parse_deout_X] for the SEIRnew model
#' @inheritParams ramp.xde::parse_deout_X
#' @return none
#' @export
parse_deout_X.SEIRnew <- function(deout, pars, i) {
  time = deout[,1]
  with(pars$ix$X[[i]],{
    S = deout[,S_ix+1]
    E = deout[,E_ix+1]
    I = deout[,I_ix+1]
    R = deout[,R_ix+1]
    
    H = S+I+R
    return(list(time=time, S=S, E=E,I=I, R=R, H=H))
  })}





#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIRnew model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIRnew <- function(varslist, pars,i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}





#' @title Compute the HTC for the SEIRnew model
#' @description Implements [HTC] for the SEIRnew model with demography.
#' @inheritParams ramp.xde::HTC
#' @return a [numeric] vector
#' @export
HTC.SEIRnew <- function(pars, i) {
  with(pars$Xpar[[i]],
       HTC <- c/r,
       return(HTC)
  )
}





#' Add lines for the density of infected individuals for the SEIRnew model
#'
#' @param XH a list with the outputs of parse_deout_X_SEIRnew
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
xde_lines_X_SEIRnew = function(XH, nStrata, clrs=c("black","darkblue","darkred","darkgreen"), llty=1){
  with(XH,{
    if(nStrata==1) {
      lines(time, S, col=clrs[1], lty = llty[1])
      lines(time, E, col=clrs[2], lty = llty[1])
      lines(time, I, col=clrs[3], lty = llty[1])
      lines(time, R, col=clrs[4], lty = llty[1])
    }
    if(nStrata>1){
      if (length(clrs)==3) clrs=matrix(clrs, 3, nStrata)
      if (length(llty)==1) llty=rep(llty, nStrata)
      
      for(i in 1:nStrata){
        lines(time, S[,i], col=clrs[1,i], lty = llty[i])
        lines(time, E[,i], col=clrs[2,i], lty = llty[i])
        lines(time, I[,i], col=clrs[3,i], lty = llty[i])
        lines(time, R[,i], col=clrs[4,i], lty = llty[i])
      }
    }
  })}






#' Plot the density of infected individuals for the SEIRnew model
#'
#' @inheritParams ramp.xde::xde_plot_X
#' @export
xde_plot_X.SEIRnew = function(pars, i=1, clrs=c("black","darkblue","darkred","darkgreen"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})
  
  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "No of. Infected", xlab = "Time"))
  
  
  xde_lines_X_SEIRnew(vars$XH[[i]], pars$Hpar[[i]]$nStrata, clrs, llty)
}
