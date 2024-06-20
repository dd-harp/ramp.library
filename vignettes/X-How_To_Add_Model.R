## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [dXdt] for the newXname model
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.newXname <- function(t, y, pars, i) {

  # do not change this
  foi <- pars$FoI[[i]]

  # attach the variables by name 
  with(list_Xvars(y, pars, i),{
    # compute H (if it isn't one of the variables) 
    H <- F_H(t, y, pars, i)

    # expose the parameters (see make_Xpar_newXname) 
    with(pars$Xpar[[i]], {
      # compute the derivatives 
      dX1 <- ... 
      dX2 <- ...
      ... 
      
      # concatenate the derivatives 
      derivs = c(dX1, dX2, ...) 
      
      # return the derivatives 
      return(derivs)
    })
  })
}

## -----------------------------------------------------------------------------
#' @title Make initial values for the newXname human model, with defaults
#' @param nStrata the number of strata in the model
#' @param Xopts a [list] to overwrite defaults
#' @param X10 the initial value for X1 
#' @param X20 the initial value for X1 
#' @return a [list]
#' @export
make_Xinits_newXname = function(nStrata, Xopts = list(), X10=NULL, X20=1){with(Xopts,{
  stopifnot(is.numeric(X10))
  stopifnot(is.numeric(X20))
  X1 = checkIt(X10, nStrata)
  X2 = checkIt(X20, nStrata)
  return(list(X1=X1, X2=X2, ...))
})}

## -----------------------------------------------------------------------------
#' @title Setup Xinits.newXname
#' @description Implements [setup_Xinits] for the newXname model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.newXname = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_newXname(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the newXname model.
#' @inheritParams make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.newXname <- function(pars, i) {with(pars,{

  X1_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(X1_ix, 1)

  X2_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(X2_ix, 1)

  ... 
  
  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(X1_ix=X1_ix, X2_ix=X2_ix, ...)
  return(pars)
})}

## -----------------------------------------------------------------------------
#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xde::list_Xvars
#' @return a [list]
#' @export
list_Xvars.newXname <- function(y, pars, i) {
  with(pars$ix$X[[i]],
      return(list(
        X1 = y[X1_ix],
        X2 = y[X2_ix], 
        ...
      )
  ))
}

## -----------------------------------------------------------------------------
#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.newXname <- function(pars, i){
  with(pars$Xinits[[i]], return(c(X1,X2)))
}

## -----------------------------------------------------------------------------
#' @title Update inits for the newXname human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X 
#' @return none
#' @export
update_inits_X.newXname <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_newXname(pars, list(), X10=X1, X20=X2, ...)
    return(pars)
})}

## -----------------------------------------------------------------------------
#' @title Make parameters for newXname human model, with defaults
#' @param nStrata is the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param p1 the first parameter 
#' @param p2 the second parameter 
#' @param p3 the third parameter 
#' @return a [list]
#' @export
make_Xpar_newXname = function(nStrata, Xopts=list(),
                         p1=1, p2=2, p3=3){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("newXname")

    Xpar$p1 = checkIt(p1, nStrata)
    Xpar$p2 = checkIt(p2, nStrata)
    Xpar$p3 = checkIt(p3, nStrata)

    return(Xpar)
  })}

## -----------------------------------------------------------------------------
#' @title Setup Xpar.Xname
#' @description Implements [setup_Xpar] for the newXname model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.newXname = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_Xname(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.newXname <- function(t, y, pars, i) {
  ########################
  # extract: 
  # VAR <- y[pars$ix$X$...] 
  ########################
  with(pars$Xpar[[i]], 
       ########################
       # compute: 
       # X <- ... F(VAR) 
       ########################
  )
  return(X)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_H] for the newXname model.
#' @inheritParams ramp.xde::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.newXname <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i), {
    H <- X1 + X2 + ...
    return(H)
})}

## -----------------------------------------------------------------------------
#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the newXname model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.newXname <- function(y, pars, i) {
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

## -----------------------------------------------------------------------------
#' @title Parse the output of deSolve and return variables for the newXname model
#' @description Implements [parse_deout_X] for the newXname model
#' @inheritParams ramp.xde::parse_deout_X
#' @return none
#' @export
parse_deout_X.newXname <- function(deout, pars, i) {
  time = deout[,1]
  with(pars$ix$X[[i]],{
    X1 = deout[,X1_ix+1]
    X2 = deout[,X2_ix+1]
    ... 
    H = X1 + X2 + ...
    return(list(time=time, X1=X1, X2=X2, ..., H=H))
})}

## -----------------------------------------------------------------------------
#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the newXname model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.newXname <- function(varslist, pars) {
  pr = compute_pr_formula 
  return(pr)
}

## -----------------------------------------------------------------------------
#' @title Compute the HTC for the newXname model
#' @description Implements [HTC] for the newXname model with demography.
#' @inheritParams ramp.xde::HTC
#' @return a [numeric] vector
#' @export
HTC.newXname <- function(pars, i) {
  with(pars$Xpar[[i]],
    #HTC <- 
    return(HTC)
  )
}

## -----------------------------------------------------------------------------
#' Add lines for the density of infected individuals for the newXname model
#'
#' @param XH a list with the outputs of parse_deout_X_newXname
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
xde_lines_X_newXname = function(XH, nStrata, clrs=c("darkblue","darkred"), llty=1){
  with(XH,{
    if(nStrata==1) {
      lines(time, S, col=clrs[1], lty = llty[1])
      lines(time, I, col=clrs[2], lty = llty[1])
    }
    if(nStrata>1){
      if (length(clrs)==2) clrs=matrix(clrs, 2, nStrata)
      if (length(llty)==1) llty=rep(llty, nStrata)

      for(i in 1:nStrata){
        lines(time, S[,i], col=clrs[1,i], lty = llty[i])
        lines(time, I[,i], col=clrs[2,i], lty = llty[i])
      }
    }
  })}

## -----------------------------------------------------------------------------
#' Plot the density of infected individuals for the newXname model
#'
#' @inheritParams ramp.xde::xde_plot_X
#' @export
xde_plot_X.newXname = function(pars, i=1, clrs=c("darkblue","darkred"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "Time"))


  xde_lines_X_newXname(vars$XH[[i]], pars$Hpar[[i]]$nStrata, clrs, llty)
}

## -----------------------------------------------------------------------------
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
      derivs = c(dS,dE, dI, dR)
      
      # return the derivatives
      return(derivs)
    })
  })
}

## -----------------------------------------------------------------------------
#' @title Make initial values for the SEIRnew human model, with defaults
#' @param nStrata the number of strata in the model
#' @param Xopts a [list] to overwrite defaults
#' @param H0 the initial value for H
#' @param S0 the initial value for S
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


## -----------------------------------------------------------------------------
#' @title Setup Xinits.SEIRnew
#' @description Implements [setup_Xinits] for the SEIRnew model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SEIRnew = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_SEIRnew(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SEIRnew model.
#' @inheritParams ramp.xde::make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.SEIRnew <- function(pars, i) {with(pars,{
  
  S_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(S_ix, 1)
  
  E_ix <- seq(from = max_ix+1, length.out=nStrata)
  max_ix <- tail(E_ix, 1)
  
  I_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(I_ix, 1)
  
  R_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(R_ix, 1)
  
  
  
  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(S_ix=S_ix,  E_ix=E_ix, I_ix=I_ix, R_ix=R_ix)
  return(pars)
})}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.SEIRnew <- function(pars, i){
  with(pars$Xinits[[i]], return(c(S,E,I,R)))
}


## -----------------------------------------------------------------------------
#' @title Update inits for the SEIRnew human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.SEIRnew <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_SEIRnew(pars, list(), S0=S,E0= E, I0=I, R0=R)
    return(pars)
  })}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Setup Xpar.SEIRnew
#' @description Implements [setup_Xpar] for the SEIRnew model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SEIRnew = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SEIRnew(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIRnew model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIRnew <- function(varslist, pars,i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}

## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------

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


## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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
      derivs = c(dS,dE, dI, dR)
      
      # return the derivatives
      return(derivs)
    })
  })
}

## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
#' @title Setup Xinits.SEIRnew
#' @description Implements [setup_Xinits] for the SEIRnew model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SEIRnew = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_SEIRnew(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SEIRnew model.
#' @inheritParams ramp.xde::make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.SEIRnew <- function(pars, i) {with(pars,{
  
  S_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(S_ix, 1)
  
  E_ix <- seq(from = max_ix+1, length.out=nStrata)
  max_ix <- tail(E_ix, 1)
  
  I_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(I_ix, 1)
  
  R_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(R_ix, 1)
  
  
  
  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(S_ix=S_ix,  E_ix=E_ix, I_ix=I_ix, R_ix=R_ix)
  return(pars)
})}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.SEIRnew <- function(pars, i){
  with(pars$Xinits[[i]], return(c(S,E,I,R)))
}


## -----------------------------------------------------------------------------
#' @title Update inits for the SEIRnew human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.SEIRnew <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_SEIRnew(pars, list(), S0=S,E0= E, I0=I, R0=R)
    return(pars)
  })}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Setup Xpar.SEIRnew
#' @description Implements [setup_Xpar] for the SEIRnew model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SEIRnew = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SEIRnew(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIRnew model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIRnew <- function(varslist, pars,i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}

## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------

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


## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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
      derivs = c(dS,dE, dI, dR)
      
      # return the derivatives
      return(derivs)
    })
  })
}

## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
#' @title Setup Xinits.SEIRnew
#' @description Implements [setup_Xinits] for the SEIRnew model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SEIRnew = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_SEIRnew(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SEIRnew model.
#' @inheritParams ramp.xde::make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.SEIRnew <- function(pars, i) {with(pars,{
  
  S_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(S_ix, 1)
  
  E_ix <- seq(from = max_ix+1, length.out=nStrata)
  max_ix <- tail(E_ix, 1)
  
  I_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(I_ix, 1)
  
  R_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(R_ix, 1)
  
  
  
  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(S_ix=S_ix,  E_ix=E_ix, I_ix=I_ix, R_ix=R_ix)
  return(pars)
})}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.SEIRnew <- function(pars, i){
  with(pars$Xinits[[i]], return(c(S,E,I,R)))
}


## -----------------------------------------------------------------------------
#' @title Update inits for the SEIRnew human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.SEIRnew <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_SEIRnew(pars, list(), S0=S,E0= E, I0=I, R0=R)
    return(pars)
  })}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Setup Xpar.SEIRnew
#' @description Implements [setup_Xpar] for the SEIRnew model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SEIRnew = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SEIRnew(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIRnew model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIRnew <- function(varslist, pars,i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}

## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------

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


## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
#' @title Setup Xinits.SEIRnew
#' @description Implements [setup_Xinits] for the SEIRnew model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SEIRnew = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_SEIRnew(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.SEIRnew <- function(pars, i){
  with(pars$Xinits[[i]], return(c(S,E,I,R)))
}


## -----------------------------------------------------------------------------
#' @title Update inits for the SEIRnew human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.SEIRnew <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_SEIRnew(pars, list(), S0=S,E0= E, I0=I, R0=R)
    return(pars)
  })}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Setup Xpar.SEIRnew
#' @description Implements [setup_Xpar] for the SEIRnew model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SEIRnew = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SEIRnew(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIRnew model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIRnew <- function(varslist, pars,i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}

## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------

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


## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [dXdt] for the SEIRV model
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SEIRV<- function(t, y, pars, i) {
  # do not change this
  foi <- pars$FoI[[i]]
  
  # attach the variables by name
  with(list_Xvars(y, pars, i),{
    # compute H (if it isn't one of the variables)
    H <- F_H(t, y, pars, i)
    
    # expose the parameters (see make_Xpar_SEIRV)
    with(pars$Xpar[[i]], {
      # compute the derivatives
      dS <- (1-tau)*Births(t, H, pars,i) - foi*S + dHdt(t, S, pars,i)+ epsilon*R
      dE <- foi*S - alpha*E + dHdt(t, E, pars,i)
      dI <- alpha*E - r*I + dHdt(t, I, pars,i)
      dR <- (1-delta)*r*I -sig*R -epislon*R+ dHdt(t, R, pars,i)
      dV <- tau**Births(t, H, pars,i) + delta*r*I + sig*R + dHdt(t, V, pars,i)
      
       # concatenate the derivatives
      derivs = c(dS, dE, dI, dR)
      
      # return the derivatives
      return(derivs)
    })
  })
}

## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
#' @title Setup Xinits.SEIRV
#' @description Implements [setup_Xinits] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SEIRV = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_SEIRV(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SEIRV model.
#' @inheritParams ramp.xde::make_indices_X
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

## -----------------------------------------------------------------------------
#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xde::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SEIRV <- function(y, pars, i) {
  with(pars$ix$X[[i]],
       return(list(
         S = y[S_ix],
         E = y[E_ix],
         I = y[I_ix],
         R = y[R_ix],
         V = y[V_ix]
       )
       ))
}

## -----------------------------------------------------------------------------
#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.SEIRV <- function(pars, i){
  with(pars$Xinits[[i]], return(c(S,E,I,R,V)))
}


## -----------------------------------------------------------------------------
#' @title Update inits for the SEIRV human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.SEIRV <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_SEIRV(pars, list(), S0=S,E0= E, I0=I, R0=R,V0 =V)
    return(pars)
  })}

## -----------------------------------------------------------------------------
#' @title Make parameters for SEIRV human model, with defaults
#' @param nStrata is the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param alpha  incubation rate
#' @param b the proportion of infective bites that cause an infection
#' @param r the the duration of an infection
#' @param c the proportion of bites on infected humans that infect a mosquito
#' @param tau  travallers with immunity from malaria
#' @param epsilon   loss of immunity rate
#' @param sig  progression rate of recovered individual to protected class
#' @return a [list]
#' @export
make_Xpar_SEIRV = function(nStrata, Xopts=list(),
                          alpha =0.1,b=0.55, r=1/180, c=0.15,tau= 0.5,epsilon,sig=0.3){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SEIRV")
    
    Xpar$alpha = checkIt(alpha, nStrata)
    Xpar$b = checkIt(b, nStrata)
    Xpar$c = checkIt(c, nStrata)
    Xpar$r = checkIt(r, nStrata)
    Xpar$tau = checkIt(tau, nStrata)
    Xpar$epsilon = checkIt(epsilon, nStrata)
    Xpar$sig = checkIt(sig, nStrata)
    
    return(Xpar)
  })}

## -----------------------------------------------------------------------------
#' @title Setup Xpar.SEIRV
#' @description Implements [setup_Xpar] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SEIRV = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SEIRV(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SEIRV <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  Y = with(pars$Xpar[[i]], c*I)
  return(Y)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SEIRV model.
#' @inheritParams ramp.xde::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SEIRV <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i), {
    H <- S +E+ I+R+V
    return(H)
  })}

## -----------------------------------------------------------------------------
#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SEIRV model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SEIRV <- function(y, pars, i) {
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

## -----------------------------------------------------------------------------
#' @title Parse the output of deSolve and return variables for the SEIRV model
#' @description Implements [parse_deout_X] for the SEIRV model
#' @inheritParams ramp.xde::parse_deout_X
#' @return none
#' @export
parse_deout_X.SEIRV <- function(deout, pars, i) {
  time = deout[,1]
  with(pars$ix$X[[i]],{
    S = deout[,S_ix+1]
    E = deout[,E_ix+1]
    I = deout[,I_ix+1]
    R = deout[,R_ix+1]
    V = deout[,V_ix+1]
    H = S+E+I+R+V
    return(list(time=time, S=S, E=E,I=I, R=R, V=V,H=H))
  })}

## -----------------------------------------------------------------------------
#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIRV model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIRV <- function(varslist, pars,i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}

## -----------------------------------------------------------------------------
#' @title Compute the HTC for the SEIRV model
#' @description Implements [HTC] for the SEIRV model with demography.
#' @inheritParams ramp.xde::HTC
#' @return a [numeric] vector
#' @export
HTC.SEIRV <- function(pars, i) {
  with(pars$Xpar[[i]],
       HTC <- c/r,
       return(HTC)
  )
}


## -----------------------------------------------------------------------------

#' Add lines for the density of infected individuals for the SEIRV model
#'
#' @param XH a list with the outputs of parse_deout_X_SEIRV
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


## -----------------------------------------------------------------------------
#' Plot the density of infected individuals for the SEIRV model
#'
#' @inheritParams ramp.xde::xde_plot_X
#' @export
xde_plot_X.SEIRV = function(pars, i=1, clrs=c("black","darkblue","darkred","darkgreen","purple"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})
  
  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "No of. Infected", xlab = "Time"))
  
  
  xde_lines_X_SEIRV(vars$XH[[i]], pars$Hpar[[i]]$nStrata, clrs, llty)
}

## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [dXdt] for the SEIRV model
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SEIRV<- function(t, y, pars, i) {
  # do not change this
  foi <- pars$FoI[[i]]
  
  # attach the variables by name
  with(list_Xvars(y, pars, i),{
    # compute H (if it isn't one of the variables)
    H <- F_H(t, y, pars, i)
    
    # expose the parameters (see make_Xpar_SEIRV)
    with(pars$Xpar[[i]], {
      # compute the derivatives
      dS <- (1-alpha)*Births(t, H, pars,i) - foi*S + dHdt(t, S, pars,i)+ epsilon*R
      dE <- foi*S - tau*E + dHdt(t, E, pars,i)
      dI <- tau*E - r*I + dHdt(t, I, pars,i)
      dR <- (1-delta)*r*I -sig*R -epislon*R+ dHdt(t, R, pars,i)
      dV <- alpha*Births(t, H, pars,i) + delta*r*I + sig*R + dHdt(t, V, pars,i)
      
       # concatenate the derivatives
      derivs = c(dS, dE, dI, dR,dV)
      
      # return the derivatives
      return(derivs)
    })
  })
}

## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
#' @title Setup Xinits.SEIRV
#' @description Implements [setup_Xinits] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SEIRV = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_SEIRV(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SEIRV model.
#' @inheritParams ramp.xde::make_indices_X
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

## -----------------------------------------------------------------------------
#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xde::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SEIRV <- function(y, pars, i) {
  with(pars$ix$X[[i]],
       return(list(
         S = y[S_ix],
         E = y[E_ix],
         I = y[I_ix],
         R = y[R_ix],
         V = y[V_ix]
       )
       ))
}

## -----------------------------------------------------------------------------
#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.SEIRV <- function(pars, i){
  with(pars$Xinits[[i]], return(c(S,E,I,R,V)))
}


## -----------------------------------------------------------------------------
#' @title Update inits for the SEIRV human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.SEIRV <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_SEIRV(pars, list(), S0=S,E0= E, I0=I, R0=R,V0 =V)
    return(pars)
  })}

## -----------------------------------------------------------------------------
#' @title Make parameters for SEIRV human model, with defaults
#' @param nStrata is the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param tau  incubation rate
#' @param b the proportion of infective bites that cause an infection
#' @param r the the duration of an infection
#' @param c the proportion of bites on infected humans that infect a mosquito
#' @param alpha  travallers with immunity from malaria
#' @param epsilon   loss of immunity rate
#' @param delta   proportion of protected after recovery
#' @param sig  progression rate of recovered individual to protected class
#' @return a [list]
#' @export
make_Xpar_SEIRV = function(nStrata, Xopts=list(),
                          alpha =0.1,b=0.55, r=1/180, c=0.15,tau= 0.5,epsilon,sig=0.3,delta =0.5){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SEIRV")
    
    Xpar$alpha = checkIt(alpha, nStrata)
    Xpar$b = checkIt(b, nStrata)
    Xpar$c = checkIt(c, nStrata)
    Xpar$r = checkIt(r, nStrata)
    Xpar$tau = checkIt(tau, nStrata)
    Xpar$epsilon = checkIt(epsilon, nStrata)
    Xpar$sig = checkIt(sig, nStrata)
    Xpar$delta = checkIt(delta, nStrata)
    
    return(Xpar)
  })}

## -----------------------------------------------------------------------------
#' @title Setup Xpar.SEIRV
#' @description Implements [setup_Xpar] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SEIRV = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SEIRV(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SEIRV <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  Y = with(pars$Xpar[[i]], c*I)
  return(Y)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SEIRV model.
#' @inheritParams ramp.xde::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SEIRV <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i), {
    H <- S +E+ I+R+V
    return(H)
  })}

## -----------------------------------------------------------------------------
#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SEIRV model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SEIRV <- function(y, pars, i) {
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

## -----------------------------------------------------------------------------
#' @title Parse the output of deSolve and return variables for the SEIRV model
#' @description Implements [parse_deout_X] for the SEIRV model
#' @inheritParams ramp.xde::parse_deout_X
#' @return none
#' @export
parse_deout_X.SEIRV <- function(deout, pars, i) {
  time = deout[,1]
  with(pars$ix$X[[i]],{
    S = deout[,S_ix+1]
    E = deout[,E_ix+1]
    I = deout[,I_ix+1]
    R = deout[,R_ix+1]
    V = deout[,V_ix+1]
    H = S+E+I+R+V
    return(list(time=time, S=S, E=E,I=I, R=R, V=V,H=H))
  })}

## -----------------------------------------------------------------------------
#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIRV model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIRV <- function(varslist, pars,i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}

## -----------------------------------------------------------------------------
#' @title Compute the HTC for the SEIRV model
#' @description Implements [HTC] for the SEIRV model with demography.
#' @inheritParams ramp.xde::HTC
#' @return a [numeric] vector
#' @export
HTC.SEIRV <- function(pars, i) {
  with(pars$Xpar[[i]],
       HTC <- c/r,
       return(HTC)
  )
}


## -----------------------------------------------------------------------------

#' Add lines for the density of infected individuals for the SEIRV model
#'
#' @param XH a list with the outputs of parse_deout_X_SEIRV
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


## -----------------------------------------------------------------------------
#' Plot the density of infected individuals for the SEIRV model
#'
#' @inheritParams ramp.xde::xde_plot_X
#' @export
xde_plot_X.SEIRV = function(pars, i=1, clrs=c("black","darkblue","darkred","darkgreen","purple"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})
  
  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "No of. Infected", xlab = "Time"))
  
  
  xde_lines_X_SEIRV(vars$XH[[i]], pars$Hpar[[i]]$nStrata, clrs, llty)
}

## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [dXdt] for the SEIRV model
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SEIRV<- function(t, y, pars, i) {
  # do not change this
  foi <- pars$FoI[[i]]
  
  # attach the variables by name
  with(list_Xvars(y, pars, i),{
    # compute H (if it isn't one of the variables)
    H <- F_H(t, y, pars, i)
    
    # expose the parameters (see make_Xpar_SEIRV)
    with(pars$Xpar[[i]], {
      # compute the derivatives
      dS <- (1-alpha)*Births(t, H, pars,i) - foi*S + dHdt(t, S, pars,i)+ eps*R
      dE <- foi*S - tau*E + dHdt(t, E, pars,i)
      dI <- tau*E - r*I + dHdt(t, I, pars,i)
      dR <- (1-delta)*r*I -sig*R -eps*R+ dHdt(t, R, pars,i)
      dV <- alpha*Births(t, H, pars,i) + delta*r*I + sig*R + dHdt(t, V, pars,i)
      
       # concatenate the derivatives
      derivs = c(dS, dE, dI, dR,dV)
      
      # return the derivatives
      return(derivs)
    })
  })
}

## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
#' @title Setup Xinits.SEIRV
#' @description Implements [setup_Xinits] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SEIRV = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_SEIRV(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SEIRV model.
#' @inheritParams ramp.xde::make_indices_X
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

## -----------------------------------------------------------------------------
#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xde::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SEIRV <- function(y, pars, i) {
  with(pars$ix$X[[i]],
       return(list(
         S = y[S_ix],
         E = y[E_ix],
         I = y[I_ix],
         R = y[R_ix],
         V = y[V_ix]
       )
       ))
}

## -----------------------------------------------------------------------------
#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.SEIRV <- function(pars, i){
  with(pars$Xinits[[i]], return(c(S,E,I,R,V)))
}


## -----------------------------------------------------------------------------
#' @title Update inits for the SEIRV human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.SEIRV <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_SEIRV(pars, list(), S0=S,E0= E, I0=I, R0=R,V0 =V)
    return(pars)
  })}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Setup Xpar.SEIRV
#' @description Implements [setup_Xpar] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SEIRV = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SEIRV(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SEIRV <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  Y = with(pars$Xpar[[i]], c*I)
  return(Y)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SEIRV model.
#' @inheritParams ramp.xde::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SEIRV <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i), {
    H <- S +E+ I+R+V
    return(H)
  })}

## -----------------------------------------------------------------------------
#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SEIRV model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SEIRV <- function(y, pars, i) {
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

## -----------------------------------------------------------------------------
#' @title Parse the output of deSolve and return variables for the SEIRV model
#' @description Implements [parse_deout_X] for the SEIRV model
#' @inheritParams ramp.xde::parse_deout_X
#' @return none
#' @export
parse_deout_X.SEIRV <- function(deout, pars, i) {
  time = deout[,1]
  with(pars$ix$X[[i]],{
    S = deout[,S_ix+1]
    E = deout[,E_ix+1]
    I = deout[,I_ix+1]
    R = deout[,R_ix+1]
    V = deout[,V_ix+1]
    H = S+E+I+R+V
    return(list(time=time, S=S, E=E,I=I, R=R, V=V,H=H))
  })}

## -----------------------------------------------------------------------------
#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIRV model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIRV <- function(varslist, pars,i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}

## -----------------------------------------------------------------------------
#' @title Compute the HTC for the SEIRV model
#' @description Implements [HTC] for the SEIRV model with demography.
#' @inheritParams ramp.xde::HTC
#' @return a [numeric] vector
#' @export
HTC.SEIRV <- function(pars, i) {
  with(pars$Xpar[[i]],
       HTC <- c/r,
       return(HTC)
  )
}


## -----------------------------------------------------------------------------

#' Add lines for the density of infected individuals for the SEIRV model
#'
#' @param XH a list with the outputs of parse_deout_X_SEIRV
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


## -----------------------------------------------------------------------------
#' Plot the density of infected individuals for the SEIRV model
#'
#' @inheritParams ramp.xde::xde_plot_X
#' @export
xde_plot_X.SEIRV = function(pars, i=1, clrs=c("black","darkblue","darkred","darkgreen","purple"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})
  
  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "No of. Infected", xlab = "Time"))
  
  
  xde_lines_X_SEIRV(vars$XH[[i]], pars$Hpar[[i]]$nStrata, clrs, llty)
}

## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [dXdt] for the SEIRV model
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SEIRV<- function(t, y, pars, i) {
  # do not change this
  foi <- pars$FoI[[i]]
  
  # attach the variables by name
  with(list_Xvars(y, pars, i),{
    # compute H (if it isn't one of the variables)
    H <- F_H(t, y, pars, i)
    
    # expose the parameters (see make_Xpar_SEIRV)
    with(pars$Xpar[[i]], {
      # compute the derivatives
      dS <- (1-alpha)*Births(t, H, pars,i) - foi*S + dHdt(t, S, pars,i)+ eps*R
      dE <- foi*S - tau*E + dHdt(t, E, pars,i)
      dI <- tau*E - r*I + dHdt(t, I, pars,i)
      dR <- (1-delta)*r*I -sig*R -eps*R+ dHdt(t, R, pars,i)
      dV <- alpha*Births(t, H, pars,i) + delta*r*I + sig*R + dHdt(t, V, pars,i)
      
       # concatenate the derivatives
      derivs = c(dS, dE, dI, dR,dV)
      
      # return the derivatives
      return(derivs)
    })
  })
}

## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
#' @title Setup Xinits.SEIRV
#' @description Implements [setup_Xinits] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SEIRV = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_SEIRV(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SEIRV model.
#' @inheritParams ramp.xde::make_indices_X
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

## -----------------------------------------------------------------------------
#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xde::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SEIRV <- function(y, pars, i) {
  with(pars$ix$X[[i]],
       return(list(
         S = y[S_ix],
         E = y[E_ix],
         I = y[I_ix],
         R = y[R_ix],
         V = y[V_ix]
       )
       ))
}

## -----------------------------------------------------------------------------
#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.SEIRV <- function(pars, i){
  with(pars$Xinits[[i]], return(c(S,E,I,R,V)))
}


## -----------------------------------------------------------------------------
#' @title Update inits for the SEIRV human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.SEIRV <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_SEIRV(pars, list(), S0=S,E0= E, I0=I, R0=R,V0 =V)
    return(pars)
  })}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Setup Xpar.SEIRV
#' @description Implements [setup_Xpar] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SEIRV = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SEIRV(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SEIRV <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  Y = with(pars$Xpar[[i]], c*I)
  return(Y)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SEIRV model.
#' @inheritParams ramp.xde::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SEIRV <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i), {
    H <- S +E+ I+R+V
    return(H)
  })}

## -----------------------------------------------------------------------------
#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SEIRV model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SEIRV <- function(y, pars, i) {
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

## -----------------------------------------------------------------------------
#' @title Parse the output of deSolve and return variables for the SEIRV model
#' @description Implements [parse_deout_X] for the SEIRV model
#' @inheritParams ramp.xde::parse_deout_X
#' @return none
#' @export
parse_deout_X.SEIRV <- function(deout, pars, i) {
  time = deout[,1]
  with(pars$ix$X[[i]],{
    S = deout[,S_ix+1]
    E = deout[,E_ix+1]
    I = deout[,I_ix+1]
    R = deout[,R_ix+1]
    V = deout[,V_ix+1]
    H = S+E+I+R+V
    return(list(time=time, S=S, E=E,I=I, R=R, V=V,H=H))
  })}

## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [dXdt] for the SEIRV model
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SEIRV<- function(t, y, pars, i) {
  # do not change this
  foi <- pars$FoI[[i]]
  
  # attach the variables by name
  with(list_Xvars(y, pars, i),{
    # compute H (if it isn't one of the variables)
    H <- F_H(t, y, pars, i)
    
    # expose the parameters (see make_Xpar_SEIRV)
    with(pars$Xpar[[i]], {
      # compute the derivatives
      dS <- (1-alpha)*Births(t, H, pars,i) - foi*S + dHdt(t, S, pars,i)+ eps*R
      dE <- foi*S - tau*E + dHdt(t, E, pars,i)
      dI <- tau*E - r*I + dHdt(t, I, pars,i)
      dR <- (1-delta)*r*I -sig*R -eps*R+ dHdt(t, R, pars,i)
      dV <- alpha*Births(t, H, pars,i) + delta*r*I + sig*R + dHdt(t, V, pars,i)
      
       # concatenate the derivatives
      derivs = c(dS, dE, dI, dR,dV)
      
      # return the derivatives
      return(derivs)
    })
  })
}

## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
#' @title Setup Xinits.SEIRV
#' @description Implements [setup_Xinits] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SEIRV = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_SEIRV(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SEIRV model.
#' @inheritParams ramp.xde::make_indices_X
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

## -----------------------------------------------------------------------------
#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xde::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SEIRV <- function(y, pars, i) {
  with(pars$ix$X[[i]],
       return(list(
         S = y[S_ix],
         E = y[E_ix],
         I = y[I_ix],
         R = y[R_ix],
         V = y[V_ix]
       )
       ))
}

## -----------------------------------------------------------------------------
#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.SEIRV <- function(pars, i){
  with(pars$Xinits[[i]], return(c(S,E,I,R,V)))
}


## -----------------------------------------------------------------------------
#' @title Update inits for the SEIRV human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.SEIRV <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_SEIRV(pars, list(), S0=S,E0= E, I0=I, R0=R,V0 =V)
    return(pars)
  })}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Setup Xpar.SEIRV
#' @description Implements [setup_Xpar] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SEIRV = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SEIRV(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SEIRV <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  Y = with(pars$Xpar[[i]], c*I)
  return(Y)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SEIRV model.
#' @inheritParams ramp.xde::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SEIRV <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i), {
    H <- S +E+ I+R+V
    return(H)
  })}

## -----------------------------------------------------------------------------
#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SEIRV model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SEIRV <- function(y, pars, i) {
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

## -----------------------------------------------------------------------------
#' @title Parse the output of deSolve and return variables for the SEIRV model
#' @description Implements [parse_deout_X] for the SEIRV model
#' @inheritParams ramp.xde::parse_deout_X
#' @return none
#' @export
parse_deout_X.SEIRV <- function(deout, pars, i) {
  time = deout[,1]
  with(pars$ix$X[[i]],{
    S = deout[,S_ix+1]
    E = deout[,E_ix+1]
    I = deout[,I_ix+1]
    R = deout[,R_ix+1]
    V = deout[,V_ix+1]
    H = S+E+I+R+V
    return(list(time=time, S=S, E=E,I=I, R=R, V=V,H=H))
  })}

## -----------------------------------------------------------------------------
#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIRV model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIRV <- function(varslist, pars,i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}

## -----------------------------------------------------------------------------
#' @title Compute the HTC for the SEIRV model
#' @description Implements [HTC] for the SEIRV model with demography.
#' @inheritParams ramp.xde::HTC
#' @return a [numeric] vector
#' @export
HTC.SEIRV <- function(pars, i) {
  with(pars$Xpar[[i]],
       HTC <- c/r,
       return(HTC)
  )
}


## -----------------------------------------------------------------------------

#' Add lines for the density of infected individuals for the SEIRV model
#'
#' @param XH a list with the outputs of parse_deout_X_SEIRV
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


## -----------------------------------------------------------------------------
#' Plot the density of infected individuals for the SEIRV model
#'
#' @inheritParams ramp.xde::xde_plot_X
#' @export
xde_plot_X.SEIRV = function(pars, i=1, clrs=c("black","darkblue","darkred","darkgreen","purple"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})
  
  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "No of. Infected", xlab = "Time"))
  
  
  xde_lines_X_SEIRV(vars$XH[[i]], pars$Hpar[[i]]$nStrata, clrs, llty)
}

## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [dXdt] for the SEIRV model
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SEIRV<- function(t, y, pars, i) {
  # do not change this
  foi <- pars$FoI[[i]]
  
  # attach the variables by name
  with(list_Xvars(y, pars, i),{
    # compute H (if it isn't one of the variables)
    H <- F_H(t, y, pars, i)
    
    # expose the parameters (see make_Xpar_SEIRV)
    with(pars$Xpar[[i]], {
      # compute the derivatives
      dS <- (1-alpha)*Births(t, H, pars,i) - foi*S + dHdt(t, S, pars,i)+ eps*R
      dE <- foi*S - tau*E + dHdt(t, E, pars,i)
      dI <- tau*E - r*I + dHdt(t, I, pars,i)
      dR <- (1-delta)*r*I -sig*R -eps*R+ dHdt(t, R, pars,i)
      dV <- alpha*Births(t, H, pars,i) + delta*r*I + sig*R + dHdt(t, V, pars,i)
      
       # concatenate the derivatives
      derivs = c(dS, dE, dI, dR,dV)
      
      # return the derivatives
      return(derivs)
    })
  })
}

## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
#' @title Setup Xinits.SEIRV
#' @description Implements [setup_Xinits] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SEIRV = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_SEIRV(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SEIRV model.
#' @inheritParams ramp.xde::make_indices_X
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

## -----------------------------------------------------------------------------
#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xde::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SEIRV <- function(y, pars, i) {
  with(pars$ix$X[[i]],
       return(list(
         S = y[S_ix],
         E = y[E_ix],
         I = y[I_ix],
         R = y[R_ix],
         V = y[V_ix]
       )
       ))
}

## -----------------------------------------------------------------------------
#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.SEIRV <- function(pars, i){
  with(pars$Xinits[[i]], return(c(S,E,I,R,V)))
}


## -----------------------------------------------------------------------------
#' @title Update inits for the SEIRV human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.SEIRV <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_SEIRV(pars, list(), S0=S,E0= E, I0=I, R0=R,V0 =V)
    return(pars)
  })}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Setup Xpar.SEIRV
#' @description Implements [setup_Xpar] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SEIRV = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SEIRV(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SEIRV <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  Y = with(pars$Xpar[[i]], c*I)
  return(Y)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SEIRV model.
#' @inheritParams ramp.xde::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SEIRV <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i), {
    H <- S +E+ I+R+V
    return(H)
  })}

## -----------------------------------------------------------------------------
#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SEIRV model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SEIRV <- function(y, pars, i) {
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

## -----------------------------------------------------------------------------
#' @title Parse the output of deSolve and return variables for the SEIRV model
#' @description Implements [parse_deout_X] for the SEIRV model
#' @inheritParams ramp.xde::parse_deout_X
#' @return none
#' @export
parse_deout_X.SEIRV <- function(deout, pars, i) {
  time = deout[,1]
  with(pars$ix$X[[i]],{
    S = deout[,S_ix+1]
    E = deout[,E_ix+1]
    I = deout[,I_ix+1]
    R = deout[,R_ix+1]
    V = deout[,V_ix+1]
    H = S+E+I+R+V
    return(list(time=time, S=S, E=E,I=I, R=R, V=V,H=H))
  })}

## -----------------------------------------------------------------------------
#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIRV model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIRV <- function(varslist, pars,i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}

## -----------------------------------------------------------------------------
#' @title Compute the HTC for the SEIRV model
#' @description Implements [HTC] for the SEIRV model with demography.
#' @inheritParams ramp.xde::HTC
#' @return a [numeric] vector
#' @export
HTC.SEIRV <- function(pars, i) {
  with(pars$Xpar[[i]],
       HTC <- c/r,
       return(HTC)
  )
}


## -----------------------------------------------------------------------------

#' Add lines for the density of infected individuals for the SEIRV model
#'
#' @param XH a list with the outputs of parse_deout_X_SEIRV
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


## -----------------------------------------------------------------------------
#' Plot the density of infected individuals for the SEIRV model
#'
#' @inheritParams ramp.xde::xde_plot_X
#' @export
xde_plot_X.SEIRV = function(pars, i=1, clrs=c("black","darkblue","darkred","darkgreen","purple"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})
  
  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "No of. Infected", xlab = "Time"))
  
  
  xde_lines_X_SEIRV(vars$XH[[i]], pars$Hpar[[i]]$nStrata, clrs, llty)
}

## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [dXdt] for the SEIRV model
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SEIRV<- function(t, y, pars, i) {
  # do not change this
  foi <- pars$FoI[[i]]
  
  # attach the variables by name
  with(list_Xvars(y, pars, i),{
    # compute H (if it isn't one of the variables)
    H <- F_H(t, y, pars, i)
    
    # expose the parameters (see make_Xpar_SEIRV)
    with(pars$Xpar[[i]], {
      # compute the derivatives
      dS <- (1-alpha)*Births(t, H, pars,i) - foi*S + dHdt(t, S, pars,i)+ eps*R
      dE <- foi*S - tau*E + dHdt(t, E, pars,i)
      dI <- tau*E - r*I + dHdt(t, I, pars,i)
      dR <- (1-delta)*r*I -sig*R -eps*R+ dHdt(t, R, pars,i)
      dV <- alpha*Births(t, H, pars,i) + delta*r*I + sig*R + dHdt(t, V, pars,i)
      
       # concatenate the derivatives
      derivs = c(dS, dE, dI, dR,dV)
      
      # return the derivatives
      return(derivs)
    })
  })
}

## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
#' @title Setup Xinits.SEIRV
#' @description Implements [setup_Xinits] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SEIRV = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_SEIRV(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SEIRV model.
#' @inheritParams ramp.xde::make_indices_X
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

## -----------------------------------------------------------------------------
#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xde::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SEIRV <- function(y, pars, i) {
  with(pars$ix$X[[i]],
       return(list(
         S = y[S_ix],
         E = y[E_ix],
         I = y[I_ix],
         R = y[R_ix],
         V = y[V_ix]
       )
       ))
}

## -----------------------------------------------------------------------------
#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.SEIRV <- function(pars, i){
  with(pars$Xinits[[i]], return(c(S,E,I,R,V)))
}


## -----------------------------------------------------------------------------
#' @title Update inits for the SEIRV human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.SEIRV <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_SEIRV(pars, list(), S0=S,E0= E, I0=I, R0=R,V0 =V)
    return(pars)
  })}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Setup Xpar.SEIRV
#' @description Implements [setup_Xpar] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SEIRV = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SEIRV(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SEIRV <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  Y = with(pars$Xpar[[i]], c*I)
  return(Y)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SEIRV model.
#' @inheritParams ramp.xde::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SEIRV <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i), {
    H <- S +E+ I+R+V
    return(H)
  })}

## -----------------------------------------------------------------------------
#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SEIRV model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SEIRV <- function(y, pars, i) {
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

## -----------------------------------------------------------------------------
#' @title Parse the output of deSolve and return variables for the SEIRV model
#' @description Implements [parse_deout_X] for the SEIRV model
#' @inheritParams ramp.xde::parse_deout_X
#' @return none
#' @export
parse_deout_X.SEIRV <- function(deout, pars, i) {
  time = deout[,1]
  with(pars$ix$X[[i]],{
    S = deout[,S_ix+1]
    E = deout[,E_ix+1]
    I = deout[,I_ix+1]
    R = deout[,R_ix+1]
    V = deout[,V_ix+1]
    H = S+E+I+R+V
    return(list(time=time, S=S, E=E,I=I, R=R, V=V,H=H))
  })}

## -----------------------------------------------------------------------------
#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIRV model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIRV <- function(varslist, pars,i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}

## -----------------------------------------------------------------------------
#' @title Compute the HTC for the SEIRV model
#' @description Implements [HTC] for the SEIRV model with demography.
#' @inheritParams ramp.xde::HTC
#' @return a [numeric] vector
#' @export
HTC.SEIRV <- function(pars, i) {
  with(pars$Xpar[[i]],
       HTC <- c/r,
       return(HTC)
  )
}


## -----------------------------------------------------------------------------

#' Add lines for the density of infected individuals for the SEIRV model
#'
#' @param XH a list with the outputs of parse_deout_X_SEIRV
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


## -----------------------------------------------------------------------------
#' Plot the density of infected individuals for the SEIRV model
#'
#' @inheritParams ramp.xde::xde_plot_X
#' @export
xde_plot_X.SEIRV = function(pars, i=1, clrs=c("black","darkblue","darkred","darkgreen","purple"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})
  
  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "No of. Infected", xlab = "Time"))
  
  
  xde_lines_X_SEIRV(vars$XH[[i]], pars$Hpar[[i]]$nStrata, clrs, llty)
}

## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [dXdt] for the SEIRV model
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SEIRV<- function(t, y, pars, i) {
  # do not change this
  foi <- pars$FoI[[i]]
  
  # attach the variables by name
  with(list_Xvars(y, pars, i),{
    # compute H (if it isn't one of the variables)
    H <- F_H(t, y, pars, i)
    
    # expose the parameters (see make_Xpar_SEIRV)
    with(pars$Xpar[[i]], {
      # compute the derivatives
      dS <- (1-alpha)*Births(t, H, pars,i) - foi*S + dHdt(t, S, pars,i)+ eps*R
      dE <- foi*S - tau*E + dHdt(t, E, pars,i)
      dI <- tau*E - r*I + dHdt(t, I, pars,i)
      dR <- (1-delta)*r*I -sig*R -eps*R+ dHdt(t, R, pars,i)
      dV <- alpha*Births(t, H, pars,i) + delta*r*I + sig*R + dHdt(t, V, pars,i)
      
       # concatenate the derivatives
      derivs = c(dS, dE, dI, dR,dV)
      
      # return the derivatives
      return(derivs)
    })
  })
}

## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
#' @title Setup Xinits.SEIRV
#' @description Implements [setup_Xinits] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SEIRV = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_SEIRV(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SEIRV model.
#' @inheritParams ramp.xde::make_indices_X
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

## -----------------------------------------------------------------------------
#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xde::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SEIRV <- function(y, pars, i) {
  with(pars$ix$X[[i]],
       return(list(
         S = y[S_ix],
         E = y[E_ix],
         I = y[I_ix],
         R = y[R_ix],
         V = y[V_ix]
       )
       ))
}

## -----------------------------------------------------------------------------
#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.SEIRV <- function(pars, i){
  with(pars$Xinits[[i]], return(c(S,E,I,R,V)))
}


## -----------------------------------------------------------------------------
#' @title Update inits for the SEIRV human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.SEIRV <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_SEIRV(pars, list(), S0=S,E0= E, I0=I, R0=R,V0 =V)
    return(pars)
  })}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Setup Xpar.SEIRV
#' @description Implements [setup_Xpar] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SEIRV = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SEIRV(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SEIRV <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  Y = with(pars$Xpar[[i]], c*I)
  return(Y)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SEIRV model.
#' @inheritParams ramp.xde::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SEIRV <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i), {
    H <- S +E+ I+R+V
    return(H)
  })}

## -----------------------------------------------------------------------------
#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SEIRV model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SEIRV <- function(y, pars, i) {
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

## -----------------------------------------------------------------------------
#' @title Parse the output of deSolve and return variables for the SEIRV model
#' @description Implements [parse_deout_X] for the SEIRV model
#' @inheritParams ramp.xde::parse_deout_X
#' @return none
#' @export
parse_deout_X.SEIRV <- function(deout, pars, i) {
  time = deout[,1]
  with(pars$ix$X[[i]],{
    S = deout[,S_ix+1]
    E = deout[,E_ix+1]
    I = deout[,I_ix+1]
    R = deout[,R_ix+1]
    V = deout[,V_ix+1]
    H = S+E+I+R+V
    return(list(time=time, S=S, E=E,I=I, R=R, V=V,H=H))
  })}

## -----------------------------------------------------------------------------
#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIRV model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIRV <- function(varslist, pars,i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}

## -----------------------------------------------------------------------------
#' @title Compute the HTC for the SEIRV model
#' @description Implements [HTC] for the SEIRV model with demography.
#' @inheritParams ramp.xde::HTC
#' @return a [numeric] vector
#' @export
HTC.SEIRV <- function(pars, i) {
  with(pars$Xpar[[i]],
       HTC <- c/r,
       return(HTC)
  )
}


## -----------------------------------------------------------------------------

#' Add lines for the density of infected individuals for the SEIRV model
#'
#' @param XH a list with the outputs of parse_deout_X_SEIRV
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


## -----------------------------------------------------------------------------
#' Plot the density of infected individuals for the SEIRV model
#'
#' @inheritParams ramp.xde::xde_plot_X
#' @export
xde_plot_X.SEIRV = function(pars, i=1, clrs=c("black","darkblue","darkred","darkgreen","purple"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})
  
  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "No of. Infected", xlab = "Time"))
  
  
  xde_lines_X_SEIRV(vars$XH[[i]], pars$Hpar[[i]]$nStrata, clrs, llty)
}

## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [dXdt] for the SEIRV model
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SEIRV<- function(t, y, pars, i) {
  # do not change this
  foi <- pars$FoI[[i]]
  
  # attach the variables by name
  with(list_Xvars(y, pars, i),{
    # compute H (if it isn't one of the variables)
    H <- F_H(t, y, pars, i)
    
    # expose the parameters (see make_Xpar_SEIRV)
    with(pars$Xpar[[i]], {
      # compute the derivatives
      dS <- (1-alpha)*Births(t, H, pars,i) - foi*S + dHdt(t, S, pars,i)+ eps*R
      dE <- foi*S - tau*E + dHdt(t, E, pars,i)
      dI <- tau*E - r*I + dHdt(t, I, pars,i)
      dR <- (1-delta)*r*I -sig*R -eps*R+ dHdt(t, R, pars,i)
      dV <- alpha*Births(t, H, pars,i) + delta*r*I + sig*R + dHdt(t, V, pars,i)
      
       # concatenate the derivatives
      derivs = c(dS, dE, dI, dR,dV)
      
      # return the derivatives
      return(derivs)
    })
  })
}

## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
#' @title Setup Xinits.SEIRV
#' @description Implements [setup_Xinits] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SEIRV = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_SEIRV(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SEIRV model.
#' @inheritParams ramp.xde::make_indices_X
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

## -----------------------------------------------------------------------------
#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xde::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SEIRV <- function(y, pars, i) {
  with(pars$ix$X[[i]],
       return(list(
         S = y[S_ix],
         E = y[E_ix],
         I = y[I_ix],
         R = y[R_ix],
         V = y[V_ix]
       )
       ))
}

## -----------------------------------------------------------------------------
#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.SEIRV <- function(pars, i){
  with(pars$Xinits[[i]], return(c(S,E,I,R,V)))
}


## -----------------------------------------------------------------------------
#' @title Update inits for the SEIRV human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.SEIRV <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_SEIRV(pars, list(), S0=S,E0= E, I0=I, R0=R,V0 =V)
    return(pars)
  })}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Setup Xpar.SEIRV
#' @description Implements [setup_Xpar] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SEIRV = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SEIRV(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SEIRV <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  Y = with(pars$Xpar[[i]], c*I)
  return(Y)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SEIRV model.
#' @inheritParams ramp.xde::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SEIRV <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i), {
    H <- S +E+ I+R+V
    return(H)
  })}

## -----------------------------------------------------------------------------
#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SEIRV model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SEIRV <- function(y, pars, i) {
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

## -----------------------------------------------------------------------------
#' @title Parse the output of deSolve and return variables for the SEIRV model
#' @description Implements [parse_deout_X] for the SEIRV model
#' @inheritParams ramp.xde::parse_deout_X
#' @return none
#' @export
parse_deout_X.SEIRV <- function(deout, pars, i) {
  time = deout[,1]
  with(pars$ix$X[[i]],{
    S = deout[,S_ix+1]
    E = deout[,E_ix+1]
    I = deout[,I_ix+1]
    R = deout[,R_ix+1]
    V = deout[,V_ix+1]
    H = S+E+I+R+V
    return(list(time=time, S=S, E=E,I=I, R=R, V=V,H=H))
  })}

## -----------------------------------------------------------------------------
#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIRV model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIRV <- function(varslist, pars,i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}

## -----------------------------------------------------------------------------
#' @title Compute the HTC for the SEIRV model
#' @description Implements [HTC] for the SEIRV model with demography.
#' @inheritParams ramp.xde::HTC
#' @return a [numeric] vector
#' @export
HTC.SEIRV <- function(pars, i) {
  with(pars$Xpar[[i]],
       HTC <- c/r,
       return(HTC)
  )
}


## -----------------------------------------------------------------------------

#' Add lines for the density of infected individuals for the SEIRV model
#'
#' @param XH a list with the outputs of parse_deout_X_SEIRV
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


## -----------------------------------------------------------------------------
#' Plot the density of infected individuals for the SEIRV model
#'
#' @inheritParams ramp.xde::xde_plot_X
#' @export
xde_plot_X.SEIRV = function(pars, i=1, clrs=c("black","darkblue","darkred","darkgreen","purple"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})
  
  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "No of. Infected", xlab = "Time"))
  
  
  xde_lines_X_SEIRV(vars$XH[[i]], pars$Hpar[[i]]$nStrata, clrs, llty)
}

## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [dXdt] for the SEIRV model
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SEIRV<- function(t, y, pars, i) {
  # do not change this
  foi <- pars$FoI[[i]]
  
  # attach the variables by name
  with(list_Xvars(y, pars, i),{
    # compute H (if it isn't one of the variables)
    H <- F_H(t, y, pars, i)
    
    # expose the parameters (see make_Xpar_SEIRV)
    with(pars$Xpar[[i]], {
      # compute the derivatives
      dS <- (1-alpha)*Births(t, H, pars,i) - foi*S + dHdt(t, S, pars,i)+ eps*R
      dE <- foi*S - tau*E + dHdt(t, E, pars,i)
      dI <- tau*E - r*I + dHdt(t, I, pars,i)
      dR <- (1-delta)*r*I -sig*R -eps*R+ dHdt(t, R, pars,i)
      dV <- alpha*Births(t, H, pars,i) + delta*r*I + sig*R + dHdt(t, V, pars,i)
      
       # concatenate the derivatives
      derivs = c(dS, dE, dI, dR,dV)
      
      # return the derivatives
      return(derivs)
    })
  })
}

## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
#' @title Setup Xinits.SEIRV
#' @description Implements [setup_Xinits] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SEIRV = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_SEIRV(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SEIRV model.
#' @inheritParams ramp.xde::make_indices_X
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

## -----------------------------------------------------------------------------
#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xde::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SEIRV <- function(y, pars, i) {
  with(pars$ix$X[[i]],
       return(list(
         S = y[S_ix],
         E = y[E_ix],
         I = y[I_ix],
         R = y[R_ix],
         V = y[V_ix]
       )
       ))
}

## -----------------------------------------------------------------------------
#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.SEIRV <- function(pars, i){
  with(pars$Xinits[[i]], return(c(S,E,I,R,V)))
}


## -----------------------------------------------------------------------------
#' @title Update inits for the SEIRV human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.SEIRV <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_SEIRV(pars, list(), S0=S,E0= E, I0=I, R0=R,V0 =V)
    return(pars)
  })}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Setup Xpar.SEIRV
#' @description Implements [setup_Xpar] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SEIRV = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SEIRV(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SEIRV <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  Y = with(pars$Xpar[[i]], c*I)
  return(Y)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SEIRV model.
#' @inheritParams ramp.xde::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SEIRV <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i), {
    H <- S +E+ I+R+V
    return(H)
  })}

## -----------------------------------------------------------------------------
#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SEIRV model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SEIRV <- function(y, pars, i) {
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

## -----------------------------------------------------------------------------
#' @title Parse the output of deSolve and return variables for the SEIRV model
#' @description Implements [parse_deout_X] for the SEIRV model
#' @inheritParams ramp.xde::parse_deout_X
#' @return none
#' @export
parse_deout_X.SEIRV <- function(deout, pars, i) {
  time = deout[,1]
  with(pars$ix$X[[i]],{
    S = deout[,S_ix+1]
    E = deout[,E_ix+1]
    I = deout[,I_ix+1]
    R = deout[,R_ix+1]
    V = deout[,V_ix+1]
    H = S+E+I+R+V
    return(list(time=time, S=S, E=E,I=I, R=R, V=V,H=H))
  })}

## -----------------------------------------------------------------------------
#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIRV model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIRV <- function(varslist, pars,i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}

## -----------------------------------------------------------------------------
#' @title Compute the HTC for the SEIRV model
#' @description Implements [HTC] for the SEIRV model with demography.
#' @inheritParams ramp.xde::HTC
#' @return a [numeric] vector
#' @export
HTC.SEIRV <- function(pars, i) {
  with(pars$Xpar[[i]],
       HTC <- c/r,
       return(HTC)
  )
}


## -----------------------------------------------------------------------------

#' Add lines for the density of infected individuals for the SEIRV model
#'
#' @param XH a list with the outputs of parse_deout_X_SEIRV
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


## -----------------------------------------------------------------------------
#' Plot the density of infected individuals for the SEIRV model
#'
#' @inheritParams ramp.xde::xde_plot_X
#' @export
xde_plot_X.SEIRV = function(pars, i=1, clrs=c("black","darkblue","darkred","darkgreen","purple"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})
  
  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "No of. Infected", xlab = "Time"))
  
  
  xde_lines_X_SEIRV(vars$XH[[i]], pars$Hpar[[i]]$nStrata, clrs, llty)
}

## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [dXdt] for the SEIRV model
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SEIRV<- function(t, y, pars, i) {
  # do not change this
  foi <- pars$FoI[[i]]
  
  # attach the variables by name
  with(list_Xvars(y, pars, i),{
    # compute H (if it isn't one of the variables)
    H <- F_H(t, y, pars, i)
    
    # expose the parameters (see make_Xpar_SEIRV)
    with(pars$Xpar[[i]], {
      # compute the derivatives
      dS <- (1-alpha)*Births(t, H, pars,i) - foi*S + dHdt(t, S, pars,i)+ eps*R
      dE <- foi*S - tau*E + dHdt(t, E, pars,i)
      dI <- tau*E - r*I + dHdt(t, I, pars,i)
      dR <- (1-delta)*r*I -sig*R -eps*R+ dHdt(t, R, pars,i)
      dV <- alpha*Births(t, H, pars,i) + delta*r*I + sig*R + dHdt(t, V, pars,i)
      
       # concatenate the derivatives
      derivs = c(dS, dE, dI, dR,dV)
      
      # return the derivatives
      return(derivs)
    })
  })
}

## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
#' @title Setup Xinits.SEIRV
#' @description Implements [setup_Xinits] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SEIRV = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_SEIRV(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SEIRV model.
#' @inheritParams ramp.xde::make_indices_X
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

## -----------------------------------------------------------------------------
#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xde::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SEIRV <- function(y, pars, i) {
  with(pars$ix$X[[i]],
       return(list(
         S = y[S_ix],
         E = y[E_ix],
         I = y[I_ix],
         R = y[R_ix],
         V = y[V_ix]
       )
       ))
}

## -----------------------------------------------------------------------------
#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.SEIRV <- function(pars, i){
  with(pars$Xinits[[i]], return(c(S,E,I,R,V)))
}


## -----------------------------------------------------------------------------
#' @title Update inits for the SEIRV human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.SEIRV <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_SEIRV(pars, list(), S0=S,E0= E, I0=I, R0=R,V0 =V)
    return(pars)
  })}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' @title Setup Xpar.SEIRV
#' @description Implements [setup_Xpar] for the SEIRV model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SEIRV = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SEIRV(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SEIRV <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  Y = with(pars$Xpar[[i]], c*I)
  return(Y)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SEIRV model.
#' @inheritParams ramp.xde::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SEIRV <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i), {
    H <- S +E+ I+R+V
    return(H)
  })}

## -----------------------------------------------------------------------------
#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SEIRV model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SEIRV <- function(y, pars, i) {
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

## -----------------------------------------------------------------------------
#' @title Parse the output of deSolve and return variables for the SEIRV model
#' @description Implements [parse_deout_X] for the SEIRV model
#' @inheritParams ramp.xde::parse_deout_X
#' @return none
#' @export
parse_deout_X.SEIRV <- function(deout, pars, i) {
  time = deout[,1]
  with(pars$ix$X[[i]],{
    S = deout[,S_ix+1]
    E = deout[,E_ix+1]
    I = deout[,I_ix+1]
    R = deout[,R_ix+1]
    V = deout[,V_ix+1]
    H = S+E+I+R+V
    return(list(time=time, S=S, E=E,I=I, R=R, V=V,H=H))
  })}

## -----------------------------------------------------------------------------
#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIRV model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIRV <- function(varslist, pars,i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}

## -----------------------------------------------------------------------------
#' @title Compute the HTC for the SEIRV model
#' @description Implements [HTC] for the SEIRV model with demography.
#' @inheritParams ramp.xde::HTC
#' @return a [numeric] vector
#' @export
HTC.SEIRV <- function(pars, i) {
  with(pars$Xpar[[i]],
       HTC <- c/r,
       return(HTC)
  )
}


## -----------------------------------------------------------------------------

#' Add lines for the density of infected individuals for the SEIRV model
#'
#' @param XH a list with the outputs of parse_deout_X_SEIRV
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


## -----------------------------------------------------------------------------
#' Plot the density of infected individuals for the SEIRV model
#'
#' @inheritParams ramp.xde::xde_plot_X
#' @export
xde_plot_X.SEIRV = function(pars, i=1, clrs=c("black","darkblue","darkred","darkgreen","purple"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})
  
  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "No of. Infected", xlab = "Time"))
  
  
  xde_lines_X_SEIRV(vars$XH[[i]], pars$Hpar[[i]]$nStrata, clrs, llty)
}

