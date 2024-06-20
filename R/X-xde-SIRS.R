#' @title Compute the derivatives for parasite infection dynamics in human population strata
#' @description Implements [dXdt] for the SIRS model
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SIRS <- function(t, y, pars, i) {
  # do not change this
  foi <- pars$FoI[[i]]
  Hpar <- pars$Hpar[[i]]

  # attach the variables by name
  with(list_Xvars(y, pars, i),{
    # compute H (if it isn't one of the variables)
    H <- F_H(t, y, pars, i)

    # expose the parameters (see make_Xpar_SIRS)
    with(pars$Xpar[[i]], {
      # compute the derivatives
      dS <- Births(t, H, Hpar) + dHdt(t, S, Hpar) - foi*S + gam*R
      dI <- foi*S- r*I + dHdt(t, I, Hpar)
      dR <- r*I - gam*R + dHdt(t, R, Hpar)

      # concatenate the derivatives
      derivs = c(dS, dI, dR)

      # return the derivatives
      return(derivs)
    })
  })
}



#' @title Make initial values for the SIRS human model, with defaults
#' @param nStrata the number of strata in the model
#' @param Xopts a [list] to overwrite defaults
#' @param H0 the initial value for H
#' @param S0 the initial value for S
#' @param I0 the initial value for I
#' @param R0 the initial values for R
#' @return a [list]
#' @export
make_Xinits_SIRS = function(nStrata, Xopts = list(), H0= NULL, S0=NULL, I0=1, R0 = 1){with(Xopts,{
  if(is.null(S0)) S0 = H0-(I0+R0)
  stopifnot(is.numeric(S0))
  S = checkIt(S0, nStrata)
  I = checkIt(I0, nStrata)
  R = checkIt(R0, nStrata)
  return(list(S=S, I=I, R =R))
})}





#' @title Setup Xinits.SIRS
#' @description Implements [setup_Xinits] for the SIRS model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SIRS = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_SIRS(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}





#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SIRS model.
#' @inheritParams ramp.xde::make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.SIRS <- function(pars, i) {with(pars,{

  S_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(S_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(I_ix, 1)

  R_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(R_ix, 1)



  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(S_ix=S_ix, I_ix=I_ix, R_ix=R_ix)
  return(pars)
})}





#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xde::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SIRS <- function(y, pars, i) {
  with(pars$ix$X[[i]],
       return(list(
         S = y[S_ix],
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
get_inits_X.SIRS <- function(pars, i){
  with(pars$Xinits[[i]], return(c(S,I,R)))
}





#' @title Update inits for the SIRS human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.SIRS <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_SIRS(pars, list(), S0=S, I0=I, R0=R)
    return(pars)
  })}



#' @title Make parameters for SIRS human model, with defaults
#' @param nStrata is the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param b the proportion of infective bites that cause an infection
#' @param r the the duration of an infection
#' @param c the proportion of bites on infected humans that infect a mosquito
#' @param gam the rate of loss of immunity
#' @return a [list]
#' @export
make_Xpar_SIRS = function(nStrata, Xopts=list(),
                          b=0.55, r=1/180, c=0.15,gam=0.5){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SIRS")

    Xpar$b = checkIt(b, nStrata)
    Xpar$c = checkIt(c, nStrata)
    Xpar$r = checkIt(r, nStrata)
    Xpar$gam = checkIt(gam, nStrata)

    return(Xpar)
  })}





#' @title Setup Xpar.SIRS
#' @description Implements [xde_setup_Xpar] for the SIRS model
#' @inheritParams ramp.xde::xde_setup_Xpar
#' @return a [list] vector
#' @export
xde_setup_Xpar.SIRS = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SIRS(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}





#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SIRS <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  Y = with(pars$Xpar[[i]], c*I)
  return(Y)
}





#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SIRS model.
#' @inheritParams ramp.xde::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SIRS <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i), {
    H <- S + I+R
    return(H)
  })}





#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SIRS model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SIRS <- function(y, pars, i) {
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





#' @title Parse the output of deSolve and return variables for the SIRS model
#' @description Implements [parse_outputs_X] for the SIRS model
#' @inheritParams ramp.xde::parse_outputs_X
#' @return none
#' @export
parse_outputs_X.SIRS <- function(outputs, pars, i) {
  time = outputs[,1]
  with(pars$ix$X[[i]],{
    S = outputs[,S_ix+1]
    I = outputs[,I_ix+1]
    R = outputs[,R_ix+1]

    H = S+I+R
    return(list(time=time, S=S, I=I, R=R, H=H))
  })}





#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SIRS model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SIRS <- function(varslist, pars,i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}





#' @title Compute the HTC for the SIRS model
#' @description Implements [HTC] for the SIRS model with demography.
#' @inheritParams ramp.xde::HTC
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
#' @param XH a list with the outputs of parse_outputs_X_SIRS
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
xde_lines_X_SIRS = function(XH, nStrata, clrs=c("darkblue","darkred","darkgreen"), llty=1){
  with(XH,{
    if(nStrata==1) {
      lines(time, S, col=clrs[1], lty = llty[1])
      lines(time, I, col=clrs[2], lty = llty[1])
      lines(time, R, col=clrs[3], lty = llty[1])
    }
    if(nStrata>1){
      if (length(clrs)==3) clrs=matrix(clrs, 3, nStrata)
      if (length(llty)==1) llty=rep(llty, nStrata)

      for(i in 1:nStrata){
        lines(time, S[,i], col=clrs[1,i], lty = llty[i])
        lines(time, I[,i], col=clrs[2,i], lty = llty[i])
        lines(time, R[,i], col=clrs[3,i], lty = llty[i])
      }
    }
  })}





#' Plot the density of infected individuals for the SIRS model
#'
#' @inheritParams ramp.xde::xds_plot_X
#' @export
xds_plot_X.SIRS = function(pars, i=1, clrs=c("darkblue","darkred","darkgreen"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "No of. Infected", xlab = "Time"))


  xde_lines_X_SIRS(vars$XH[[i]], pars$Hpar[[i]]$nStrata, clrs, llty)
}
