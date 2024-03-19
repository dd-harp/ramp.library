# specialized methods for the human SIR model

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIR model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SIR <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  X = with(pars$Xpar[[i]], c*I)
  return(X)
}


#' @title Size of human population
#' @description Implements [F_H] for the SIR model.
#' @inheritParams ramp.xde::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SIR <- function(t, y, pars, i){
  S = y[pars$ix$X[[i]]$S_ix]
  I = y[pars$ix$X[[i]]$I_ix]
  R = y[pars$ix$X[[i]]$R_ix]
  return(S+I+R)
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SIR model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SIR <- function(varslist, pars, i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}


#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SIR model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SIR <- function(y, pars, i) {
  with(pars$Xpar[[i]], b)
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIR model with demography.
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SIR <- function(t, y, pars, i) {

  foi <- pars$FoI[[i]]

  with(pars$ix$X[[i]],{
    S = y[S_ix]
    I = y[I_ix]
    R = y[R_ix]
    H <- F_H(t, y, pars, i)

    with(pars$Xpar[[i]],{

      dS <- Births(t, H, pars) - foi*S + dHdt(t, S, pars)
      dI <- foi*S - r*I + dHdt(t, I, pars)
      dR <- r*I + dHdt(t, R, pars)

      return(c(dS, dI, dR))
    })
  })
}

#' @title Setup Xpar.SIR
#' @description Implements [setup_Xpar] for the SIR model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SIR = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SIR(pars$nStrata, Xopts)
  return(pars)
}

#' @title Setup Xinits.SIR
#' @description Implements [setup_Xinits] for the SIR model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SIR = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars,make_Xinits_SIR(nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}

#' @title Make parameters for SIR human model, with defaults
#'
#' @param nStrata is the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @return a [list]
#'
#' @export
make_Xpar_SIR = function(nStrata, Xopts=list(),
                         b=0.55, r=1/180, c=0.15){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SIR")

    Xpar$b = checkIt(b, nStrata)
    Xpar$c = checkIt(c, nStrata)
    Xpar$r = checkIt(r, nStrata)

    return(Xpar)
  })}

#' @title Make initial values for the SIR human model, with defaults
#' @param nStrata the number of strata in the model
#' @param Xopts a [list] to overwrite defaults
#' @param H0 the initial human population density
#' @param S0 the initial values of the parameter S
#' @param I0 the initial values of the parameter I
#' @param R0 the initial values of the parameter R
#' @return a [list]
#' @export
make_Xinits_SIR = function(nStrata, Xopts = list(), H0=NULL, S0=NULL, I0=1, R0=0){with(Xopts,{
  if(is.null(S0)) S0 = H0 - I0 - R0
  stopifnot(is.numeric(S0))
  S = checkIt(S0, nStrata)
  I = checkIt(I0, nStrata)
  R = checkIt(R0, nStrata)
  return(list(S=S, I=I, R=R))
})}

#' @title Parse the output of deSolve and return variables for the SIR model
#' @description Implements [parse_deout_X] for the SIR model
#' @inheritParams ramp.xde::parse_deout_X
#' @return none
#' @export
parse_deout_X.SIR <- function(deout, pars,i) {
  time = deout[,1]
  with(pars$ix$X[[i]],{
    S = deout[,S_ix+1]
    I = deout[,I_ix+1]
    R = deout[,R_ix+1]
    return(list(time=time, S=S, I=I, R=R, H=H))
})}

#' @title Compute the HTC for the SIR model
#' @description Implements [HTC] for the SIR model with demography.
#' @inheritParams ramp.xde::HTC
#' @return a [numeric] vector
#' @export
HTC.SIR <- function(pars) {
  with(pars$Xpar,
       return(c/r)
  )
}

#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SIR model.
#' @inheritParams ramp.xde::make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.SIR <- function(pars, i) {with(pars,{

  S_ix <- seq(from = max_ix+1, length.out=nStrata)
  max_ix <- tail(S_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=nStrata)
  max_ix <- tail(I_ix, 1)

  R_ix <- seq(from = max_ix+1, length.out=nStrata)
  max_ix <- tail(R_ix, 1)

  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(S_ix=S_ix, I_ix=I_ix, R_ix=R_ix)

  return(pars)
})}

#' @title Update inits for the SIR human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.SIR <- function(pars, y0, i) {
  with(pars$ix$X[[i]],{
  S = y0[S_ix]
  I = y0[I_ix]
  R = y0[R_ix]
  pars = make_Xinits_SIR(pars, list(), S, I, R)
  return(pars)
})}


#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xde::get_inits_X
#' @return none
#' @export
get_inits_X.SIR <- function(pars, i){
  with(pars$Xinits[[i]], c(S, I, R))
}

#' Plot the density of infected individuals for the SIR model
#'
#' @inheritParams ramp.xde::xde_plot_X
#' @export
xde_plot_X.SIR = function(pars, i, clrs=c("black", "darkred", "darkblue"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "Time"))

  xde_lines_X_SIR(vars$XH[[i]], pars, clrs, llty)
}


#' Add lines for the density of infected individuals for the SIR model
#'
#' @param XH a list with the outputs of parse_deout_X_SIS
#' @param pars a list that defines an `ramp.xde` model (*e.g.*,  generated by [ramp.xde::xde_setup])
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
#'
xde_lines_X_SIR = function(XH, pars, clrs=c("black", "darkred", "darkblue"), llty=1){
  with(XH,{
    if(pars$nStrata==1) {
      lines(time, S, col=clrs[1], lty = llty[1])
      lines(time, I, col=clrs[2], lty = llty[1])
      lines(time, R, col=clrs[3], lty = llty[1])
    }
    if(pars$nStrata>1){
      if (length(clrs)==1) clrs=matrix(clrs, 3, pars$nStrata)
      if (length(llty)==1) llty=rep(llty, pars$nStrata)
      for(i in 1:pars$nStrata){
        lines(time, S[,i], col=clrs[1,i], lty = llty[i])
        lines(time, I[,i], col=clrs[2,i], lty = llty[i])
        lines(time, R[,i], col=clrs[3,i], lty = llty[i])
      }
    }
  })}
