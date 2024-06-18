# specialized methods for the human SIPw dts model

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIPw_dts model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SIPw_dts <- function(t, y, pars, i) {
  with(list_Xvars(y, pars, i),
    with(pars$Xpar[[i]],{
      X = c*I
      return(X)
    })
  )
}

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIPw_dts model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SIPw_dts <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i),return(S+I+P) )
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SIPw_dts model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SIPw_dts <- function(varslist, pars, i) {
  pr = with(varslist$XH[[i]], I/H)
  return(pr)
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SIPw_dts model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SIPw_dts <- function(y, pars,i) {
  with(pars$Xpar[[i]], b)
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIPw_dts model.
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SIPw_dts <- function(t, y, pars, i){

  attack <- pars$AR[[i]]

  with(list_Xvars(y, pars, i),{
    H <- F_H(t, y, pars, i)
    with(pars$Xpar[[i]], {

      St <- (1-attack)*(S+r*I) + eta*P - xi*S
      It <- (1-r)*I + attack*(1-rho)*(S+r*I) - xi*I
      Pt <- xi*(S+I) + attack*rho*(S+r*I) + (1-eta)*P
      wt <- w + ar

      St <- dHdt(t, St, i) + Births(t, H, pars, i)
      It <- dHdt(t, It, i)
      Pt <- dHdt(t, Pt, i)
      wt <- dHdt(t, wt, i)

      return(c(St, It, Pt, wt))
    })
  })
}

#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xde::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SIPw_dts <- function(y, pars, i) {
  with(pars$ix$X[[i]],
      return(list(
        S = y[S_ix],
        I = y[I_ix],
        P = y[P_ix],
        w = y[w_ix]
      )
  ))
}

#' @title Compute the HTC for the SIPw_dts model
#' @description Implements [HTC] for the SIPw_dts model with demography.
#' @inheritParams ramp.xde::HTC
#' @return a [numeric] vector
#' @export
HTC.SIPw_dts <- function(pars, i) {
  with(pars$Xpar[[i]],
       return((1-rho)*b/(r+xi)*xi/(eta+xi))
  )
}

#' @title Setup Xpar.SIPw_dts
#' @description Implements [setup_Xpar] for the SIPw_dts model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SIPw_dts = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SIPw_dts(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}

#' @title Setup Xinits.SIPw_dts
#' @description Implements [setup_Xinits] for the SIPw_dts model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SIPw_dts = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = make_Xinits_SIPw_dts(pars$Hpar[[i]]$nStrata, Xopts, H0=pars$Hpar[[i]]$H)
  return(pars)
}

#' @title Make parameters for SIPw_dts human model, with defaults
#' @param nStrata the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @param rho probability of successful treatment upon infection
#' @param eta prophylaxis waning rate
#' @param xi background treatment rate
#' @return a [list]
#' @export
make_Xpar_SIPw_dts = function(nStrata, Xopts=list(),
                         b=0.55, r=1/180, c=0.15,
                         rho=.1, eta=1/25, xi=1/365){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SIPw_dts")

    Xpar$b = checkIt(b, nStrata)
    Xpar$c = checkIt(c, nStrata)
    Xpar$r = checkIt(r, nStrata)
    Xpar$rho = checkIt(rho, nStrata)
    Xpar$eta = checkIt(eta, nStrata)
    Xpar$xi = checkIt(xi, nStrata)

    return(Xpar)
})}


#' @title Make initial values for the SIPw_dts human model, with defaults
#' @param nStrata the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param H0 the initial human population density
#' @param S0 the initial values of the variable S
#' @param I0 the initial values of the variable I
#' @param P0 the initial values of the variable P
#' @param w0 the initial values of the tracking variable w
#' @return a [list]
#' @export
make_Xinits_SIPw_dts = function(nStrata, Xopts = list(), H0=NULL, S0=NULL,
                           I0=1, P0=0, w0=0){with(Xopts,{
  if(is.null(S0)) S0=H0-I0-P0
  S = checkIt(S0, nStrata)
  I = checkIt(I0, nStrata)
  P = checkIt(P0, nStrata)
  w = checkIt(w0, nStrata)
  return(list(S=S, I=I, P=P, w=w))
})}

#' @title Parse the output of deSolve and return variables for the SIPw_dts model
#' @description Implements [parse_outputs_X] for the SIPw_dts model
#' @inheritParams ramp.xde::parse_outputs_X
#' @return none
#' @export
parse_outputs_X.SIPw_dts <- function(outputs, pars, i) {
  time = outputs[,1]
  with(pars$ix$X[[i]],{
    S = outputs[,S_ix+1]
    I = outputs[,I_ix+1]
    P = outputs[,P_ix+1]
    w = outputs[,w_ix+1]
    H = S+I+P
    return(list(time=time,S=S,I=I,P=P,H=H,w=w))
})}

#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SIPw_dts model.
#' @inheritParams ramp.xde::make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.SIPw_dts <- function(pars, i) {with(pars,{

  S_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(S_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(I_ix, 1)

  P_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(P_ix, 1)

  w_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(w_ix, 1)

  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(S_ix=S_ix, I_ix=I_ix, P_ix=P_ix, w_ix=w_ix)
  return(pars)
})}

#' @title Update inits for the SIPw_dts human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.SIPw_dts <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
  pars$Xinits[[i]] = make_Xinits_SIPw_dts(pars$Hpar[[i]]$nStrata, list(), S0=S, I0=I, P0=P, w0=w)
  return(pars)
})}


#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return none
#' @export
get_inits_X.SIPw_dts <- function(pars, i){with(pars$Xinits[[i]],{
  c(S, I, P, w)
})}


#' Plot the density of infected individuals for the SIPw_dts model
#'
#' @inheritParams ramp.xde::xds_plot_X
#' @export
xds_plot_X.SIPw_dts = function(pars, i=1, clrs=c("darkblue", "darkred", "darkgreen"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "Time"))

  add_lines_X_SIPw_dts(vars$XH[[i]], pars, clrs, llty)
}


#' Add lines for the density of infected individuals for the SIPw_dts model
#'
#' @param XH a list with the outputs of parse_outputs_X.SIPw_dts
#' @param pars a list that defines an `ramp.dts` model (*e.g.*,  generated by `dts_setup()`)
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
add_lines_X_SIPw_dts = function(XH, pars, clrs=c("darkblue", "darkred", "darkgreen"), llty=1){
  with(XH,{
    if(pars$Hpar[[1]]$nStrata==1) {
      lines(time, S, col=clrs[1], lty = llty[1])
      lines(time, I, col=clrs[2], lty = llty[1])
      lines(time, P, col=clrs[3], lty = llty[1])
    }
    if(pars$Hpar[[1]]$nStrata>1){
      if (length(clrs)==1) clrs=matrix(clrs, 3, pars$Hpar[[i]]$nStrata)
      if (length(llty)==1) llty=rep(llty, pars$Hpar[[i]]$nStrata)
      for(i in 1:pars$Hpar[[1]]$nStrata){
        lines(time, S[,i], col=clrs[1,i], lty = llty[i])
        lines(time, I[,i], col=clrs[2,i], lty = llty[i])
        lines(time, P[,i], col=clrs[3,i], lty = llty[i])
      }
    }
  })}

