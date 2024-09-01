# specialized methods for the human SIPw model, formulated for differential and difference equations


#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIPw-xde model.
#' @inheritParams ramp.xds::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SIPw <- function(t, y, pars, i){

  foi <- pars$FoI[[i]]
  Hpar <- pars$Hpar[[i]]

  with(list_Xvars(y, pars, i),{
    with(pars$Xpar[[i]], {

      dS <- -foi*S         - xi*S  + r*I  + eta*P + dHdt(t, S, Hpar) + Births(t, H, Hpar)
      dI <-  foi*(1-rho)*S - (xi+sigma)*I  - r*I          + dHdt(t, I, Hpar)
      dP <-  foi*rho*S     + xi*(S+I) + sigma*I    - eta*P + dHdt(t, P, Hpar)
      dw <- foi + dHdt(t, w, Hpar)

      return(c(dS, dI, dP, dw))
    })
  })
}

#' @title Derivatives for human population
#' @description Implements [Update_Xt] for the SIPw dts model
#' @inheritParams ramp.xds::Update_Xt
#' @return a [numeric] vector
#' @export
Update_Xt.SIPw <- function(t, y, pars, i){

  ar <- pars$AR[[i]]
  Hpar <- pars$Hpar[[i]]

  with(list_Xvars(y, pars, i),{
    with(pars$Xpar[[i]], {

      St <- (1-ar)*(1-xi)*S       + (1-ar)*r*(1-xisig)*I             + eta*P
      It <- ar*(1-rho)*(1-xi)*S    + (1-r)*(1-xisig)*I + ar*(1-rho)*r*(1-xisig)*I
      Pt <- ar*rho*(1-xi)*S+ xi*S  + xisig*I + ar*rho*r*(1-xisig)*I  + (1-eta)*P
      wt <- w + ar

      St <- St+dHdt(t, St, Hpar) + Births(t, H, Hpar)
      It <- It+dHdt(t, It, Hpar)
      Pt <- Pt+dHdt(t, Pt, Hpar)
      wt <- wt+dHdt(t, wt, Hpar)

      return(c(St, It, Pt, wt))
    })
  })
}

#' @title Setup the Xpar for the SIPw_xde model
#' @description implements [make_Xpar] for the SIPw model
#' @inheritParams ramp.xds::make_Xpar
#' @return a [list] vector
#' @export
make_Xpar.SIPw = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = create_Xpar_SIPw(pars$nStrata[i], Xopts)
  return(pars)
}

#' @title Make parameters for SIPw_xde human model, with defaults
#' @param nStrata the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @param rho probability of successful treatment upon infection
#' @param sigma probability of treatment for treated, above background
#' @param xi background treatment rate
#' @param eta rate of loss of chemo-protection
#' @return a [list]
#' @export
create_Xpar_SIPw= function(nStrata, Xopts=list(), b=0.55, c=0.15, r=1/180,
                             rho=.1, sigma=1/730, xi=1/365,  eta=1/25){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SIPw")

    Xpar$b = checkIt(b, nStrata)
    Xpar$c = checkIt(c, nStrata)
    Xpar$r = checkIt(r, nStrata)
    Xpar$rho = checkIt(rho, nStrata)
    Xpar$eta = checkIt(eta, nStrata)
    Xpar$sigma = checkIt(sigma, nStrata)
    Xpar$xi = checkIt(xi, nStrata)

    return(Xpar)
})}


#' @title Size of effective infectious human population
#' @description Implements [F_X] for SIPw models
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SIPw <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  X = with(pars$Xpar[[i]], c*I)
  return(X)
}

#' @title Size of effective infectious human population
#' @description Implements [F_X] for SIPw models
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SIPw <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i),return(H))
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for SIPw models
#' @inheritParams ramp.xds::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SIPw <- function(vars, Xpar) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for SIPw models
#' @inheritParams ramp.xds::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SIPw <- function(y, pars,i) {
  with(pars$Xpar[[i]], b)
}


#' @title Return the SIPw model variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xds::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SIPw <- function(y, pars, i) {
  with(pars$ix$X[[i]],{
       S = y[S_ix]
       I = y[I_ix]
       P = y[P_ix]
       w = y[w_ix]
       H = S+I+P
       return(list(S=S,I=I,P=P,w=w,H=H))})
}

#' @title Return the SIPw model variables as a list, returned from Update_Xt.SIPw
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xds::put_Xvars
#' @return a [list]
#' @export
put_Xvars.SIPw <- function(Xvars, y, pars, i) {
  with(pars$ix$X[[i]],
    with(as.list(Xvars),{
      y[S_ix] = S
      y[I_ix] = I
      y[P_ix] = P
      y[w_ix] = w
      return(y)
}))}

#' @title Compute the HTC for the SIPw_xde model
#' @description Implements [HTC] for the SIPw_xde model with demography.
#' @inheritParams ramp.xds::HTC
#' @return a [numeric] vector
#' @export
HTC.SIPw <- function(pars, i) {
  with(pars$Xpar[[i]],
       return((1-rho)*b/(r+xi)*xi/(eta+xi))
  )
}





#' @title Setup Xinits.SIPw
#' @description Implements [make_Xinits] for the SIPw models
#' @inheritParams ramp.xds::make_Xinits
#' @return a [list] vector
#' @export
make_Xinits.SIPw = function(pars,H, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars,create_Xinits_SIPw(pars$nStrata[i],H, Xopts))
  return(pars)
}

#' @title Make initial values for a SIPw human model, with defaults
#' @param nStrata the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param H the initial human population density
#' @param I the initial values of the variable I
#' @param P the initial values of the variable P
#' @param w the initial values of the tracking variable w
#' @return a [list]
#' @export
create_Xinits_SIPw = function(nStrata, H, Xopts = list(),
                           I=1, P=0, w=0){with(Xopts, {
    S = unname(as.vector(checkIt(H -I-P, nStrata)))
    I = unname(as.vector(checkIt(I, nStrata)))
    P = unname(as.vector(checkIt(P, nStrata)))
    w = unname(as.vector(checkIt(w, nStrata)))
    return(list(S=S, I=I, P=P, w=w))
})}

#' @title Parse the output of deSolve and return variables for SIPw models
#' @description Implements [parse_Xorbits] for SIPw models
#' @inheritParams ramp.xds::parse_Xorbits
#' @return none
#' @export
parse_Xorbits.SIPw <- function(outputs, pars, i) {with(pars$ix$X[[i]],{
    S <-outputs[,S_ix]
    I <- outputs[,I_ix]
    P <-outputs[,P_ix]
    w <- outputs[,w_ix]
    H <- S+I+P
    ni <- pars$Xpar[[i]]$c*I/H
    true_pr <- I/H
    return(list(time=time,S=S,I=I,P=P,H=H,w=w,ni=ni, true_pr= true_pr))
  })}

#' @title Add indices for human population to parameter list
#' @description Implements [make_X_indices] for SIPw models
#' @inheritParams ramp.xds::make_X_indices
#' @return none
#' @importFrom utils tail
#' @export
make_X_indices.SIPw <- function(pars, i) {with(pars,{

  S_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(S_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(I_ix, 1)

  P_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(P_ix, 1)

  w_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(w_ix, 1)

  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(S_ix=S_ix, I_ix=I_ix, P_ix=P_ix, w_ix=w_ix)
  return(pars)
})}

#' @title Update inits for SIPw models from a vector of states
#' @inheritParams ramp.xds::update_Xinits
#' @return none
#' @export
update_Xinits.SIPw <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars$Xinits[[i]] = create_Xinits_SIPw(pars$nStrata[i], pars$H0,list(), I=I, P=P, w=w)
    return(pars)
  })}


#' @title Return initial values as a vector from a SIPw model
#' @description This method dispatches on the type of `pars$Xpar[[i]]`
#' @inheritParams ramp.xds::get_Xinits
#' @return none
#' @export
get_Xinits.SIPw <- function(pars, i){pars$Xinits[[i]]}


#' Plot the density of infected individuals for the SIPw model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @export
xds_plot_X.SIPw = function(pars, i=1, clrs=c("darkblue", "darkred", "darkgreen"), llty=1, add=FALSE){
  XH = pars$outputs$orbits$XH[[i]]
  times = pars$outputs$time

  if(add==FALSE)
    plot(times, 0*times, type = "n", ylim = c(0, max(XH$H)),
         ylab = "# Infected", xlab = "Time")

  xds_lines_X_SIPw(times, XH, pars$nStrata[i], clrs, llty)
}


#' Add lines for the density of infected individuals for the SIP model
#'
#' @param times time points for the observations
#' @param XH a list with the outputs of parse_Xorbits.SIP
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
xds_lines_X_SIPw = function(times, XH, nStrata, clrs=c("darkblue", "darkred", "darkgreen"), llty=1){
  if (length(llty)< nStrata) llty = rep(llty, nStrata)
  with(XH,{
    if(nStrata == 1){
      lines(times, S, col=clrs[1], lty = llty)
      lines(times, I, col=clrs[2], lty = llty)
      lines(times, P, col=clrs[3], lty = llty)
    } else {
      for(i in 1:nStrata)
        lines(times, S[,i], col=clrs[1], lty = llty[i])
        lines(times, I[,i], col=clrs[2], lty = llty[i])
        lines(times, P[,i], col=clrs[3], lty = llty[i])
    }})
}


#' @title Compute the steady states for the SIP model as a function of the daily foi
#' @description Compute the steady state of the SIP model as a function of the daily eir.
#' @inheritParams ramp.xds::xde_steady_state_X
#' @return the steady states as a named vector
#' @export
xde_steady_state_X.SIPw = function(foi, H, Xpar){with(Xpar,{
  Ieq = (foi*H*eta*(1-rho))/((foi+r+xi)*(eta+xi) +foi*(r-eta)*rho)
  Peq  = (H*xi*(foi+r+xi) + (foi*H*r*rho))/((foi+r+xi)*(eta+xi) +foi*(r-eta)*rho)
  Seq = H -Ieq - Peq
  return(list(S=as.vector(Seq), I=as.vector(Ieq), P = as.vector(Peq)))
})}

