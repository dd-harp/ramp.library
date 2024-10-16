# specialized methods for the human SIPd model

#' @title \eqn{\cal X} Component Derivatives for the `SIPd` Model
#' @description Compute the derivatives for SIPd compartmental model:
#' \deqn{
#' \begin{array}{rrrrcc}
#' dS/dt =& - (h +\xi) S &+ r I & + \eta P & + d{\cal H}(S) &+ B(t, H) \\
#' dI/dt =& (1-\rho) h S & - (r+\xi) I &&+ d{\cal H}(I) \\
#' dP/dt =& (\rho h +\xi)S & +\xi I & - \eta P &+ d{\cal H}(P)
#' \end{array}
#' }
#' where \eqn{H = S+I+P}; \eqn{B(t, H)} is the
#' time-dependent birth rate; and the \eqn{d{\cal H}}
#' operator computes derivatives for the demographic model \eqn{\cal H}.
#' @seealso The parameters are defined in [make_Xpar_SIPd]
#' @inheritParams ramp.xds::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SIPd <- function(t, y, pars, i){

  foi <- pars$FoI[[i]]

  with(list_Xvars(y, pars, i),{
    Hpar <- pars$Hpar[[i]]
    with(pars$Xpar[[i]], {
      if (t <= eta) {
        S_eta = 0
        I_eta = 0
        foi_eta = foi
      } else {
        S_eta   = lagvalue(t=t-eta, nr=pars$ix$X[[i]]$S_ix)
        I_eta   = lagvalue(t=t-eta, nr=pars$ix$X[[i]]$I_ix)
        foi_eta = lagderiv(t=t-eta, nr=pars$ix$X[[i]]$cfoi_ix)
      }

      dS <- Births(t, H, Hpar) - foi*S -xi*S + r*I  + (rho*foi_eta + xi)*S_eta + xi*I_eta + dHdt(t, S, Hpar)
      dI <- (1-rho)*foi*S - (r+xi)*I + dHdt(t, I, Hpar)
      dP <- rho*foi*S + xi*(S+I) - (rho*foi_eta + xi)*S_eta  - xi*I_eta + dHdt(t, P, Hpar)

      return(c(dS, dI, dP, foi))
    })
  })
}

#' @title Setup `Xpar` for an `SIPd`
#' @description Implements [setup_Xpar] for the SIPd model
#' @inheritParams ramp.xds::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SIPd = function(Xname, pars, i, Xopts=list()){
  pars = xds_dde(pars)
  pars$Xpar[[i]] = make_Xpar_SIPd(pars$nStrata[i], Xopts)
  return(pars)
}


#' @title Make parameters for SIPd human model, with defaults
#' @param nStrata the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @param rho probability of successful treatment upon infection
#' @param eta duration of chemo-protection
#' @param xi background treatment rate
#' @return a [list]
#' @export
make_Xpar_SIPd = function(nStrata, Xopts=list(),
                            b=0.55, r=1/180, c=0.15,
                            rho=.1, eta=25, xi=1/365){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SIPd")

    Xpar$b = checkIt(b, nStrata)
    Xpar$c = checkIt(c, nStrata)
    Xpar$r = checkIt(r, nStrata)
    Xpar$rho = checkIt(rho, nStrata)
    Xpar$eta = checkIt(eta, nStrata)
    Xpar$xi = checkIt(xi, nStrata)

    return(Xpar)
  })}


#' @title Return the parameters as a list
#' @description This method dispatches on the type of `pars$Xpar[[i]]`.
#' @inheritParams ramp.xds::set_Xpars
#' @return an **`xds`** object
#' @export
set_Xpars.SIPd <- function(pars, i=1, Xopts=list()) {
  nHabitats <- pars$nHabitats
  with(pars$Xpar[[i]], with(Xopts,{
    pars$Xpar[[i]]$b <- b
    pars$Xpar[[i]]$c <- c
    pars$Xpar[[i]]$r <- r
    pars$Xpar[[i]]$rho <- rho
    pars$Xpar[[i]]$eta <- eta
    pars$Xpar[[i]]$xi <- xi
    return(pars)
  }))}


#' @title Return the parameters as a list
#' @description This method dispatches on the type of `pars$Xpar[[i]]`.
#' @inheritParams ramp.xds::set_Xinits
#' @return an **`xds`** object
#' @export
set_Xinits.SIPd <- function(pars, i=1, Xopts=list()) {
  with(get_Xinits(pars, i), with(Xopts,{
    pars$Xinits[[i]]$S = get_H(pars,i)-I-P
    pars$Xinits[[i]]$I = I
    pars$Xinits[[i]]$P = P
    return(pars)
  }))}


#' @title Return the parameters as a list
#' @description Parameter values for the \eqn{i^{th}} host are
#' stored as `pars$Xpar[[i]]`. This returns the stored parameter
#' values as a list.
#' @inheritParams ramp.xds::get_Xpars
#' @return a [list]
#' @seealso [make_Xpar_SIPd]
#' @export
get_Xpars.SIPd <- function(pars, i=1) {
  with(pars$Xpar[[i]],list(b=b, c=c, r=r, rho=rho, eta=eta, xi=xi))
}

#' @title Derivatives for human population
#' @description Implements [Update_Xt] for the SIPd model.
#' @inheritParams ramp.xds::Update_Xt
#' @return a [numeric] vector
#' @export
Update_Xt.SIPd <- function(t, y, pars, i){

  attack <- pars$AR[[i]]

  with(list_Xvars(y, pars, i),{
    with(pars$Xpar[[i]], {

      St <- (1-attack)*(S+r*I) + eta*P - xi*S
      It <- (1-r)*I + attack*(1-rho)*(S+r*I) - xi*I
      Pt <- xi*(S+I) + attack*rho*(S+r*I) + (1-eta)*P

      St <- dHdt(t, St, pars$Hpar[[i]]) + Births(t, H, pars$Hpar[[i]])
      It <- dHdt(t, It, pars$Hpar[[i]])
      Pt <- dHdt(t, Pt, pars$Hpar[[i]])

      return(c(St, It, Pt))
    })
  })
}

#' @title Compute the steady states for the  dts SIPd model as a function of the daily EIR
#' @description Compute the steady state of the  dts SIPd model as a function of the daily eir.
#' @inheritParams ramp.xds::dts_steady_state_X
#' @return the steady states as a named vector
#' @export
dts_steady_state_X.SIPd = function(ar, H, Xpar){with(Xpar,{
  Iteq =(ar*H*eta*(1-rho))/((eta+xi)*(ar*(1-r)+(r+xi)+ar*(eta*(r-1)+r+xi)))
  Pteq = (H*(xi*(ar+r*(1-ar)+xi)+ar*(1+xi)*rho))/((eta+xi)*(ar*(1-r)+(r+xi)+ar*(eta*(r-1)+r+xi)))
  Steq = H-Iteq-Pteq
  return(c(S=Steq, I=Iteq,P =Pteq))
})}

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIPd model.
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SIPd <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  X = with(pars$Xpar[[i]], c*I)
  return(X)
}

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIPd model.
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SIPd <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i), return(H))
}


#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SIPd model.
#' @inheritParams ramp.xds::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SIPd <- function(y, pars,i) {
  with(pars$Xpar[[i]], b)
}



#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xds::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SIPd <- function(y, pars, i) {
  with(pars$ix$X[[i]],{
    S = y[S_ix]
    I = y[I_ix]
    P = y[P_ix]
    H = S+I+P
    return(list(S=S,I=I,P=P,H=H))})
}


#' @title Compute the HTC for the SIPd model
#' @description Implements [HTC] for the SIPd model with demography.
#' @inheritParams ramp.xds::HTC
#' @return a [numeric] vector
#' @export
HTC.SIPd <- function(pars, i) {
  with(pars$Xpar[[i]],
       return((1-rho)*b/(r+xi)*xi/(eta+xi))
  )
}


#' @title Setup Xinits.SIPd
#' @description Implements [setup_Xinits] for the SIPd model
#' @inheritParams ramp.xds::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SIPd = function(pars, H, i, Xopts=list()){
  pars$Xinits[[i]] = make_Xinits_SIPd(pars$nStrata[i], H, Xopts)
  return(pars)
}

#' @title Make initial values for the SIPd human model, with defaults
#' @param nStrata the number of population strata
#' @param H the initial human population density
#' @param Xopts a [list] that could overwrite defaults
#' @param I the initial values of the parameter I
#' @param P the initial values of the parameter P
#' @return a [list]
#' @export
make_Xinits_SIPd = function(nStrata, H, Xopts = list(),
                             I=1, P=0){with(Xopts,{
                               S = checkIt(H-I-P, nStrata)
                               I = checkIt(I, nStrata)
                               P = checkIt(P, nStrata)
                               return(list(S=S, I=I, P=P, cfoi=P*0))
  })}

#' @title Parse the output of deSolve and return variables for the SIPd model
#' @description Implements [parse_Xorbits] for the SIPd model
#' @inheritParams ramp.xds::parse_Xorbits
#' @return none
#' @export
parse_Xorbits.SIPd <- function(outputs, pars, i) {
  with(pars$ix$X[[i]],{
    S = outputs[,S_ix]
    I = outputs[,I_ix]
    P = outputs[,P_ix]
    H = S+I+P
    vars <- list(S=S, I=I, P=P, H=H)
    vars$ni <- F_ni(vars, pars$Xpar[[i]])
    vars$true_pr <- F_pr(vars, pars$Xpar[[i]])
    return(vars)
  })}

#' @title Add indices for human population to parameter list
#' @description Implements [setup_Xix] for the SIPd model.
#' @inheritParams ramp.xds::setup_Xix
#' @return none
#' @importFrom utils tail
#' @export
setup_Xix.SIPd <- function(pars, i) {with(pars,{

  S_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(S_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(I_ix, 1)

  P_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(P_ix, 1)

  cfoi_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(cfoi_ix, 1)

  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(S_ix=S_ix, I_ix=I_ix, P_ix=P_ix, cfoi_ix=cfoi_ix)
  return(pars)
})}


#' @title Update inits for the SIPd human model from a vector of states
#' @inheritParams ramp.xds::update_Xinits
#' @return none
#' @export
update_Xinits.SIPd <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars$Xinits[[i]] = make_Xinits_SIPd(pars$nStrata[i], H, list(), I=I, P=P)
    return(pars)
  })}


#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xds::get_Xinits
#' @return a named [vector]
#' @export
get_Xinits.SIPd <- function(pars, i){
  pars$Xinits[[i]]
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_ni] for the SIPd model.
#' @inheritParams ramp.xds::F_ni
#' @return a [numeric] vector of length `nStrata`
#' @export
F_ni.SIPd <- function(vars, Xpar) {
  ni = with(vars, Xpar$c*I/H)
  return(ni)
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SIPd model.
#' @inheritParams ramp.xds::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SIPd <- function(vars, Xpar) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute the prevalence of infection by light microscopy
#' @description Implements [F_pr] for the SIPd model.
#' @inheritParams ramp.xds::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr_by_lm.SIPd <- function(vars, Xpar) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute the prevalence of infection by RDT
#' @description Implements [F_pr] for the SIPd model.
#' @inheritParams ramp.xds::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr_by_rdt.SIPd <- function(vars, Xpar) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute the prevalence of infection by pcr
#' @description Implements [F_pr] for the SIPd model.
#' @inheritParams ramp.xds::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr_by_pcr.SIPd <- function(vars, Xpar) {
  pr = with(vars, I/H)
  return(pr)
}

#' Plot the density of infected individuals for the SIPd model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @export
xds_plot_X.SIPd = function(pars, i=1, clrs=c("darkblue", "darkred", "darkgreen"), llty=1, add=FALSE){
  XH = pars$outputs$orbits$XH[[i]]
  time = pars$outputs$time

  if(add==FALSE)
    plot(time, 0*time, type = "n", ylim = c(0, max(XH$H)),
         ylab = "# Infected", xlab = "Time")

  add_lines_X_SIS(time, XH, pars$nStrata[i], clrs, llty)
}


#' Add lines for the density of infected individuals for the SIPd model
#'
#' @param time time points for the observations
#' @param XH a list with the outputs of parse_Xorbits.SIPd
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
add_lines_X_SIPd = function(time, XH, nStrata, clrs=c("darkblue", "darkred", "darkgreen"), llty=1){
  if (length(llty)< nStrata) llty = rep(llty, nStrata)
  with(XH,{
    if(nStrata == 1){
      lines(time, S, col=clrs[1], lty = llty)
      lines(time, I, col=clrs[2], lty = llty)
      lines(time, P, col=clrs[3], lty = llty)
    } else {
      for(i in 1:nStrata)
        lines(time, S[,i], col=clrs[1], lty = llty[i])
      lines(time, I[,i], col=clrs[2], lty = llty[i])
      lines(time, P[,i], col=clrs[3], lty = llty[i])
    }})
}

#' @title Compute the steady states for the SIPd model as a function of the daily foi
#' @description Compute the steady state of the SIPd model as a function of the daily eir.
#' @inheritParams ramp.xds::xde_steady_state_X
#' @return the steady states as a named vector
#' @export
xde_steady_state_X.SIPd = function(foi, H, Xpar){with(Xpar,{
  Ieq = (foi*H*eta*(1-rho))/((foi+r+xi)*(eta+xi) +foi*(r-eta)*rho)
  Peq  = (H*xi*(foi+r+xi) + (foi*H*r*rho))/((foi+r+xi)*(eta+xi) +foi*(r-eta)*rho)
  Seq = H -Ieq - Peq
  return(list(S=as.vector(Seq), I=as.vector(Ieq), P = as.vector(Peq)))
})}

