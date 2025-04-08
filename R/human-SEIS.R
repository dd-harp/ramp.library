# The SEIS

#' @title \eqn{\cal X} Component Derivatives for the `SEIS` Model
#' @description Implements [dXdt] for the SEIS model
#' @inheritParams ramp.xds::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SEIS <- function(t, y, pars, i) {

  foi <- pars$FoI[[i]]
  Hpar <- pars$Hpar[[i]]
  with(list_Xvars(y, pars, i),{
    with(pars$Xpar[[i]], {
      dH <- Births(t, H, Hpar) + dHdt(t, H, Hpar)
      dE <- foi*S - nu*E + dHdt(t, E, Hpar)
      dI <- nu*E - r*I + dHdt(t, I, Hpar)
      return(c(dH, dE, dI))
    })
  })
}

#' @title DTS updating for the SEIS model for human / vertebrate host infections
#' @description Implements [Update_Xt] for the SEIS model
#' @inheritParams ramp.xds::Update_Xt
#' @return a [numeric] vector
#' @export
Update_Xt.SEIS <- function(t, y, pars, i) {

  ar <- pars$AR[[i]]
  Hpar <- pars$Hpar[[i]]
  with(list_Xvars(y, pars, i),{
    with(pars$Xpar[[i]], {

      Ht <- dHdt(t, H, Hpar) + Births(t, H, Hpar)
      Et <- ar*S + (1-nr)*ar*I + nu*E + dHdt(t, E, Hpar)
      It <- nr*I + (1-nu)*E + dHdt(t, I, Hpar)

      return(c(H=unname(Ht), E=unname(Et), I=unname(It)))
    })
  })
}

#' @title Compute the steady states for the  dts SEIS model as a function of the daily EIR
#' @description Compute the steady state of the  dts SIS model as a function of the daily eir.
#' @inheritParams ramp.xds::dts_steady_state_X
#' @return the steady states as a named vector
#' @export
dts_steady_state_X.SEIS = function(ar, H, Xpar){with(Xpar,{
  Iteq = (ar*H*(1-nu))/((1-nr)+ar-nu+ (nr*nu*(1-ar)))
  Eteq = (ar*H*(1-nr))/((1-nr)+ar-nu+ (nr*nu*(1-ar)))
  Steq = H-Iteq-Eteq

  return(c(S=Steq, E=Eteq, I=Iteq))
})}

#' @title Setup Xpar.SEIS
#' @description Implements [setup_Xpar] for the SEIS model
#' @inheritParams ramp.xds::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SEIS = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SEIS(pars$nStrata[i], Xopts)
  return(pars)
}

#' @title Make parameters for SEIS xde human model, with defaults
#' @param nStrata is the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param b transmission probability (efficiency) from mosquito to human
#' @param r recovery rate
#' @param nu 1/latent period
#' @param c transmission probability (efficiency) from human to mosquito
#' @return a [list]
#' @export
make_Xpar_SEIS = function(nStrata, Xopts=list(),
                            b=0.55, r=1/180, nu=1/14, c=0.15){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- "SEIS"

    Xpar$b = checkIt(b, nStrata)
    Xpar$c = checkIt(c, nStrata)
    Xpar$r = checkIt(r, nStrata)
    Xpar$nu = checkIt(nu, nStrata)

    return(Xpar)
  })}


#' @title Return the parameters as a list
#' @description This method dispatches on the type of `pars$Xpar[[i]]`.
#' @inheritParams ramp.xds::set_Xpars
#' @return an **`xds`** object
#' @export
set_Xpars.SEIS <- function(pars, i=1, Xopts=list()) {
  with(pars$Xpar[[i]], with(Xopts,{
    pars$Xpar[[i]]$b <- b
    pars$Xpar[[i]]$c <- c
    pars$Xpar[[i]]$r <- r
    pars$Xpar[[i]]$nu <- nu
    return(pars)
  }))}


#' @title Return the parameters as a list
#' @description This method dispatches on the type of `pars$Xpar[[i]]`.
#' @inheritParams ramp.xds::set_Xinits
#' @return an **`xds`** object
#' @export
set_Xinits.SEIS <- function(pars, i=1, Xopts=list()) {
  with(get_Xinits(pars, i), with(Xopts,{
    pars$Xinits[[i]]$H = H
    pars$Xinits[[i]]$E = E
    pars$Xinits[[i]]$I = I
    pars$Xinits[[i]]$S = H-E-I
    return(pars)
  }))}

#' @title Return the parameters as a list
#' @description Parameter values for the \eqn{i^{th}} host are
#' stored as `pars$Xpar[[i]]`. This returns the stored parameter
#' values as a list.
#' @inheritParams ramp.xds::get_Xpars
#' @return a [list]
#' @seealso [make_Xpar_SEIS]
#' @export
get_Xpars.SEIS <- function(pars, i=1) {
  with(pars$Xpar[[i]],list(b=b, c=c, r=r, nu=nu))
}


# specialized methods for the human SEIS model

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SEIS model.
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SEIS <- function(t, y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  X = with(pars$Xpar[[i]], c*I)
  return(X)
}

#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SEIS model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SEIS <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i), return(H))
}


#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SEIS model.
#' @inheritParams ramp.xds::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SEIS <- function(y, pars, i) {
  with(pars$Xpar[[i]], b)
}

#' @title Make initial values for the SEIS xde human model, with defaults
#' @param nStrata the number of strata in the model
#' @param H the initial human population density
#' @param Xopts a [list] to overwrite defaults
#' @param E the initial values of the parameter E
#' @param I the initial values of the parameter I
#' @return a [list]
#' @export
make_Xinits_SEIS = function(nStrata, H, Xopts = list(), E=0, I=1){with(Xopts,{
  H = checkIt(H, nStrata)
  E = checkIt(E, nStrata)
  I = checkIt(I, nStrata)
  return(list(H=H, E=E, I=I))
})}




#' @title Return the SEIS model variables as a list, returned from Update_Xt.SISd
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xds::put_Xvars
#' @return a [list]
#' @export
put_Xvars.SEIS <- function(Xvars, y, pars, i) {
  with(pars$ix$X[[i]],
       with(as.list(Xvars),{
         y[H_ix] = H
         y[E_ix] = E
         y[I_ix] = I
         return(y)
       }))}

#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xds::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SEIS <- function(y, pars, i) {
  with(pars$ix$X[[i]],{
    H = y[H_ix]
    E = y[E_ix]
    I = y[I_ix]
    S = H-E-I
    return(list(S=S,E=E,I=I,H=H))})
}


#' @title Setup Xinits.SEIS
#' @description Implements [setup_Xinits] for the SEIS model
#' @inheritParams ramp.xds::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SEIS = function(pars, H, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_SEIS(pars$nStrata[i], H, Xopts))
  return(pars)
}

#' @title Add indices for human population to parameter list
#' @description Implements [setup_Xix] for the SEIS model.
#' @inheritParams ramp.xds::setup_Xix
#' @return none
#' @importFrom utils tail
#' @export
setup_Xix.SEIS <- function(pars, i) {with(pars,{

  H_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(H_ix, 1)

  E_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(E_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(I_ix, 1)

  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(H_ix=H_ix, E_ix=E_ix, I_ix=I_ix)
  return(pars)
})}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xds::get_Xinits
#' @return a [numeric] vector
#' @export
get_Xinits.SEIS <- function(pars, i=1){pars$Xinits[[i]]}

#' @title Update inits for the SEIS xde human model from a vector of states
#' @inheritParams ramp.xds::update_Xinits
#' @return none
#' @export
update_Xinits.SEIS <- function(pars, y, i) {
  with(list_Xvars(y, pars, i),{
    pars$Xinits[[i]] = make_Xinits_SEIS(pars$nStrata[i], H, list(), E=E, I=I)
    return(pars)
  })}


#' @title Parse the output of deSolve and return variables for the SEIS model
#' @description Implements [parse_Xorbits] for the SEIS model
#' @inheritParams ramp.xds::parse_Xorbits
#' @return none
#' @export
parse_Xorbits.SEIS <- function(outputs, pars, i) {
  with(pars$ix$X[[i]],{
    H = outputs[,H_ix]
    E = outputs[,E_ix]
    I = outputs[,I_ix]
    S = H-E-I
    vars <- list(S=S, E=E, I=I, H=H)
    vars$ni <- F_ni(vars, pars$Xpar[[i]])
    vars$pr <- F_pr(vars, pars$Xpar[[i]])
    return(vars)
  })}

#' @title Compute the HTC for the SEIS model
#' @description Implements [HTC] for the SEIS model with demography.
#' @inheritParams ramp.xds::HTC
#' @return a [numeric] vector
#' @export
HTC.SEIS <- function(pars, i) {
  with(pars$Xpar[[i]],
       return(c/r)
  )
}

#' @title Compute the NI
#' @description Implements [F_ni] for the SEIS model.
#' @inheritParams ramp.xds::F_ni
#' @return a [numeric] vector of length `nStrata`
#' @export
F_ni.SEIS <- function(vars, Xpar) {
  ni = with(vars, Xpar$c*I/H)
  return(ni)
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIS model.
#' @inheritParams ramp.xds::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIS <- function(vars, Xpar) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute the prevalence of infection by light microscopy
#' @description Implements [F_pr] for the SEIS model.
#' @inheritParams ramp.xds::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr_by_lm.SEIS <- function(vars, Xpar) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute the prevalence of infection by RDT
#' @description Implements [F_pr] for the SEIS model.
#' @inheritParams ramp.xds::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr_by_rdt.SEIS <- function(vars, Xpar) {
  pr = with(vars, I/H)
  return(pr)
}

#' @title Compute the prevalence of infection by pcr
#' @description Implements [F_pr] for the SEIS model.
#' @inheritParams ramp.xds::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr_by_pcr.SEIS <- function(vars, Xpar) {
  pr = with(vars, I/H)
  return(pr)
}


#' Plot the density of infected individuals for the SEIS model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @export
xds_plot_X.SEIS = function(pars, i=1, clrs=c("darkblue","darkred"), llty=1, add=FALSE){
  XH = pars$outputs$orbits$XH[[i]]
  time = pars$outputs$time

  if(add==FALSE)
    plot(time, 0*time, type = "n", ylim = c(0, max(XH$H)),
         ylab = "# Infected", xlab = "Time")

  add_lines_X_SIS(time, XH, pars$nStrata[i], clrs, llty)
}

#' Add lines for the density of infected individuals for the SEIS model
#'
#' @param time time points for the observations
#' @param XH a list with the outputs of parse_Xorbits_SEIS
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
add_lines_X_SEIS = function(time, XH, nStrata, clrs=c("darkblue","darkred"), llty=1){
  if (length(llty)< nStrata) llty = rep(llty, nStrata)
  with(XH,{
    if(nStrata == 1){
      lines(time, S+E, col=clrs[1], lty = llty)
      lines(time, I, col=clrs[2], lty = llty)
    } else {
      for(i in 1:nStrata)
        lines(time, S[,i] + E[,i], col=clrs[1], lty = llty[i])
      lines(time, I[,i], col=clrs[2], lty = llty[i])
    }})
}

#' @title Compute the steady states for the SEIS model as a function of the daily EIR
#' @description Compute the steady state of the SEIS model as a function of the daily eir.
#' @inheritParams ramp.xds::xde_steady_state_X
#' @return the steady states as a named vector
#' @export
xde_steady_state_X.SEIS = function(foi, H, Xpar){with(Xpar,{
  Ieq = (foi*H*nu)/(foi*(r+nu) +r*nu)
  Eeq = (foi*H*r)/(foi*(r+nu) +r*nu)
  Seq = H- Ieq-Eeq
  return(c(S=Seq, E=Eeq, I=Ieq, H=H))
})}

