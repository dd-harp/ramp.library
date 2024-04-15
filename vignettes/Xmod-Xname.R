#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIS model, no demography.
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.workhorse <- function(t, y, pars, i) {

  foi <- pars$FoI[[i]]

  with(list_Xvars(y, pars, i),{
    H <- F_H(t, y, pars, i)

    with(pars$Xpar[[i]], {
      peak = function(w){return(1)}
      fever1 = function(w){return(1)}
      trt_A = fever1(w)*zeta_A
      trt_1 = zeta_1 + zeta_U
      trt_2 = zeta_2 + zeta_U
      trt_3 = zeta_3 + zeta_U
      trt_4 = zeta_4 + zeta_U

      foi = F_foi(a-tau)

      AA = A0+A1+A2+A3+A4
      zA = trt_1*A1 + trt_2*A2 + trt_3*A3 + trt_4*A4

      dU  =  r*(I1+I2+I3+I4) + sigma*(P+G) - (foi + zeta_U)*U
      dA0 = foi*U + r*(A1+A2+A3+A4) - (phi + trt_A + zeta_U)*A0

      dP =  zeta_U*(U+A0) + zA + trt_2*I2 + trt_3*I3 + trt_4*I4 - sigma*P
      dG =  trt_1*I1 + trt_A*AA - sigma*G

      dI1 =  peak(w)*phi*AA               - (foi+r+xi_1+trt_1)*I1
      dI2 =  (1-peak(w))*phi*AA + xi_1*I1 - (foi+r+xi_2+trt_2)*I2
      dI3 =                       xi_2*I2 - (foi+r+xi_3+trt_3)*I3
      dI4 =                       xi_3*I3 - (foi+r     +trt_4)*I4

      dA1 = foi*I1 - (phi + r + trt_1 + trt_A)*A1
      dA2 = foi*I2 - (phi + r + trt_2 + trt_A)*A2
      dA3 = foi*I3 - (phi + r + trt_3 + trt_A)*A3
      dA4 = foi*I4 - (phi + r + trt_4 + trt_A)*A4

      dw = foi_tau

      return(c(dU, dP, dG, dI1, dI2, dI3, dI4, dA0, dA1, dA2, dA3, dA4, dw))
    })
  })
}





#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xde::list_Xvars
#' @return a [list]
#' @export
list_Xvars.workhorse <- function(y, pars, i) {
  vars = with(pars$ix$X[[i]], list(
        U =  y[U_ix],
        A0 = y[A0_ix],
        P =  y[P_ix],
        G =  y[G_ix],
        I1 = y[I1_ix],
        I2 = y[I2_ix],
        I3 = y[I3_ix],
        I4 = y[I4_ix],
        A1 = y[A1_ix],
        A2 = y[A2_ix],
        A3 = y[A3_ix],
        A4 = y[A4_ix],
        w =  y[w_ix]
  ))
  return(vars)
}


#' @title Make parameters for workhorse human model, with defaults
#' @param pars a [list]
#' @param Xopts a [list] that could overwrite defaults
#' @param tau the incubation period
#' @param r the clearance rate
#' @param sigma the rate individuals lose chemoprotection
#' @param phi the rate of transition from the acute to the chronic phase
#' @param xi_1 the rate of transition from stage 1 to stage 2
#' @param xi_2 the rate of transition from stage 2 to stage 3
#' @param xi_3 the rate of transition from stage 3 to stage 4
#' @param zeta_U the treatment rate for uninfected individuals
#' @param zeta_A the treatment rate for acutely infected individuals
#' @param zeta_1 the treatment rate for stage 1
#' @param zeta_2 the treatment rate for stage 2
#' @param zeta_3 the treatment rate for stage 3
#' @param zeta_4 the treatment rate for stage 4
#' @return a [list]
#' @export
make_Xpar_workhorse = function(pars, Xopts=list(),
                               tau = 12,
                               b = 0.55,
                               r = 1/200,
                               sigma = 1/30,
                               phi = 1/10,
                               xi_1 = 1/20,
                               xi_2 = 1/40,
                               xi_3 = 1/80,
                               zeta_U = 1/365,
                               zeta_A = 1/3,
                               zeta_1 = 1/30,
                               zeta_2 = 1/180,
                               zeta_3 = 1/180,
                               zeta_4 = 1/365){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("workhorse")

    Xpar$tau = checkIt(tau, pars$nStrata)
    Xpar$b = checkIt(b, pars$nStrata)
    Xpar$r = checkIt(r, pars$nStrata)
    Xpar$sigma = checkIt(sigma, pars$nStrata)
    Xpar$phi = checkIt(phi, pars$nStrata)
    Xpar$xi_1 = checkIt(xi_1, pars$nStrata)
    Xpar$xi_2 = checkIt(xi_2, pars$nStrata)
    Xpar$xi_3 = checkIt(xi_3, pars$nStrata)
    Xpar$zeta_U = checkIt(zeta_U, pars$nStrata)
    Xpar$zeta_A = checkIt(zeta_A, pars$nStrata)
    Xpar$zeta_1 = checkIt(zeta_1, pars$nStrata)
    Xpar$zeta_2 = checkIt(zeta_2, pars$nStrata)
    Xpar$zeta_3 = checkIt(zeta_3, pars$nStrata)
    Xpar$zeta_4 = checkIt(zeta_4, pars$nStrata)

    pars$Xpar = Xpar
    return(pars)
  })}



#' @title Setup Xpar.workhorse
#' @description Implements [setup_Xpar] for the workhorse model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.workhorse = function(workhorse, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_workhorse(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}





#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.workhorse <- function(t, y, pars, i) {
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



#' @title Size of effective infectious human population
#' @description Implements [F_H] for the workhorse model.
#' @inheritParams ramp.xde::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.workhorse <- function(t, y, pars, i){
  with(list_Xpars(y, pars, i),{
    H <- X1 + X2 + ...
    return(H)
})}





#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the workhorse model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.workhorse <- function(y, pars, i) {
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





#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the workhorse model.
#' @inheritParams ramp.xde::make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.workhorse <- function(pars, i) {with(pars,{

  X1_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(X1_ix, 1)

  X2_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(X2_ix, 1)

  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(X1_ix=X1_ix, X2_ix=X2_ix, ...)
  return(pars)
})}



#' @title Make initial values for the workhorse human model, with defaults
#' @param nStrata the number of strata in the model
#' @param Xopts a [list] to overwrite defaults
#' @param X10 the initial value for X1
#' @param X20 the initial value for X1
#' @return a [list]
#' @export
make_Xinits_workhorse = function(nStrata, Xopts = list(), X10=NULL, X20=1){with(Xopts,{
  stopifnot(is.numeric(X10))
  stopifnot(is.numeric(X20))
  X1 = checkIt(X10, nStrata)
  X2 = checkIt(X20, nStrata)
  return(list(X1=X1, X2=X2, ...))
})}





#' @title Setup Xinits.workhorse
#' @description Implements [setup_Xinits] for the workhorse model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.workhorse = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_workhorse(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}





#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.workhorse <- function(pars, i){
  with(pars$Xinits[[i]], return(c(X1,X2)))
}





#' @title Update inits for the workhorse human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.workhorse <- function(pars, y0, i) {
  with(pars$ix$X[[i]],{
    X10 = y0[X1_ix]
    X20 = y0[X2_ix]
    ...
    pars = make_Xinits_workhorse(pars, list(), X10, X20, ...)
    return(pars)
})}





#' @title Parse the output of deSolve and return variables for the workhorse model
#' @description Implements [parse_deout_X] for the workhorse model
#' @inheritParams ramp.xde::parse_deout_X
#' @return none
#' @export
parse_deout_X.workhorse <- function(deout, pars, i) {
  time = deout[,1]
  with(pars$ix$X[[i]],{
    X1 = deout[,X1_ix+1]
    X2 = deout[,X2_ix+1]
    ...
    H = X1 + X2 + ...
    return(list(time=time, X1=X1, X2=X2, ..., H=H))
})}





#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the workhorse model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.workhorse<- function(varslist, pars) {
  pr = compute_pr_formula
  return(pr)
}





#' @title Compute the HTC for the workhorse model
#' @description Implements [HTC] for the workhorse model with demography.
#' @inheritParams ramp.xde::HTC
#' @return a [numeric] vector
#' @export
HTC.workhorse <- function(pars, i) {
  with(pars$Xpar[[i]],
    #HTC <-
    return(HTC)
  )
}



#' Add lines for the density of infected individuals for the workhorse model
#'
#' @param XH a list with the outputs of parse_deout_X_workhorse
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
xde_lines_X_workhorse = function(XH, nStrata, clrs=c("darkblue","darkred"), llty=1){
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





#' Plot the density of infected individuals for the workhorse model
#'
#' @inheritParams ramp.xde::xde_plot_X
#' @export
xde_plot_X.workhorse = function(pars, i=1, clrs=c("darkblue","darkred"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "Time"))


  xde_lines_X_workhorse(vars$XH[[i]], pars$Hpar[[i]]$nStrata, clrs, llty)
}
