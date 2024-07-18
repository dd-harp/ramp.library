#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIS model, no demography.
#' @inheritParams ramp.xds::dXdt
#' @return a [numeric] vector
#' @export
dXdt.workhorse <- function(t, y, pars, i) {

  foi <- pars$FoI[[i]]
  with(list_Xvars(y, pars, i), {

    with(pars$Xpar[[i]], {
      peak = function(w){return(1)}
      fever1 = function(w){return(1)}
      trt_A = fever1(w)*zeta_A
      trt_1 = zeta_1 + zeta_U
      trt_2 = zeta_2 + zeta_U
      trt_3 = zeta_3 + zeta_U
      trt_4 = zeta_4 + zeta_U



      if(t <= tau) {
        foi_tau = foi
      } else {
        foi_tau = lagderiv(t = t-tau, nr = pars$ix$X[[i]]$cfoi_ix)
      }


      AA = A0+A1+A2+A3+A4
      zA = trt_1*A1 + trt_2*A2 + trt_3*A3 + trt_4*A4

      dU  =  r*(I1+I2+I3+I4) + sigma*(P+G) - (foi_tau + zeta_U)*U
      dA0 = foi_tau*U + r*(A1+A2+A3+A4) - (phi + trt_A + zeta_U)*A0

      dP =  zeta_U*(U+A0) + zA + trt_2*I2 + trt_3*I3 + trt_4*I4 - sigma*P
      dG =  trt_1*I1 + trt_A*AA - sigma*G

      dI1 =  peak(w)*phi*AA               - (foi_tau+r+xi_1+trt_1)*I1
      dI2 =  (1-peak(w))*phi*AA + xi_1*I1 - (foi_tau+r+xi_2+trt_2)*I2
      dI3 =                       xi_2*I2 - (foi_tau+r+xi_3+trt_3)*I3
      dI4 =                       xi_3*I3 - (foi_tau+r     +trt_4)*I4

      dA1 = foi_tau*I1 - (phi + r + trt_1 + trt_A)*A1
      dA2 = foi_tau*I2 - (phi + r + trt_2 + trt_A)*A2
      dA3 = foi_tau*I3 - (phi + r + trt_3 + trt_A)*A3
      dA4 = foi_tau*I4 - (phi + r + trt_4 + trt_A)*A4

      dw = foi_tau
      #dH = sum(c(dU, dA0, dP, dG, dI1, dI2, dI3, dI4, dA1, dA2, dA3, dA4))
      #browser()
      return(c(dU, dA0, dP, dG, dI1, dI2, dI3, dI4, dA1, dA2, dA3, dA4, dw, foi))
    })
  })
}


#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xds::list_Xvars
#' @return a [list]
#' @export
list_Xvars.workhorse <- function(y, pars, i) {
  vars = with(pars$ix$X[[i]],{
    U =  y[U_ix]
    A0 = y[A0_ix]
    P =  y[P_ix]
    G =  y[G_ix]
    I1 = y[I1_ix]
    I2 = y[I2_ix]
    I3 = y[I3_ix]
    I4 = y[I4_ix]
    A1 = y[A1_ix]
    A2 = y[A2_ix]
    A3 = y[A3_ix]
    A4 = y[A4_ix]
    w =  y[w_ix]
    H <- U + A0 + P + G + I1 + I2 + I3 + I4 + A1 + A2 + A3 + A4
    return(list(U=U,A0=A0,P=P,G=G,I1=I1,I2=I2,I3=I3,I4=I4,A1=A1,A2=A2,A3=A3,A4=A4,w=w,H=H))})
}



#' @title Make parameters for workhorse human model, with defaults
#' @param nStrata the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param tau the incubation period
#' @param b the probability of infection, per infectious bite
#' @param r the clearance rate
#' @param sigma the rate individuals lose chemoprotection
#' @param phi the rate of transition from the acute to the chronic phase
#' @param c1 the probability of infecting a mosquito, stage 1
#' @param c2 the probability of infecting a mosquito, stage 2
#' @param c3 the probability of infecting a mosquito, stage 3
#' @param c4 the probability of infecting a mosquito, stage 4
#' @param cG the probability of infecting a mosquito, prophylaxed
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
make_Xpar_workhorse = function(nStrata, Xopts=list(),
                               tau = 12,
                               b = 0.55,
                               r = 1/200,
                               sigma = 1/30, phi = 1/10,
                               c1 = .3, c2 =.1, c3 = 0.03, c4 = 0.01, cG=0.15,
                               xi_1 = 1/20, xi_2 = 1/40, xi_3 = 1/80,
                               zeta_U = 1/365,
                               zeta_A = 1/100,
                               zeta_1 = 1/100,
                               zeta_2 = 1/180,
                               zeta_3 = 1/180,
                               zeta_4 = 1/365){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("workhorse")

    Xpar$tau = checkIt(tau, nStrata)
    Xpar$b = checkIt(b, nStrata)
    Xpar$r = checkIt(r, nStrata)
    Xpar$sigma = checkIt(sigma, nStrata)
    Xpar$c1 = checkIt(c1, nStrata)
    Xpar$c2 = checkIt(c2, nStrata)
    Xpar$c3 = checkIt(c3, nStrata)
    Xpar$c4 = checkIt(c4, nStrata)
    Xpar$cG = checkIt(cG, nStrata)
    Xpar$phi = checkIt(phi, nStrata)
    Xpar$xi_1 = checkIt(xi_1, nStrata)
    Xpar$xi_2 = checkIt(xi_2, nStrata)
    Xpar$xi_3 = checkIt(xi_3, nStrata)
    Xpar$zeta_U = checkIt(zeta_U, nStrata)
    Xpar$zeta_A = checkIt(zeta_A, nStrata)
    Xpar$zeta_1 = checkIt(zeta_1, nStrata)
    Xpar$zeta_2 = checkIt(zeta_2, nStrata)
    Xpar$zeta_3 = checkIt(zeta_3, nStrata)
    Xpar$zeta_4 = checkIt(zeta_4, nStrata)

    return(Xpar)
  })}

#' @title xde_setup Xpar.workhorse
#' @description Implements [xde_setup_Xpar] for the workhorse model
#' @inheritParams ramp.xds::xde_setup_Xpar
#' @return a [list] vector
#' @export
xde_setup_Xpar.workhorse = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_workhorse(pars$Hpar[[i]]$nStrata, Xopts)
  pars$xde = "dde"
  return(pars)
}

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.workhorse <- function(y, pars, i) {
  with(list_Xvars(y, pars,i),
  with(pars$Xpar[[i]],{
    X = c1*(A1+I1) + c2*(A2+I2) + c3*(A3+I3) + c4*(A4+I4) + cG*G
    return(X)
}))}


#' @title Size of effective infectious human population
#' @description Implements [F_H] for the workhorse model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.workhorse <- function(y, pars, i){
  with(list_Xvars(y, pars, i), return(H))
}


#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the workhorse model.
#' @inheritParams ramp.xds::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.workhorse <- function(y, pars, i) {
  with(pars$Xpar[[i]],{
    b = pars$Xpar[[i]]$b
    return(b)
})}

#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the workhorse model.
#' @inheritParams ramp.xds::make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.workhorse <- function(pars, i) {with(pars,{

  U_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(U_ix, 1)

  A0_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(A0_ix, 1)

  P_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(P_ix, 1)

  G_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(G_ix, 1)

  I1_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(I1_ix, 1)

  I2_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(I2_ix, 1)

  I3_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(I3_ix, 1)

  I4_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(I4_ix, 1)

  A1_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(A1_ix, 1)

  A2_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(A2_ix, 1)

  A3_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(A3_ix, 1)

  A4_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(A4_ix, 1)

  w_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(w_ix, 1)

  cfoi_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(cfoi_ix, 1)

  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(U_ix=U_ix,  A0_ix=A0_ix,  P_ix=P_ix,  G_ix=G_ix,
                        I1_ix=I1_ix,  I2_ix=I2_ix,  I3_ix=I3_ix,  I4_ix=I4_ix,
                        A1_ix=A1_ix,  A2_ix=A2_ix,  A3_ix=A3_ix,  A4_ix=A4_ix,
                        w_ix=w_ix,    cfoi_ix=cfoi_ix)
  return(pars)
})}



#' @title Make initial values for the workhorse human model, with defaults
#' @param nStrata the number of strata in the model
#' @param Xopts a [list] to overwrite defaults
#' @param H0 the initial population size
#' @param U the initial number uninfected
#' @param A0 the initial number in U and A
#' @param P the initial number chemoprotected and not infectious
#' @param G the initial number chemoprotected and infectious
#' @param I1 the initial number with the youngest chronic phase infection in stage 1
#' @param I2 the initial number with the youngest chronic phase infection in stage 2
#' @param I3 the initial number with the youngest chronic phase infection in stage 3
#' @param I4 the initial number with the youngest chronic phase infection in stage $
#' @param A1 the initial number in I1 and A
#' @param A2 the initial number in I2 and A
#' @param A3 the initial number in I3 and A
#' @param A4 the initial number in I4 and A
#' @param w exposure tracking variable
#' @return a [list]
#' @export
make_Xinits_workhorse = function(nStrata, Xopts = list(), H0 = NULL,
                                 U = 180,  A0 = 13, P = 11,  G = 3,
                                 I1 = 20, I2 = 40, I3 = 80, I4 = 100,
                                 A1 = 2, A2 = 4, A3 = 8, A4 = 12,
                                 w = 0
                                 ){
  with(Xopts,{
    stopifnot(is.numeric(c(U, A0, P, G, I1, I2, I3, I4, A1, A2, A3, A4, w)))
    U = checkIt(U, nStrata)
    A0 = checkIt(A0, nStrata)
    P = checkIt(P, nStrata)
    G = checkIt(G, nStrata)
    I1 = checkIt(I1, nStrata)
    I2 = checkIt(I2, nStrata)
    I3 = checkIt(I3, nStrata)
    I4 = checkIt(I4, nStrata)
    A1 = checkIt(A1, nStrata)
    A2 = checkIt(A2, nStrata)
    A3 = checkIt(A3, nStrata)
    A4 = checkIt(A4, nStrata)
    w = checkIt(w, nStrata)
    return(list(
         U=U,    A0=A0,  P=P,    G=G,
         I1=I1,  I2=I2,  I3=I3,  I4=I4,
         A1=A1,  A2=A2,  A3=A3,  A4=A4,
         w=w, rep(0, nStrata)))
})}





#' @title Setup Xinits.workhorse
#' @description Implements [setup_Xinits] for the workhorse model
#' @inheritParams ramp.xds::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.workhorse = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_workhorse(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}


#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xds::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.workhorse <- function(pars, i){
  with(pars$Xinits[[i]], return(c(U, A0, P, G, I1, I2, I3, I4, A1, A2, A3, A4, w, 0*w)))
}


#' @title Update inits for the workhorse human model from a vector of states
#' @inheritParams ramp.xds::update_inits_X
#' @return none
#' @export
update_inits_X.workhorse <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_workhorse(pars, list(), list(
      U=U,    A0=A0,  P=P,    G=G,
      I1=I1,  I2=I2,  I3=I3,  I4=I4,
      A1=A1,  A2=A2,  A3=A3,  A4=A4,
      w=w, cfoi = 0*w))
    return(pars)
})}


#' @title Parse the output of deSolve and return variables for the workhorse model
#' @description Implements [parse_outputs_X] for the workhorse model
#' @inheritParams ramp.xds::parse_outputs_X
#' @return none
#' @export
parse_outputs_X.workhorse <- function(outputs, pars, i) {
  time = outputs[,1]
  with(pars$ix$X[[i]],{
    U = outputs[,U_ix+1]
    A0 = outputs[,A0_ix+1]
    P = outputs[,P_ix+1]
    G = outputs[,G_ix+1]
    I1 = outputs[,I1_ix+1]
    I2 = outputs[,I2_ix+1]
    I3 = outputs[,I3_ix+1]
    I4 = outputs[,I4_ix+1]
    A1 = outputs[,A1_ix+1]
    A2 = outputs[,A2_ix+1]
    A3 = outputs[,A3_ix+1]
    A4 = outputs[,A4_ix+1]
    w = outputs[,w_ix+1]
    H = U+A0+P+G+I1+I2+I3+I4+A1+A2+A3+A4
    return(list(time=time, U=U, A0=A0, P=P, G=G,
                I1=I1, I2=I2, I3=I3, I4=I4,
                A1=A1, A2=A2, A3=A3, A4=A4,
                w=w, H=H))
})}





#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the workhorse model.
#' @inheritParams ramp.xds::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.workhorse <- function(vars, Xpar) {
  pr = with(vars, (A0+A1+A2+A3+A4+I1+I2+I3+I4)/H)
  return(pr)
}

#' @title Compute the HTC for the workhorse model
#' @description Implements [HTC] for the workhorse model with demography.
#' @inheritParams ramp.xds::HTC
#' @return a [numeric] vector
#' @export
HTC.workhorse <- function(pars, i) {
  HTC <- with(pars$Xpar[[i]],
            c1*I1/(r+xi_1) + xi_1/(r+xi_1)*c2/(r+xi_2) + xi_2/(r+xi_2)*c3/(r+xi_3) + xi_3/(r+xi_3)*c4/r
  )
    return(HTC)
}



#' Add lines for the density of infected individuals for the workhorse model
#'
#' @param XH a list with the outputs of parse_outputs_X_workhorse
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
xde_lines_X_workhorse = function(XH, nStrata, clrs=c("darkblue","darkgreen", "darkred", "purple"), llty=1){
  with(XH,{
    PG = P+G
    I = A0+A1+A2+A3+A4+I1+I2+I3+I4
    A = A0+A1+A2+A3+A4
    if(nStrata==1) {
      lines(time, U, col=clrs[1], lty = llty[1])
      lines(time, PG, col=clrs[2], lty = llty[1])
      lines(time, A, col=clrs[3], lty = llty[1])
      lines(time, I, col=clrs[4], lty = llty[1])
    }
    if(nStrata>1){
      if (length(clrs)==2) clrs=matrix(clrs, 2, nStrata)
      if (length(llty)==1) llty=rep(llty, nStrata)

      for(i in 1:nStrata){
        lines(time, U[,i], col=clrs[1,i], lty = llty[i])
        lines(time, PG[,i], col=clrs[2,i], lty = llty[i])
        lines(time, A[,i], col=clrs[3,i], lty = llty[i])
        lines(time, I[,i], col=clrs[4,i], lty = llty[i])
      }
    }
  })}





#' Plot the density of infected individuals for the workhorse model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @export
xds_plot_X.workhorse = function(pars, i=1, clrs=c("darkblue","darkgreen", "darkred", "purple"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "Time"))


  xde_lines_X_workhorse(vars$XH[[i]], pars$Hpar[[i]]$nStrata, clrs, llty)
}
