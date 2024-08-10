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

      return(c(S=unname(St), I=unname(It), P=unname(Pt), w=unname(wt)))
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

#' @title Make parameters for SIPw_dts human model, with defaults
#' @param nStrata the number of population strata
#' @param D the operating time step
#' @param Xopts options to overwrite defaults
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @param rho probability of successful treatment upon infection
#' @param sigma probability of treatment for treated, above background
#' @param xi background treatment rate
#' @param eta rate of loss of chemo-protection
#' @return a [list]
#' @export
dts_make_Xpar_SIPw = function(nStrata, D=1, Xopts=list(), b=0.55, c=0.15, r=1/180,
                              rho=0.1, sigma=1/730, xi=1/365, eta=1/25){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SIPw")

    Xpar$D     = checkIt(D, 1)
    Xpar$b     = checkIt(b, nStrata)
    Xpar$c     = checkIt(c, nStrata)
    Xpar$r     = 1-exp(-checkIt(r, nStrata)*D)
    Xpar$rho   = checkIt(rho, nStrata)
    Xpar$xi    = 1-exp(-checkIt(xi, nStrata)*D)
    Xpar$xisig = 1-exp(-checkIt(sigma+xi, nStrata)*D)
    Xpar$eta   = 1-exp(-checkIt(eta, nStrata)*D)

    return(Xpar)
})}

#' @title Size of effective infectious human population
#' @description Implements [F_X] for SIPw models
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SIPw <- function(y, pars, i) {
  with(list_Xvars(y, pars, i),
       with(pars$Xpar[[i]],{
         X = c*I
         return(X)
       })
  )
}

#' @title Size of effective infectious human population
#' @description Implements [F_X] for SIPw models
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SIPw <- function(y, pars, i){
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
make_Xinits.SIPw = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = create_Xinits_SIPw(pars$nStrata[i], Xopts, H0=pars$Hpar[[i]]$H)
  return(pars)
}

#' @title Make initial values for a SIPw human model, with defaults
#' @param nStrata the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param H0 the initial human population density
#' @param S0 the initial values of the variable S
#' @param I0 the initial values of the variable I
#' @param P0 the initial values of the variable P
#' @param w0 the initial values of the tracking variable w
#' @return a [list]
#' @export
create_Xinits_SIPw = function(nStrata, Xopts = list(),
                            H0=NULL, S0=NULL, I0=1, P0=0, w0=0){
  with(Xopts,{
    if(is.null(S0)) S0=H0-I0-P0
    S = checkIt(S0, nStrata)
    I = checkIt(I0, nStrata)
    P = checkIt(P0, nStrata)
    w = checkIt(w0, nStrata)
    return(list(S=S, I=I, P=P, w=w))
})}

#' @title Parse the output of deSolve and return variables for SIPw models
#' @description Implements [parse_outputs_X] for SIPw models
#' @inheritParams ramp.xds::parse_outputs_X
#' @return none
#' @export
parse_outputs_X.SIPw <- function(outputs, pars, i) {
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
    pars$Xinits[[i]] = create_Xinits_SIPw(pars$nStrata[i], list(), S0=S, I0=I, P0=P, w0=w)
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
xds_plot_X.SIPw = function(pars, i=1, clrs=c("darkblue", "darkred", "darkgreen"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "Time"))

  add_lines_X_SIPw(vars$XH[[i]], pars, clrs, llty)
}


#' Add lines for the density of infected individuals for a SIPw model
#'
#' @param XH a list with the outputs of parse_outputs_X.SIPw
#' @param pars a list that defines an `ramp.xds` model (*e.g.*,  generated by `make()`)
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
add_lines_X_SIPw = function(XH, pars, clrs=c("darkblue", "darkred", "darkgreen"), llty=1){
  with(XH,{
    if(pars$Hpar[[1]]$nStrata==1) {
      lines(time, S, col=clrs[1], lty = llty[1])
      lines(time, I, col=clrs[2], lty = llty[1])
      lines(time, P, col=clrs[3], lty = llty[1])
    }
    if(pars$Hpar[[1]]$nStrata>1){
      if (length(clrs)==1) clrs=matrix(clrs, 3, pars$nStrata[i])
      if (length(llty)==1) llty=rep(llty, pars$nStrata[i])
      for(i in 1:pars$Hpar[[1]]$nStrata){
        lines(time, S[,i], col=clrs[1,i], lty = llty[i])
        lines(time, I[,i], col=clrs[2,i], lty = llty[i])
        lines(time, P[,i], col=clrs[3,i], lty = llty[i])
      }
    }
  })}

