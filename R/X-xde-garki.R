# specialized methods for the human garki_xde model

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the garki_xde model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.garki_xde <- function(t, y, pars, i){
  y1 <- y[pars$Xpar[[i]]$y1_ix]
  X = with(pars$Xpar[[i]], y1)
  return(X)
}

#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SIS model.
#' @inheritParams ramp.xde::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.garki_xde <- function(t, y, pars, i){
  with(list_Xvars(y, pars,i),{
    H = x1+x2+x3+x4+y1+y2+y3
    return(H)
  })
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the garki_xde model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.garki_xde <- function(y, pars, i) {
  with(pars$Xpar[[i]], b)
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the garki_xde model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.garki_xde <- function(varslist, pars, i) {
#  pr = with(pars$Xpar,with(varslist$XH, (q1*y1+q2*y2+q3*y3)/H))
  pr = with(pars$Xpar,with(varslist$XH[[i]], (y1+y2+y3)/H))
  return(pr)
}
#' @title Derivatives for human population
#' @description Implements [dXdt] for the garki_xde model.
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.garki_xde = function(t, y, pars, i){

  foi <- pars$FoI[[i]]

  with(list_Xvars(y, pars, i),{

    with(pars$Xpar[[i]],{
      R1 = foi/(exp(foi/r1) - 1)
      R2 = foi/(exp(foi/r2) - 1)

      #dx1 = Births(t, x1, pars) - FoI*x1 + R1*y2 + dHdt(t, x1, pars)
      dx1 = - foi*x1 + R1*y2 + dHdt(t, x1, pars)
      dx2 = foi*x1 - nu*x2 + dHdt(t, x2, pars)
      dy1 = nu*x2 - alpha1*y1  + dHdt(t, y1, pars)
      dy2 = alpha1*y1 - R1*y2 - alpha2*y2 + dHdt(t, y2, pars)
      dy3 = alpha2*y2 + nu*x4 - R2*y3 + dHdt(t, y3, pars)
      dx3 = R2*y3 - foi*x3 + dHdt(t, x3, pars)
      dx4 = foi*x3 - nu*x4 + dHdt(t, x4, pars)
      #dH =  dHdt(t, H, pars)

      return(c(dx1, dx2, dy1, dy2, dy3, dx3, dx4, dH))
    })
  })
}

#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xde::list_Xvars
#' @return a [list]
#' @export
list_Xvars.garki_xde <- function(y, pars, i) {
  with(pars$ix$X[[i]],
       return(list(
         x1 <- y[x1_ix],
         x2 <- y[x2_ix],
         x3 <- y[x3_ix],
         x4 <- y[x4_ix],
         y1 <- y[y1_ix],
         y2 <- y[y2_ix],
         y3 <- y[y3_ix]
       )))
}

#' @title Setup Xpar.garki_xde
#' @description Implements [setup_Xpar] for the garki_xde model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.garki_xde = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_garki_xde(pars$nStrata, Xopts)
  return(pars)
}

#' @title Setup Xinits.garki_xde
#' @description Implements [setup_Xinits] for the garki_xde model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.garki_xde = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars,make_Xinits_garki_xde(nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}

#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the garki_xde model.
#' @inheritParams ramp.xde::make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.garki_xde <- function(pars, i) {with(pars,{

  x1_ix <- seq(from = max_ix+1, length.out=nx1trata)
  max_ix <- tail(x1_ix, 1)

  x2_ix <- seq(from = max_ix+1, length.out=nx2trata)
  max_ix <- tail(x2_ix, 1)

  y1_ix <- seq(from = max_ix+1, length.out=ny1trata)
  max_ix <- tail(y1_ix, 1)

  y2_ix <- seq(from = max_ix+1, length.out=ny2trata)
  max_ix <- tail(y2_ix, 1)

  y3_ix <- seq(from = max_ix+1, length.out=ny3trata)
  max_ix <- tail(y3_ix, 1)

  x3_ix <- seq(from = max_ix+1, length.out=nx3trata)
  max_ix <- tail(x3_ix, 1)

  x4_ix <- seq(from = max_ix+1, length.out=nx4trata)
  max_ix <- tail(x4_ix, 1)

  pars$max_ix <- max_ix
  pars$ix$X[[i]] = list(x1_ix = x1_ix,  x2_ix = x2_ix,
                        y1_iy = y1_iy,  y2_iy = y2_iy,  y3_iy = y3_iy,
                        x3_ix = x3_ix,  x4_ix = x4_ix)

  return(pars)
})}

#' @title Make parameters for garki_xde human model
#' @param nStrata is the number of population strata
#' @param Xopts an [list]
#' @param b transmission probability (efficiency) from mosquito to human
#' @param r1 a [numeric] recovery rate for non-immunes
#' @param r2 a [numeric]recovery rate for partially immmunes
#' @param nu a [numeric]incubation period
#' @param alpha1 a [numeric] rate of losing infectivity
#' @param alpha2 a [numeric] rate of acquiring immunity
#' @param q1 a [numeric] detection of y1
#' @param q2 a [numeric] detection of y2
#' @param q3 a [numeric] detection of y3
#' @param mu a [numeric] the death rate
#' @return a [list]
#' @export
make_Xpar_garki_xde = function(nStrata, Xopts=list(), b=0.55,
                           r1=.0023, r2=.023, nu=1/15,
                           alpha1=.002, alpha2=.00019,
                           q1=.7, q2=0.5, q3=0.3, mu=1/65/365){
  with(Xopts,{
    xde <- 'ode'
    class(xde) <- 'ode'

    garki_xde = list()
    class(garki_xde) <- "garki_xde"
    garki_xde$xde <- xde
    garki_xde$b=checkIt(b, nStrata)
    garki_xde$r1=checkIt(r1, nStrata)
    garki_xde$r2=checkIt(r2, nStrata)
    garki_xde$nu=checkIt(nu, nStrata)
    garki_xde$alpha1=checkIt(alpha1, nStrata)
    garki_xde$alpha2=checkIt(alpha2, nStrata)
    garki_xde$mu=checkIt(mu, nStrata)
    garki_xde$q1=checkIt(q1, nStrata)
    garki_xde$q2=checkIt(q2, nStrata)
    garki_xde$q3=checkIt(q3, nStrata)

    return(garki_xde)
  })}

#' @title Make inits for garki_xde human model. Note that the variables should sum up to H, so the initial value of x1 is not set. The values are passed in the same order as they are presented in the original paper.
#' @param nStrata is the number of population strata
#' @param Xopts a [list] with values to override default values
#' @param H0 a [numeric] initial value for total human population density
#' @param x1 a [numeric] initial value for the variable x1
#' @param x2 a [numeric] initial value for the variable x2
#' @param y1 a [numeric] initial value for the variable y1
#' @param y2 a [numeric] initial value for the variable y2
#' @param y3 a [numeric] initial value for the variable y3
#' @param x3 a [numeric] initial value for the variable x3
#' @param x4 a [numeric] initial value for the variable x4
#' @return none
#' @export
make_Xinits_garki_xde <- function(nStrata, Xopts = list(), H0=NULL, x1=NULL, x2=0, y1=0, y2=0, y3=0, x3=0, x4=0) {
  stopifnot(is.numeric(x2))
  stopifnot(is.numeric(y1))
  stopifnot(is.numeric(y2))
  stopifnot(is.numeric(y3))
  stopifnot(is.numeric(x3))
  stopifnot(is.numeric(x4))
  if(is.null(x1)) x1 = H0 - x2 - y1 - y2 - y3 - x3 - x4
  stopifnot(x1>0)

  x1 = checkIt(x1, nStrata)
  x2 = checkIt(x2, nStrata)
  y1 = checkIt(y1, nStrata)
  y2 = checkIt(y2, nStrata)
  y3 = checkIt(y3, nStrata)
  x3 = checkIt(x3, nStrata)
  x4 = checkIt(x4, nStrata)
  return(list(x1=x1,x2=x2,y1=y1,y2=y2,y3=y3,x3=x3,x4=x4))
}

#' @title Return initial values as a vector for the garki_xde model
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a named [list]
#' @export
get_inits_X.garki_xde <- function(pars, i){
  with(pars$Xinits[[i]], c(x1, x2, y1, y2, y3, x3, x4))
}

#' @title Update Xinits for the garki_xde model
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.garki_xde <- function(pars, y0, i){
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_garki_xde(pars, x1=x1, x2=x2, y1=y1, y2=y2, y3=y3, x3=x3, x4=x4)
    return(pars)
})}

#' Plot the density of infected individuals for the garki_xde model
#'
#' @inheritParams ramp.xde::xds_plot_X
#' @export
xds_plot_X.garki_xde = function(pars, i, clrs=viridisLite::turbo(7), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "Time"))

  xde_lines_X_garki_xde(vars$XH, pars, clrs, llty)
}

#' Add lines for the density of infected individuals for the garki_xde model
#'
#' @param XH a list with the outputs of parse_outputs_X_SIS
#' @param pars a list that defines an `ramp.xde` model (*e.g.*,  generated by `xde_setup()`)
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
xde_lines_X_garki_xde= function(XH, pars, clrs=viridisLite::turbo(7), llty=1){
  with(XH,{
    if(pars$nStrata==1){
      lines(time, x1, col=clrs[1], lty = llty[1])
      lines(time, x2, col=clrs[2], lty = llty[1])
      lines(time, y1, col=clrs[3], lty = llty[1])
      lines(time, y2, col=clrs[4], lty = llty[1])
      lines(time, y3, col=clrs[5], lty = llty[1])
      lines(time, x3, col=clrs[6], lty = llty[1])
      lines(time, x4, col=clrs[7], lty = llty[1])
    }
    if(pars$nStrata>1){
      if (length(clrs)==1) clrs=rep(clrs, pars$nStrata)
      if (length(llty)==1) llty=rep(llty, pars$nStrata)
      for(i in 1:pars$nStrata){
        lines(time, x1[,i], col=clrs[1], lty = llty[i])
        lines(time, x2[,i], col=clrs[2], lty = llty[i])
        lines(time, y1[,i], col=clrs[3], lty = llty[i])
        lines(time, y2[,i], col=clrs[4], lty = llty[i])
        lines(time, y3[,i], col=clrs[5], lty = llty[i])
        lines(time, x3[,i], col=clrs[6], lty = llty[i])
        lines(time, x4[,i], col=clrs[7], lty = llty[i])
      }
    }
  })}

#' @title Parse the output of deSolve and return variables for the garki_xde model
#' @description Implements [parse_outputs_X] for the garki_xde model
#' @inheritParams ramp.xde::parse_outputs_X
#' @return none
#' @export
parse_outputs_X.garki_xde <- function(outputs, pars, i) {
  time = outputs[,1]
  with(pars$ix$X[[i]],{
    x1 = outputs[,x1_ix+1]
    x2 = outputs[,x2_ix+1]
    y1 = outputs[,y1_ix+1]
    y2 = outputs[,y2_ix+1]
    y3 = outputs[,y3_ix+1]
    x3 = outputs[,x3_ix+1]
    x4 = outputs[,x4_ix+1]
    return(list(time=time, x1=x1, x2=x2, y1=y1, y2=y2, y3=y3, x3=x3, x4=x4, H=H))
  })}

#' @title Compute the HTC for the garki_xde model
#' @description Implements [HTC] for the garki_xde model
#' @inheritParams ramp.xde::HTC
#' @return a [numeric] vector
#' @export
HTC.garki_xde <- function(pars, i) {
  with(pars$Xpar[[i]],
       return(1/r1)
  )
}
