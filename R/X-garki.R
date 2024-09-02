# specialized methods for the human garki model

#' @title Derivatives for human population
#' @description Implements [dXdt] for the garki model.
#' @inheritParams ramp.xds::dXdt
#' @return a [numeric] vector
#' @export
dXdt.garki = function(t, y, pars, i){

  foi <- pars$FoI[[i]]
  Hpar <- pars$Hpar[[i]]

  with(list_Xvars(y, pars, i),{

    with(pars$Xpar[[i]],{
      R1 = foi/(exp(foi/r1) - 1)
      R2 = foi/(exp(foi/r2) - 1)

      #dx1 = Births(t, x1, pars) - FoI*x1 + R1*y2 + dHdt(t, x1, pars)
      dx1 = - foi*x1 + R1*y2 + dHdt(t, x1, Hpar)
      dx2 = foi*x1 - nu*x2 + dHdt(t, x2, Hpar)
      dy1 = nu*x2 - alpha1*y1  + dHdt(t, y1, Hpar)
      dy2 = alpha1*y1 - R1*y2 - alpha2*y2 + dHdt(t, y2, Hpar)
      dy3 = alpha2*y2 + nu*x4 - R2*y3 + dHdt(t, y3, Hpar)
      dx3 = R2*y3 - foi*x3 + dHdt(t, x3, Hpar)
      dx4 = foi*x3 - nu*x4 + dHdt(t, x4, Hpar)
      #dH =  dHdt(t, H, pars)

      return(c(dx1, dx2, dy1, dy2, dy3, dx3, dx4, dH))
    })
  })
}


#' @title make Xpar.garki
#' @description Implements [make_Xpar] for the garki model
#' @inheritParams ramp.xds::make_Xpar
#' @return a [list] vector
#' @export
make_Xpar.garki = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = create_Xpar_garki(pars$nStrata, Xopts)
  return(pars)
}

#' @title Make parameters for garki human model
#' @param nStrata is the number of population strata
#' @param Xopts an [list]
#' @param b transmission probability (efficiency) from mosquito to human
#' @param r1 a [numeric] recovery rate for non-immunes
#' @param r2 a [numeric]recovery rate for partially immunes
#' @param nu a [numeric]incubation period
#' @param alpha1 a [numeric] rate of losing infectivity
#' @param alpha2 a [numeric] rate of acquiring immunity
#' @param q1 a [numeric] detection of y1
#' @param q2 a [numeric] detection of y2
#' @param q3 a [numeric] detection of y3
#' @param mu a [numeric] the death rate
#' @return a [list]
#' @export
create_Xpar_garki = function(nStrata, Xopts=list(), b=0.55,
                             r1=.0023, r2=.023, nu=1/15,
                             alpha1=.002, alpha2=.00019,
                             q1=.7, q2=0.5, q3=0.3, mu=1/65/365){
  with(Xopts,{
    xde <- 'ode'
    class(xde) <- 'ode'

    garki = list()
    class(garki) <- "garki_xde"
    garki$xde <- xde
    garki$b=checkIt(b, nStrata)
    garki$r1=checkIt(r1, nStrata)
    garki$r2=checkIt(r2, nStrata)
    garki$nu=checkIt(nu, nStrata)
    garki$alpha1=checkIt(alpha1, nStrata)
    garki$alpha2=checkIt(alpha2, nStrata)
    garki$mu=checkIt(mu, nStrata)
    garki$q1=checkIt(q1, nStrata)
    garki$q2=checkIt(q2, nStrata)
    garki$q3=checkIt(q3, nStrata)

    return(garki)
  })}


#' @title Size of effective infectious human population
#' @description Implements [F_X] for the garki model.
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.garki <- function(t, y, pars, i){
  y1 <- y[pars$Xpar[[i]]$y1_ix]
  X = with(pars$Xpar[[i]], y1)
  return(X)
}

#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SIS model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.garki <- function(t, y, pars, i){
  with(list_Xvars(y, pars,i),{
    H = x1+x2+x3+x4+y1+y2+y3
    return(H)
  })
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the garki model.
#' @inheritParams ramp.xds::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.garki <- function(y, pars, i) {
  with(pars$Xpar[[i]], b)
}


#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Compute the true pr for the garki model.
#' @inheritParams ramp.xds::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.garki <- function(vars, Xpar) {
  pr = with(Xpar, with(vars, (y1+y2+y3)/H))
  return(pr)
}


#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xds::list_Xvars
#' @return a [list]
#' @export
list_Xvars.garki <- function(y, pars, i) {
  with(pars$ix$X[[i]],{
    x1 <- y[x1_ix]
    x2 <- y[x2_ix]
    x3 <- y[x3_ix]
    x4 <- y[x4_ix]
    y1 <- y[y1_ix]
    y2 <- y[y2_ix]
    y3 <- y[y3_ix]
    H = x1+x2+x3+x4+y1+y2+y3
    return(list(x1=x1,x2=x2,x3=x3,x4=x4,y1=y1,y2=y2,y3=y3))})
}


#' @title Make inits for garki human model. Note that the variables should sum up to H, so the initial value of x1 is not set. The values are passed in the same order as they are presented in the original paper.
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
create_Xinits_garki <- function(nStrata, Xopts = list(), H0=NULL, x1=NULL, x2=0, y1=0, y2=0, y3=0, x3=0, x4=0) {with(Xopts,{
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
})}

#' @title Setup Xinits.garki
#' @description Implements [make_Xinits] for the garki model
#' @inheritParams ramp.xds::make_Xinits
#' @return a [list] vector
#' @export
make_Xinits.garki = function(pars, H, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars,create_Xinits_garki(pars$nStrata[i], Xopts, H=H))
  return(pars)
}

#' @title Return initial values as a vector for the garki model
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xds::get_Xinits
#' @return a named [list]
#' @export
get_Xinits.garki <- function(pars, i){
  pars$Xinits[[i]]
}

#' @title Add indices for human population to parameter list
#' @description Implements [make_X_indices] for the garki model.
#' @inheritParams ramp.xds::make_X_indices
#' @return none
#' @importFrom utils tail
#' @export
make_X_indices.garki <- function(pars, i) {with(pars,{

  x1_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(x1_ix, 1)

  x2_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(x2_ix, 1)

  y1_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(y1_ix, 1)

  y2_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(y2_ix, 1)

  y3_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(y3_ix, 1)

  x3_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(x3_ix, 1)

  x4_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(x4_ix, 1)

  pars$max_ix <- max_ix
  pars$ix$X[[i]] = list(x1_ix = x1_ix,  x2_ix = x2_ix,
                        y1_iy = y1_iy,  y2_iy = y2_iy,  y3_iy = y3_iy,
                        x3_ix = x3_ix,  x4_ix = x4_ix)

  return(pars)
})}

#' @title Update inits for the garki model
#' @inheritParams ramp.xds::update_Xinits
#' @return an 'xds' object
#' @export
update_Xinits.garki <- function(pars, y0, i){
  with(list_Xvars(y0, pars, i),{
    pars = create_Xinits_garki(pars, x1=x1, x2=x2, y1=y1, y2=y2, y3=y3, x3=x3, x4=x4)
    return(pars)
  })}


#' @title Parse the output of deSolve and return variables for the garki model
#' @description Implements [parse_Xorbits] for the garki model
#' @inheritParams ramp.xds::parse_Xorbits
#' @return none
#' @export
parse_Xorbits.garki <- function(outputs, pars, i) {
  time = outputs[,1]
  with(pars$ix$X[[i]],{
    x1 = outputs[,x1_ix+1]
    x2 = outputs[,x2_ix+1]
    y1 = outputs[,y1_ix+1]
    y2 = outputs[,y2_ix+1]
    y3 = outputs[,y3_ix+1]
    x3 = outputs[,x3_ix+1]
    x4 = outputs[,x4_ix+1]
    H =x1+x2+y1+y2+y3+x3+x4
    vars <- list(time=time, x1=x1, x2=x2, y1=y1, y2=y2, y3=y3, x3=x3, x4=x4, H=H)
    return(vars)
  })}


#' Plot the density of infected individuals for the garki model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @export
xds_plot_X.garki = function(pars, i, clrs=viridisLite::turbo(7), llty=1, add=FALSE){
  XH = pars$outputs$orbits$XH[[i]]
  times = pars$outputs$time

  if(add==FALSE)
    with(XH,
         plot(times, 0*times, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "times"))

  xds_lines_X_garki(times, XH, pars, clrs, llty)
}

#' Add lines for the density of infected individuals for the garki model
#'
#' @param times time points for the observations
#' @param XH a list with the outputs of parse_outputs_X_SIS
#' @param pars a list that defines an `ramp.xds` model (*e.g.*,  generated by `make()`)
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
xds_lines_X_garki= function(times, XH, pars, clrs=viridisLite::turbo(7), llty=1){
  with(XH,{
    if(pars$nStrata==1){
      lines(times, x1, col=clrs[1], lty = llty[1])
      lines(times, x2, col=clrs[2], lty = llty[1])
      lines(times, y1, col=clrs[3], lty = llty[1])
      lines(times, y2, col=clrs[4], lty = llty[1])
      lines(times, y3, col=clrs[5], lty = llty[1])
      lines(times, x3, col=clrs[6], lty = llty[1])
      lines(times, x4, col=clrs[7], lty = llty[1])
    }
    if(pars$nStrata>1){
      if (length(clrs)==1) clrs=rep(clrs, pars$nStrata)
      if (length(llty)==1) llty=rep(llty, pars$nStrata)
      for(i in 1:pars$nStrata){
        lines(times, x1[,i], col=clrs[1], lty = llty[i])
        lines(times, x2[,i], col=clrs[2], lty = llty[i])
        lines(times, y1[,i], col=clrs[3], lty = llty[i])
        lines(times, y2[,i], col=clrs[4], lty = llty[i])
        lines(times, y3[,i], col=clrs[5], lty = llty[i])
        lines(times, x3[,i], col=clrs[6], lty = llty[i])
        lines(times, x4[,i], col=clrs[7], lty = llty[i])
      }
    }
  })}


#' @title Parse the output of deSolve and return variables for the garki model
#' @description Implements [parse_Xorbits] for the garki model
#' @inheritParams ramp.xds::parse_Xorbits
#' @return none
#' @export
parse_Xorbits.garki <- function(outputs, pars, i) {with(pars$ix$X[[i]],{
  x1 = outputs[,x1_ix]
  x2 = outputs[,x2_ix]
  y1 = outputs[,y1_ix]
  y2 = outputs[,y2_ix]
  y3 = outputs[,y3_ix]
  x3 = outputs[,x3_ix]
  x4 = outputs[,x4_ix]
  return(list(x1=x1, x2=x2, y1=y1, y2=y2, y3=y3, x3=x3, x4=x4, H=H))
})}


#' @title Compute the HTC for the garki model
#' @description Implements [HTC] for the garki model
#' @inheritParams ramp.xds::HTC
#' @return a [numeric] vector
#' @export
HTC.garki <- function(pars, i) {
  with(pars$Xpar[[i]],
       return(1/r1)
  )
}
