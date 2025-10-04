
# specialized methods for the human garki model

#' @title The **XH** Module Skill Set
#'
#' @description The **XH** skill set is a list of
#' an module's capabilities.
#'
#' @note This method dispatches on `class(xds_obj$XH_obj)`
#'
#' @inheritParams ramp.xds::skill_set_XH
#'
#' @return the skill set, as a list
#'
#' @export
skill_set_XH.garki = function(Xname = "garki"){
  return(list(
    H_dynamics = TRUE,
    mda        = FALSE,
    msat       = FALSE,
    malaria    = TRUE,
    pr_obs     = TRUE,
    pf_lm      = TRUE,
    pf_rdt     = FALSE,
    pf_pcr     = FALSE
  ))
}

#' Check / update before solving
#'
#' @inheritParams ramp.xds::check_XH
#'
#' @returns an **`xds`** model object
#' @export
check_XH.garki = function(xds_obj, i){
  return(xds_obj)
}


#' @title Derivatives for human population
#' @importFrom Rdpack reprompt
#' @description Implements a continuous time version of the Garki model
#' @inheritParams ramp.xds::dXHdt
#' @return a [numeric] vector
#' @references{This implements a version of the model
#' developed for the Garki Project
#' \insertRef{DietzK1974GarkiModel}{ramp.library}}
#' @export
dXHdt.garki = function(t, y, xds_obj, i){

  foi <- xds_obj$terms$FoI[[i]]

  with(get_XH_vars(y, xds_obj, i),{
    with(xds_obj$XH_obj[[i]],{
      R1 = foi/(exp(foi/r1) - 1)
      R2 = foi/(exp(foi/r2) - 1)

      dH = Births(t, H, births) + D_matrix %*% H
#     dx1 = Births(t, H, Hpar) -foi*x1 + R1*y2 + dHdt(t, H, Hpar)
      dx2 = foi*x1 - nu*x2 + D_matrix %*% x2
      dy1 = nu*x2 - alpha1*y1  + D_matrix %*% y1
      dy2 = alpha1*y1 - R1*y2 - alpha2*y2 + D_matrix %*% y2
      dy3 = alpha2*y2 + nu*x4 - R2*y3 + D_matrix %*% y3
      dx3 = R2*y3 - foi*x3 + D_matrix %*% x3
      dx4 = foi*x3 - nu*x4 + D_matrix %*% x4
      return(c(dH, dx2, dy1, dy2, dy3, dx3, dx4))
    })
  })
}


#' @title make XH_obj for the Garki model
#' @description Implements [setup_XH_obj] for the garki model
#' @inheritParams ramp.xds::setup_XH_obj
#' @return a [list] vector
#' @export
setup_XH_obj.garki = function(Xname, xds_obj, i, options=list()){
  xds_obj$XH_obj[[i]] = make_XH_obj_garki(xds_obj$nStrata, options)
  return(xds_obj)
}

#' @title Make parameters for garki human model
#' @param nStrata is the number of population strata
#' @param options an [list]
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
make_XH_obj_garki = function(nStrata, options=list(), b=0.55,
                             r1=.0023, r2=.023, nu=1/15,
                             alpha1=.002, alpha2=.00019,
                             q1=.7, q2=0.5, q3=0.3, mu=1/65/365){
  with(options,{
    xde <- 'ode'
    class(xde) <- 'ode'

    garki = list()
    class(garki) <- "garki"
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

    garki$D_matrix = diag(0, nStrata)
    births = "zero"
    class(births) = births
    garki$births = births
    garki$mda = F_zero
    garki$msat = F_zero

    return(garki)
  })}


#' @title Return the parameters as a list
#' @description This method dispatches on the type of `xds_obj$XH_obj[[i]]`.
#' @inheritParams ramp.xds::change_XH_pars
#' @return an **`xds`** object
#' @export
change_XH_pars.garki <- function(xds_obj, i=1, options=list()) {
  nHabitats <- xds_obj$nHabitats
  with(xds_obj$XH_obj[[i]], with(options,{
    xds_obj$XH_obj[[i]]$b <- b
    xds_obj$XH_obj[[i]]$r1 <- r1
    xds_obj$XH_obj[[i]]$r2 <- r2
    xds_obj$XH_obj[[i]]$nu <- nu
    xds_obj$XH_obj[[i]]$alpha1 <- alpha1
    xds_obj$XH_obj[[i]]$alpha2 <- alpha2
    xds_obj$XH_obj[[i]]$mu <- mu
    xds_obj$XH_obj[[i]]$q1 <- q1
    xds_obj$XH_obj[[i]]$q2 <- q2
    xds_obj$XH_obj[[i]]$q3 <- q3
    return(xds_obj)
  }))}


#' @title Size of effective infectious human population
#' @description Implements [F_I] for the garki model.
#' @inheritParams ramp.xds::F_I
#' @return a [numeric] vector of length `nStrata`
#' @export
F_I.garki <- function(t, y, xds_obj, i){
  y1 <- y[xds_obj$XH_obj[[i]]$ix$y1_ix]
  return(y1)
}

#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SIS model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.garki <- function(t, y, xds_obj, i){
  with(get_XH_vars(y, xds_obj,i),{
    return(H)
  })
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_infectivity] for the garki model.
#' @inheritParams ramp.xds::F_infectivity
#' @return a [numeric] vector of length `nStrata`
#' @export
F_infectivity.garki <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]], b)
}


#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Compute the true pr for the garki model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @export
F_prevalence.garki <- function(vars, XH_obj) {
  pr = with(XH_obj, with(vars, (y1+y2+y3)/H))
  return(pr)
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_ni] for the SIR model.
#' @inheritParams ramp.xds::F_ni
#' @return a [numeric] vector of length `nStrata`
#' @export
F_ni.garki <- function(vars, XH_obj) {
  with(vars, with(XH_obj, (q1*y1+q2*y2+q3*y3)/H))
}

#' @title Return the variables as a list
#' @description This method dispatches on the type of `xds_obj$XH_obj`
#' @inheritParams ramp.xds::get_XH_vars
#' @return a [list]
#' @export
get_XH_vars.garki <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
    H <- y[H_ix]
    x2 <- y[x2_ix]
    x3 <- y[x3_ix]
    x4 <- y[x4_ix]
    y1 <- y[y1_ix]
    y2 <- y[y2_ix]
    y3 <- y[y3_ix]
    x1 <- H - x2 - x3 - x4 - y1 - y2 - y3
    return(list(H=H,x1=x1,x2=x2,x3=x3,x4=x4,y1=y1,y2=y2,y3=y3,H=H))})
}


#' @title Make inits for garki human model.
#' @note We use H instead of x, but other variables are passed in the same order as they are presented in the original paper.
#' @param nStrata is the number of population strata
#' @param options a [list] with values to override default values
#' @param H a [numeric] initial value for total human population density
#' @param x2 a [numeric] initial value for the variable x2
#' @param y1 a [numeric] initial value for the variable y1
#' @param y2 a [numeric] initial value for the variable y2
#' @param y3 a [numeric] initial value for the variable y3
#' @param x3 a [numeric] initial value for the variable x3
#' @param x4 a [numeric] initial value for the variable x4
#' @return none
#' @export
make_XH_inits_garki <- function(nStrata, H, options = list(), x2=0, y1=0, y2=0, y3=0, x3=0, x4=0) {
  stopifnot(is.numeric(x2))
  stopifnot(is.numeric(y1))
  stopifnot(is.numeric(y2))
  stopifnot(is.numeric(y3))
  stopifnot(is.numeric(x3))
  stopifnot(is.numeric(x4))
  stopifnot(H>0)

  H = checkIt(H, nStrata)
  x2 = checkIt(x2, nStrata)
  y1 = checkIt(y1, nStrata)
  y2 = checkIt(y2, nStrata)
  y3 = checkIt(y3, nStrata)
  x3 = checkIt(x3, nStrata)
  x4 = checkIt(x4, nStrata)
  return(list(H=H,x2=x2,y1=y1,y2=y2,y3=y3,x3=x3,x4=x4))
}

#' @title Set XH_inits.garki
#' @description Implements [change_XH_inits] for the garki model
#' @inheritParams ramp.xds::change_XH_inits
#' @return a [list] vector
#' @export
change_XH_inits.garki = function(xds_obj, i, options=list()){
  with(xds_obj$XH_obj[[i]]$inits,
    with(options,{
      xds_obj$XH_inits[[i]]$H = H
      xds_obj$XH_inits[[i]]$x2 = x2
      xds_obj$XH_inits[[i]]$y1 = y1
      xds_obj$XH_inits[[i]]$y2 = y2
      xds_obj$XH_inits[[i]]$y3 = y3
      xds_obj$XH_inits[[i]]$x3 = x3
      xds_obj$XH_inits[[i]]$x4 = x4
      return(xds_obj)
}))}

#' @title Setup XH_inits.garki
#' @description Implements [setup_XH_inits] for the garki model
#' @inheritParams ramp.xds::setup_XH_inits
#' @return a [list] vector
#' @export
setup_XH_inits.garki = function(xds_obj, H, i, options=list()){
  xds_obj$XH_obj[[i]]$inits = with(xds_obj, make_XH_inits_garki(xds_obj$nStrata[i], H, options))
  return(xds_obj)
}


#' @title Add indices for human population to parameter list
#' @description Implements [setup_XH_ix] for the garki model.
#' @inheritParams ramp.xds::setup_XH_ix
#' @return none
#' @importFrom utils tail
#' @export
setup_XH_ix.garki <- function(xds_obj, i) {with(xds_obj,{

  H_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(H_ix, 1)

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

  xds_obj$max_ix <- max_ix
  xds_obj$XH_obj[[i]]$ix = list(H_ix = H_ix,  x2_ix = x2_ix,
                        y1_ix = y1_ix,  y2_ix = y2_ix,  y3_ix = y3_ix,
                        x3_ix = x3_ix,  x4_ix = x4_ix)

  return(xds_obj)
})}


#' @title parse the output of deSolve and return variables for the garki model
#' @description Implements [parse_XH_orbits] for the garki model
#' @inheritParams ramp.xds::parse_XH_orbits
#' @return none
#' @export
parse_XH_orbits.garki <- function(outputs, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
    H = outputs[,H_ix]
    x2 = outputs[,x2_ix]
    y1 = outputs[,y1_ix]
    y2 = outputs[,y2_ix]
    y3 = outputs[,y3_ix]
    x3 = outputs[,x3_ix]
    x4 = outputs[,x4_ix]
    x1 <- H - x2 - x3 - x4 - y1 - y2 - y3
    return(list(H=H, x1=x1, x2=x2, y1=y1, y2=y2, y3=y3, x3=x3, x4=x4))
})}


#' Plot the density of infected individuals for the garki model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @export
xds_plot_X.garki = function(xds_obj, i=1, clrs=viridisLite::turbo(7), llty=1, add=FALSE){
  XH = xds_obj$outputs$orbits$XH[[i]]
  times = xds_obj$outputs$time

  if(add==FALSE)
    with(XH,
         plot(times, 0*times, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "times"))

  xds_lines_X_garki(times, XH, xds_obj, clrs, llty)
}

#' Add lines for the density of infected individuals for the garki model
#'
#' @param times time points for the observations
#' @param XH a list with the outputs of parse_outputs_X_SIS
#' @param xds_obj a list that defines an `ramp.xds` model (*e.g.*,  generated by `make()`)
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
xds_lines_X_garki= function(times, XH, xds_obj, clrs=viridisLite::turbo(8), llty=1){
  with(XH,{
    if(xds_obj$nStrata==1){
      lines(times, x1, col=clrs[1], lty = llty[1])
      lines(times, x2, col=clrs[2], lty = llty[1])
      lines(times, y1, col=clrs[3], lty = llty[1])
      lines(times, y2, col=clrs[4], lty = llty[1])
      lines(times, y3, col=clrs[5], lty = llty[1])
      lines(times, x3, col=clrs[6], lty = llty[1])
      lines(times, x4, col=clrs[7], lty = llty[1])
      lines(times, H, col=clrs[8], lty = llty[1])
    }
    if(xds_obj$nStrata>1){
      if (length(clrs)==1) clrs=rep(clrs, xds_obj$nStrata)
      if (length(llty)==1) llty=rep(llty, xds_obj$nStrata)
      for(i in 1:xds_obj$nStrata){
        lines(times, x1[,i], col=clrs[1], lty = llty[i])
        lines(times, x2[,i], col=clrs[2], lty = llty[i])
        lines(times, y1[,i], col=clrs[3], lty = llty[i])
        lines(times, y2[,i], col=clrs[4], lty = llty[i])
        lines(times, y3[,i], col=clrs[5], lty = llty[i])
        lines(times, x3[,i], col=clrs[6], lty = llty[i])
        lines(times, x4[,i], col=clrs[7], lty = llty[i])
        lines(times, H[,i], col=clrs[8], lty = llty[i])
      }
    }
  })}


#' @title Compute the HTC for the garki model
#' @description Implements [HTC] for the garki model
#' @inheritParams ramp.xds::HTC
#' @return a [numeric] vector
#' @export
HTC.garki <- function(xds_obj, i) {
  with(xds_obj$XH_obj[[i]],
       return(1/r1)
  )
}
