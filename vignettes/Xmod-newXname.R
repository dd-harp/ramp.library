#' @title Compute the derivatives for parasite infection dynamics in human population strata
#' @description Implements [dXdt] for the newXname model
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.newXname <- function(t, y, pars, i) {

  # do not change this
  foi <- pars$FoI[[i]]

  # attach the variables by name
  with(list_Xvars(y, pars, i),{
    # compute H (if it isn't one of the variables)
    H <- F_H(t, y, pars, i)

    # expose the parameters (see make_Xpar_Xname)
    with(pars$Xpar[[i]], {
      # compute the derivatives
      dX1 <- ...
      dX2 <- ...
      ...

      # concatenate the derivatives
      derivs = c(dX1, dX2, ...)

      # return the derivatives
      return(derivs)
    })
  })
}



#' @title Make initial values for the Xname human model, with defaults
#' @param nStrata the number of strata in the model
#' @param Xopts a [list] to overwrite defaults
#' @param X10 the initial value for X1
#' @param X20 the initial value for X1
#' @return a [list]
#' @export
make_Xinits_newXname = function(nStrata, Xopts = list(), X10=NULL, X20=1){with(Xopts,{
  stopifnot(is.numeric(X10))
  stopifnot(is.numeric(X20))
  X1 = checkIt(X10, nStrata)
  X2 = checkIt(X20, nStrata)
  return(list(X1=X1, X2=X2, ...))
})}





#' @title Setup Xinits.newXname
#' @description Implements [setup_Xinits] for the newXname model
#' @inheritParams ramp.xde::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.newXname = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, make_Xinits_newXname(pars$Hpar[[i]]$nStrata, Xopts, H0=Hpar[[i]]$H))
  return(pars)
}





#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the newXname model.
#' @inheritParams make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.newXname <- function(pars, i) {with(pars,{

  X1_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(X1_ix, 1)

  X2_ix <- seq(from = max_ix+1, length.out=Hpar[[i]]$nStrata)
  max_ix <- tail(X2_ix, 1)

  ...

  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(X1_ix=X1_ix, X2_ix=X2_ix, ...)
  return(pars)
})}





#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xde::list_Xvars
#' @return a [list]
#' @export
list_Xvars.newXname <- function(y, pars, i) {
  with(pars$ix$X[[i]],{
      X1 = y[X1_ix]
      X2 = y[X2_ix]
      ...
      H = X1+X2+...
      return(list(X1=X1,X2=X2,...,H=H))})
}





#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xde::get_inits_X
#' @return a [numeric] vector
#' @export
get_inits_X.newXname <- function(pars, i){
  with(pars$Xinits[[i]], return(c(X1,X2)))
}





#' @title Update inits for the newXname human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X
#' @return none
#' @export
update_inits_X.newXname <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = make_Xinits_newXname(pars, list(), X10=X1, X20=X2, ...)
    return(pars)
})}



#' @title Make parameters for newXname human model, with defaults
#' @param nStrata is the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param p1 the first parameter
#' @param p2 the second parameter
#' @param p3 the third parameter
#' @return a [list]
#' @export
make_Xpar_newXname = function(nStrata, Xopts=list(),
                         p1=1, p2=2, p3=3){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("newXname")

    Xpar$p1 = checkIt(p1, nStrata)
    Xpar$p2 = checkIt(p2, nStrata)
    Xpar$p3 = checkIt(p3, nStrata)

    return(Xpar)
  })}



#' @title Setup Xpar.newXname
#' @description Implements [setup_Xpar] for the newXname model
#' @inheritParams ramp.xde::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.newXname = function(newXname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_newXname(pars$Hpar[[i]]$nStrata, Xopts)
  return(pars)
}





#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.newXname <- function(y, pars, i) {
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
#' @description Implements [F_H] for the newXname model.
#' @inheritParams ramp.xde::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.newXname <- function(y, pars, i){
  with(list_Xvars(y, pars, i), return(H))
}





#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the newXname model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.newXname <- function(y, pars, i) {
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





#' @title Parse the output of deSolve and return variables for the newXname model
#' @description Implements [parse_deout_X] for the newXname model
#' @inheritParams ramp.xde::parse_deout_X
#' @return none
#' @export
parse_deout_X.newXname <- function(deout, pars, i) {
  time = deout[,1]
  with(pars$ix$X[[i]],{
    X1 = deout[,X1_ix+1]
    X2 = deout[,X2_ix+1]
    ...
    H = X1 + X2 + ...
    return(list(time=time, X1=X1, X2=X2, ..., H=H))
})}





#' @title Compute the true prevalence of infection / parasite rate
#' @description Implements F_pr for the newXname model.
#' @inheritParams F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.newXname <- function(vars, Xpar) {
  pr = with(vars, with(Xpar,(return(X1/H))))
  return(pr)
}





#' @title Compute the HTC for the newXname model
#' @description Implements [HTC] for the newXname model with demography.
#' @inheritParams ramp.xde::HTC
#' @return a [numeric] vector
#' @export
HTC.newXname <- function(pars, i) {
  with(pars$Xpar[[i]],
    #HTC <-
    return(HTC)
  )
}



#' Add lines for the density of infected individuals for the newXname model
#'
#' @param XH a list with the outputs of parse_deout_X_newXname
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
xds_lines_X_newXname = function(XH, nStrata, clrs=c("darkblue","darkred"), llty=1){
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





#' Plot the density of infected individuals for the newXname model
#'
#' @inheritParams ramp.xde::xds_plot_X
#' @export
xds_plot_X.newXname = function(pars, i=1, clrs=c("darkblue","darkred"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "Time"))

  xds_lines_X_newXname(vars$XH[[i]], pars$Hpar[[i]]$nStrata, clrs, llty)
}
