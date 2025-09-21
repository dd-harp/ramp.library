## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [dXHdt] for the newXname model
#' @inheritParams ramp.xds::dXHdt
#' @return a [numeric] vector
#' @export
dXHdt.newXname <- function(t, y, xds_obj, i) {

  # do not change this
  foi <- xds_obj$terms$FoI[[i]]

  # do not change this
  # attach the variables by name 
  with(get_XH_vars(y, xds_obj, i),{
    # do not change this
    # attach the parameters by name 
    with(xds_obj$XH_obj[[i]], {
      # Do not change this 
      dH <- Births(t, H, births) + D_matrix %*% H
      # Change this
      dX1 <- F_1(foi, X) + D_matrix %*% X1 
      dX2 <- F_2(foi, X) + D_matrix %*% X2 
      ... 
      
      # Change this 
      derivs = c(dH, dX1, dX2, ...) 
      
      # return the derivatives 
      return(derivs)
    })
  })
}

## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [UpdateXHt] for the newXname model
#' @inheritParams ramp.xds::UpdateXHt
#' @return a [numeric] vector
#' @export
Update_MYt.newXname<- function(t, y, xds_obj, s) {
  ar = xds_obj$terms$AR[[s]]

   # do not change this
  # attach the variables by name 
  with(get_XH_vars(y, xds_obj, i),{
    # do not change this
    # attach the parameters by name 
    with(xds_obj$XH_obj[[i]], {
      # Do not change this 
      Ht <- H + Births(t, H, births) + D_matrix %*% H
      # Change this
      X1t <- X1 + F_1(foi, X) + D_matrix %*% X1 
      X2t <- X2 + F_2(foi, X) + D_matrix %*% X2 
      ... 
      
      # Change this 
      states = c(Ht, X1t, X2t, ...) 
      
      # return the derivatives 
      return(states)
    })
  })
}

## -----------------------------------------------------------------------------
#' @title Make parameters for newXname human model, with defaults
#' @param nStrata is the number of population strata
#' @param options a [list] that could overwrite defaults
#' @param p1 the first parameter 
#' @param p2 the second parameter 
#' @param p3 the third parameter 
#' @return a [list]
#' @export
make_XH_obj_newXname = function(nStrata, options=list(),
                         p1=1, p2=2, p3=3){
  with(options,{
    XH_obj = list()
    class(XH_obj) <- c("newXname")

    # Change this
    XH_obj$p1 = checkIt(p1, nStrata)
    XH_obj$p2 = checkIt(p2, nStrata)
    XH_obj$p3 = checkIt(p3, nStrata)

    # Don't change this
    # Ports for human / host demography 
    XH_obj$D_matrix = diag(0, nStrata) 
    births = "zero"
    class(births) = births
    XH_obj$births = births 
    
    # Maybe delete this 
    # Ports for mass treatment 
    XH_obj$mda = F_zero 
    XH_obj$msat = F_zero 
    
    return(XH_obj)
  })}

## -----------------------------------------------------------------------------
#' @title Setup XH_obj.newXname
#' @description Implements [setup_XH_obj] for the newXname model
#' @inheritParams ramp.xds::setup_XH_obj
#' @return a [list] vector
#' @export
setup_XH_obj.newXname = function(newXname, xds_obj, i, options=list()){
  xds_obj$XH_obj[[i]] = make_XH_obj_newXname(xds_obj$nStrata[i], options)
  return(xds_obj)
}

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [setup_XH_ix] for the newXname model.
#' @inheritParams ramp.xds::setup_XH_ix
#' @return none
#' @importFrom utils tail
#' @export
setup_XH_ix.newXname <- function(xds_obj, i) {with(xds_obj,{
  
  X1_ix <- seq(from = max_ix+1, length.out=xds_obj$nStrata[i])
  max_ix <- tail(X1_ix, 1)

  X2_ix <- seq(from = max_ix+1, length.out=xds_obj$nStrata[i])
  max_ix <- tail(X2_ix, 1)

  ... 
  
  xds_obj$max_ix = max_ix
  xds_obj$XH_obj[[i]]$ix = list(X1_ix=X1_ix, X2_ix=X2_ix, ...)
  return(xds_obj)
})}

## -----------------------------------------------------------------------------
#' @title Return the variables as a list
#' @description This method dispatches on the type of `xds_obj$XH_obj`
#' @inheritParams ramp.xds::get_XH_vars
#' @return a [list]
#' @export
get_XH_vars.newXname <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
      H = y[H_ix]
      X1 = y[X1_ix]
      X2 = y[X2_ix]
      ...
      X0 = H-X1-X2 
      return(list(H=H, X1=X1, X2=X2,...,X0=X0))})
}

## -----------------------------------------------------------------------------
#' @title parse the output of deSolve and return variables for the newXname model
#' @description Implements [parse_XH_orbits] for the newXname model
#' @inheritParams ramp.xds::parse_XH_orbits
#' @return a named list  
#' @export
parse_XH_orbits.newXname <- function(outputs, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
    H <- outputs[,H_ix]
    X1 <- outputs[,X1_ix]
    X2 <- outputs[,X2_ix]
    X0 <- H-X1-X2
    vars <- list(H=H, X1=X1, X2=X2, X0=X0)
    return(vars)
})}

## -----------------------------------------------------------------------------
#' @title Setup initial values for *newXname*
#' 
#' @inheritParams ramp.xds::setup_XH_inits
#' 
#' @return a **`ramp.xds`** object 
#' 
#' @export
setup_XH_inits.newXname = function(xds_obj, H, i, options=list()){
  xds_obj$XH_obj[[i]]$inits = make_XH_inits_newXname(xds_obj$nStrata[i], H, options)
  return(xds_obj)
}

## -----------------------------------------------------------------------------
#' @title Make initial values for the Xname human model, with defaults
#' @param nStrata the number of strata in the model
#' @param options a [list] to overwrite defaults
#' @param X10 the initial value for X1 
#' @param X20 the initial value for X1 
#' @return a [list]
#' @export
make_XH_inits_newXname = function(nStrata, H, options = list(), X10=NULL, X20=1){with(options,{
  stopifnot(is.numeric(X10))
  stopifnot(is.numeric(X20))
  X1 = checkIt(X10, nStrata)
  X2 = checkIt(X20, nStrata)
  return(list(X1=X1, X2=X2, ...))
})}

## -----------------------------------------------------------------------------
#' @title Return the parameters as a list
#' 
#' @inheritParams ramp.xds::change_XH_inits 
#' 
#' @return an **`xds`** object
#' @export
change_XH_inits.SIS <- function(xds_obj, i=1, options=list()) {
  with(get_XH_inits(xds_obj, i), with(options,{
    xds_obj = change_H(H, xds_obj, i)
    xds_obj$Xinits[[i]]$X1 = X1 
    xds_obj$Xinits[[i]]$X2 = X2 
    ... 
    return(xds_obj)
  }))}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.newXname <- function(y, xds_obj, i) {
  with(get_XH_vars(y, xds_obj, i), 
    with(xds_obj$XH_obj[[i]], {
      X = ...  
      return(X)
}))}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_H] for the newXname model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.newXname <- function(y, xds_obj, i){
  return(get_XH_vars(y, xds_obj, i)$H)
}

## -----------------------------------------------------------------------------
#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_infectivity] for the newXname model.
#' @inheritParams ramp.xds::F_infectivity
#' @return a [numeric] vector of length `nStrata`
#' @export
F_infectivity.newXname <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]],{ 
    ########################
    # retrieve or compute it 
    ########################
    b = xds_obj$XH_obj[[i]]$b
    ########################
    # return it 
    ########################
    return(b)
  })
}

## -----------------------------------------------------------------------------
#' @title Compute the true prevalence of infection / parasite rate
#' @description Implements F_pr for the newXname model.
#' @inheritParams F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_prevalence.newXname <- function(vars, XH_obj) {
  pr = with(vars, with(XH_obj,(return(X1/H))))
  return(pr)
}

## -----------------------------------------------------------------------------
#' @title Compute the HTC for the newXname model
#' @description Implements [HTC] for the newXname model with demography.
#' @inheritParams ramp.xds::HTC
#' @return a [numeric] vector
#' @export
HTC.newXname <- function(xds_obj, i) {
  with(xds_obj$XH_obj[[i]],
    #HTC <- 
    return(HTC)
  )
}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
#' Plot the density of infected individuals for the newXname model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @export
xds_plot_X.newXname = function(xds_obj, i=1, clrs=c("darkblue","darkred"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(xds_obj$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "Time"))

  xds_lines_X_newXname(vars$XH[[i]], xds_obj$Hpar[[i]]$nStrata, clrs, llty)
}

