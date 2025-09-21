## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [dMYdt] for the newMYname model
#' @inheritParams ramp.xds::dMYdt
#' @return a [numeric] vector
#' @export
dMYdt.newMYname <- function(t, y, xds_obj, i) {

  # do not change this
  Lambda <- xds_obj$terms$Lambda[[i]]
  kappa <- xds_obj$terms$kappa[[i]]

  # do not change this
  # attach the variables by name 
  with(get_MY_vars(y, xds_obj, i),{
    # do not change this
    # attach the parameters by name 
    with(xds_obj$MY_obj[[i]], {
      # Do not change this 
      dM <- F_M(Lambda, M) 
      dY <- F_Y(kappa, M, Y) 
      ... 
      
      # Change this 
      derivs = c(dM, dY, ...) 
      
      # return the derivatives 
      return(derivs)
    })
  })
}

## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [UpdateMYt] for the newMYname model
#' @inheritParams ramp.xds::UpdateMYt
#' @return a [numeric] vector
#' @export
Update_MYt.newMYname<- function(t, y, xds_obj, s) {
  Lambda <- xds_obj$terms$Lambda[[i]]
  kappa <- xds_obj$terms$kappa[[i]]

   # do not change this
  # attach the variables by name 
  with(get_MY_vars(y, xds_obj, i),{
    # do not change this
    # attach the parameters by name 
    with(xds_obj$MY_obj[[i]], {
      # Change this
      Mt <- M + Lambda + F_M(M) 
      Yt <- Y + F_Y(kappa, M, Y) 
      ... 
      
      # Change this 
      states = c(Mt, Yt, ...) 
      
      # return the derivatives 
      return(states)
    })
  })
}

## -----------------------------------------------------------------------------
#' @title Make parameters for newMYname human model, with defaults
#' @param nPatches is the number of population strata
#' @param options a [list] that could overwrite defaults
#' @param p1 the first parameter 
#' @param p2 the second parameter 
#' @param p3 the third parameter 
#' @return a [list]
#' @export
make_MY_obj_newMYname = function(nPatches, options=list(),
                         p1=1, p2=2, p3=3){
  with(options,{
    MY_obj = list()
    class(MY_obj) <- c("newMYname")

    # Change this
    MY_obj$p1 = checkIt(p1, nPatches)
    MY_obj$p2 = checkIt(p2, nPatches)
    MY_obj$p3 = checkIt(p3, nPatches)

    # Don't change this
     MY_obj <- setup_K_obj(nPatches, MY_obj) 
    
    return(MY_obj)
  })}

## -----------------------------------------------------------------------------
#' @title Setup MY_obj.newMYname
#' @description Implements [setup_MY_obj] for the newMYname model
#' @inheritParams ramp.xds::setup_MY_obj
#' @return a [list] vector
#' @export
setup_MY_obj.newMYname = function(newMYname, xds_obj, i, options=list()){
  xds_obj$MY_obj[[i]] = make_MY_obj_newMYname(xds_obj$nPatches, options)
  return(xds_obj)
}

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [setup_MY_ix] for the newMYname model.
#' @inheritParams ramp.xds::setup_MY_ix
#' @return none
#' @importFrom utils tail
#' @export
setup_MY_ix.newMYname <- function(xds_obj, i) {with(xds_obj,{
  
  M_ix <- seq(from = max_ix+1, length.out=xds_obj$nPatches[i])
  max_ix <- tail(M_ix, 1)

  Y_ix <- seq(from = max_ix+1, length.out=xds_obj$nPatches[i])
  max_ix <- tail(Y_ix, 1)

  ... 
  
  xds_obj$max_ix = max_ix
  xds_obj$MY_obj[[i]]$ix = list(M_ix=M_ix, Y_ix=Y_ix, ...)
  return(xds_obj)
})}

## -----------------------------------------------------------------------------
#' @title Return the variables as a list
#' @description This method dispatches on the type of `xds_obj$MY_obj`
#' @inheritParams ramp.xds::get_MY_vars
#' @return a [list]
#' @export
get_MY_vars.newMYname <- function(y, xds_obj, i) {
  with(xds_obj$MY_obj[[i]]$ix,{
      H = y[H_ix]
      M = y[M_ix]
      Y = y[Y_ix]
      ...
      X0 = H-M-Y 
      return(list(H=H, M=M, Y=Y,...,X0=X0))})
}

## -----------------------------------------------------------------------------
#' @title parse the output of deSolve and return variables for the newMYname model
#' @description Implements [parse_MY_orbits] for the newMYname model
#' @inheritParams ramp.xds::parse_MY_orbits
#' @return a named list  
#' @export
parse_MY_orbits.newMYname <- function(outputs, xds_obj, i) {
  with(xds_obj$MY_obj[[i]]$ix,{
    H <- outputs[,H_ix]
    M <- outputs[,M_ix]
    Y <- outputs[,Y_ix]
    X0 <- H-M-Y
    vars <- list(H=H, M=M, Y=Y, X0=X0)
    return(vars)
})}

## -----------------------------------------------------------------------------
#' @title Setup initial values for *newMYname*
#' 
#' @inheritParams ramp.xds::setup_MY_inits
#' 
#' @return a **`ramp.xds`** object 
#' 
#' @export
setup_MY_inits.newMYname = function(xds_obj, H, i, options=list()){
  xds_obj$MY_obj[[i]]$inits = make_MY_inits_newMYname(xds_obj$nPatches[i], H, options)
  return(xds_obj)
}

## -----------------------------------------------------------------------------
#' @title Make initial values for the Xname human model, with defaults
#' @param nPatches the number of strata in the model
#' @param options a [list] to overwrite defaults
#' @param M0 the initial value for M 
#' @param Y0 the initial value for M 
#' @return a [list]
#' @export
make_MY_inits_newMYname = function(nPatches, H, options = list(), M0=NULL, Y0=1){with(options,{
  stopifnot(is.numeric(M0))
  stopifnot(is.numeric(Y0))
  M = checkIt(M0, nPatches)
  Y = checkIt(Y0, nPatches)
  return(list(M=M, Y=Y, ...))
})}

## -----------------------------------------------------------------------------
#' @title Return the parameters as a list
#' 
#' @inheritParams ramp.xds::change_MY_inits 
#' 
#' @return an **`xds`** object
#' @export
change_MY_inits.SIS <- function(xds_obj, i=1, options=list()) {
  with(get_MY_inits(xds_obj, i), with(options,{
    xds_obj = change_H(H, xds_obj, i)
    xds_obj$Xinits[[i]]$M = M 
    xds_obj$Xinits[[i]]$Y = Y 
    ... 
    return(xds_obj)
  }))}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_fqZ] for the SIS model.
#' @inheritParams ramp.xds::F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ.newMYname <- function(y, xds_obj, i) {
  with(get_MY_vars(y, xds_obj, i), 
    with(xds_obj$MY_obj[[i]], {
      fqZ = ...  
      return(fqZ)
}))}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_fqM] for the SIS model.
#' @inheritParams ramp.xds::F_fqM
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqM.newMYname <- function(y, xds_obj, i) {
  with(get_MY_vars(y, xds_obj, i), 
    with(xds_obj$MY_obj[[i]], {
      fqM = ...  
      return(fqM)
}))}

## -----------------------------------------------------------------------------
#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_eggs] for the newMYname model.
#' @inheritParams ramp.xds::F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.newMYname <- function(y, xds_obj, i) {
  with(xds_obj$MY_obj[[i]],{ 
    ########################
    # retrieve or compute it 
    ########################
    eggs = ... 
    ########################
    # return it 
    ########################
    return(eggs)
  })
}

## -----------------------------------------------------------------------------
#' Add lines for the density of infected individuals for the newMYname model
#'
#' @param MY a list with the outputs of parse_deout_X_newMYname
#' @param nPatches the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
xds_lines_X_newMYname = function(MY, nPatches, clrs=c("darkblue","darkred"), llty=1){
  with(MY,{
    if(nPatches==1) {
      lines(time, S, col=clrs[1], lty = llty[1])
      lines(time, I, col=clrs[2], lty = llty[1])
    }
    if(nPatches>1){
      if (length(clrs)==2) clrs=matrix(clrs, 2, nPatches)
      if (length(llty)==1) llty=rep(llty, nPatches)

      for(i in 1:nPatches){
        lines(time, S[,i], col=clrs[1,i], lty = llty[i])
        lines(time, I[,i], col=clrs[2,i], lty = llty[i])
      }
    }
  })}

## -----------------------------------------------------------------------------
#' Plot the density of infected individuals for the newMYname model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @export
xds_plot_X.newMYname = function(xds_obj, i=1, clrs=c("darkblue","darkred"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(xds_obj$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$MY[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "Time"))

  xds_lines_X_newMYname(vars$MY[[i]], xds_obj$Hpar[[i]]$nPatches, clrs, llty)
}

