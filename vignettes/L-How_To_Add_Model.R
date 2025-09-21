## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [dLdt] for the newLname model
#' @inheritParams ramp.xds::dLdt
#' @return a [numeric] vector
#' @export
dLdt.newLname <- function(t, y, xds_obj, i) {

  # do not change this
  eta <- xds_obj$terms$eta[[i]]

  # do not change this
  # attach the variables by name 
  with(get_L_vars(y, xds_obj, i),{
    # do not change this
    # attach the parameters by name 
    with(xds_obj$L_obj[[i]], {
      # Do not change this 
      dL <- F_L(eta, L) 
      ... 
      
      # Change this 
      derivs = c(dL, ...) 
      
      # return the derivatives 
      return(derivs)
    })
  })
}

## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [UpdateLt] for the newLname model
#' @inheritParams ramp.xds::UpdateLt
#' @return a [numeric] vector
#' @export
Update_Lt.newLname<- function(t, y, xds_obj, s) {
  
  # do not change this
  eta <- xds_obj$terms$eta[[i]]

   # do not change this
  # attach the variables by name 
  with(get_L_vars(y, xds_obj, i),{
    # do not change this
    # attach the parameters by name 
    with(xds_obj$L_obj[[i]], {
      # Change this
      Lt <- L + Lambda + F_L(L) 
      ... 
      
      # Change this 
      states = c(Lt, ...) 
      
      # return the derivatives 
      return(states)
    })
  })
}

## -----------------------------------------------------------------------------
#' @title Lake parameters for newLname human model, with defaults
#' @param nPatches is the number of population strata
#' @param options a [list] that could overwrite defaults
#' @param p1 the first parameter 
#' @param p2 the second parameter 
#' @param p3 the third parameter 
#' @return a [list]
#' @export
make_L_obj_newLname = function(nPatches, options=list(),
                         p1=1, p2=2, p3=3){
  with(options,{
    L_obj = list()
    class(L_obj) <- c("newLname")

    # Change this
    L_obj$p1 = checkIt(p1, nPatches)
    L_obj$p2 = checkIt(p2, nPatches)
    L_obj$p3 = checkIt(p3, nPatches)

    # Don't change this
     L_obj <- setup_K_obj(nPatches, L_obj) 
    
    return(L_obj)
  })}

## -----------------------------------------------------------------------------
#' @title Setup L_obj.newLname
#' @description Implements [setup_L_obj] for the newLname model
#' @inheritParams ramp.xds::setup_L_obj
#' @return a [list] vector
#' @export
setup_L_obj.newLname = function(newLname, xds_obj, i, options=list()){
  xds_obj$L_obj[[i]] = make_L_obj_newLname(xds_obj$nPatches, options)
  return(xds_obj)
}

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [setup_L_ix] for the newLname model.
#' @inheritParams ramp.xds::setup_L_ix
#' @return none
#' @importFrom utils tail
#' @export
setup_L_ix.newLname <- function(xds_obj, i) {with(xds_obj,{
  
  L_ix <- seq(from = max_ix+1, length.out=xds_obj$nPatches[i])
  max_ix <- tail(L_ix, 1)

  Y_ix <- seq(from = max_ix+1, length.out=xds_obj$nPatches[i])
  max_ix <- tail(Y_ix, 1)

  ... 
  
  xds_obj$max_ix = max_ix
  xds_obj$L_obj[[i]]$ix = list(L_ix=L_ix, Y_ix=Y_ix, ...)
  return(xds_obj)
})}

## -----------------------------------------------------------------------------
#' @title Return the variables as a list
#' @description This method dispatches on the type of `xds_obj$L_obj`
#' @inheritParams ramp.xds::get_L_vars
#' @return a [list]
#' @export
get_L_vars.newLname <- function(y, xds_obj, i) {
  with(xds_obj$L_obj[[i]]$ix,{
      H = y[H_ix]
      L = y[L_ix]
      Y = y[Y_ix]
      ...
      X0 = H-L-Y 
      return(list(H=H, L=L, Y=Y,...,X0=X0))})
}

## -----------------------------------------------------------------------------
#' @title parse the output of deSolve and return variables for the newLname model
#' @description Implements [parse_L_orbits] for the newLname model
#' @inheritParams ramp.xds::parse_L_orbits
#' @return a named list  
#' @export
parse_L_orbits.newLname <- function(outputs, xds_obj, i) {
  with(xds_obj$L_obj[[i]]$ix,{
    H <- outputs[,H_ix]
    L <- outputs[,L_ix]
    Y <- outputs[,Y_ix]
    X0 <- H-L-Y
    vars <- list(H=H, L=L, Y=Y, X0=X0)
    return(vars)
})}

## -----------------------------------------------------------------------------
#' @title Setup initial values for *newLname*
#' 
#' @inheritParams ramp.xds::setup_L_inits
#' 
#' @return a **`ramp.xds`** object 
#' 
#' @export
setup_L_inits.newLname = function(xds_obj, H, i, options=list()){
  xds_obj$L_obj[[i]]$inits = make_L_inits_newLname(xds_obj$nPatches[i], H, options)
  return(xds_obj)
}

## -----------------------------------------------------------------------------
#' @title Lake initial values for the Xname human model, with defaults
#' @param nPatches the number of strata in the model
#' @param options a [list] to overwrite defaults
#' @param L0 the initial value for L 
#' @param Y0 the initial value for L 
#' @return a [list]
#' @export
make_L_inits_newLname = function(nPatches, H, options = list(), L0=NULL, Y0=1){with(options,{
  stopifnot(is.numeric(L0))
  stopifnot(is.numeric(Y0))
  L = checkIt(L0, nPatches)
  Y = checkIt(Y0, nPatches)
  return(list(L=L, Y=Y, ...))
})}

## -----------------------------------------------------------------------------
#' @title Return the parameters as a list
#' 
#' @inheritParams ramp.xds::change_L_inits 
#' 
#' @return an **`xds`** object
#' @export
change_L_inits.SIS <- function(xds_obj, i=1, options=list()) {
  with(get_L_inits(xds_obj, i), with(options,{
    xds_obj = change_H(H, xds_obj, i)
    xds_obj$Xinits[[i]]$L = L 
    xds_obj$Xinits[[i]]$Y = Y 
    ... 
    return(xds_obj)
  }))}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_emerge] for the SIS model.
#' @inheritParams ramp.xds::F_emerge
#' @return a [numeric] vector of length `nPatches`
#' @export
F_emerge.newLname <- function(y, xds_obj, i) {
  with(get_L_vars(y, xds_obj, i), 
    with(xds_obj$L_obj[[i]], {
      emerge = ...  
      return(emerge)
}))}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_L] for the SIS model.
#' @inheritParams ramp.xds::F_L
#' @return a [numeric] vector of length `nPatches`
#' @export
F_L.newLname <- function(y, xds_obj, i) {
  with(get_L_vars(y, xds_obj, i), 
    with(xds_obj$L_obj[[i]], {
      L = ...  
      return(L)
}))}

## -----------------------------------------------------------------------------
#' @title Compute the true prevalence of infection / parasite rate
#' @description Implements F_pr for the newLname model.
#' @inheritParams F_pr
#' @return a [numeric] vector of length `nPatches`
#' @export
F_capacity.newLname <- function(vars, L_obj) {
  with(vars, 
    with(L_obj,
       return(...)
  ))
}

## -----------------------------------------------------------------------------
#' Add lines for the density of infected individuals for the newLname model
#'
#' @param L a list with the outputs of parse_deout_X_newLname
#' @param nPatches the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
xds_lines_X_newLname = function(L, nPatches, clrs=c("darkblue","darkred"), llty=1){
  with(L,{
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
#' Plot the density of infected individuals for the newLname model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @export
xds_plot_X.newLname = function(xds_obj, i=1, clrs=c("darkblue","darkred"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(xds_obj$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$L[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "Time"))

  xds_lines_X_newLname(vars$L[[i]], xds_obj$Hpar[[i]]$nPatches, clrs, llty)
}

