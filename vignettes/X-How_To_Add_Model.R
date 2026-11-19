## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [dXHdt] for the newXHname model
#' @inheritParams ramp.xds::dXHdt
#' @return a [numeric] vector
#' @export
dXHdt.newXHname <- function(t, y, xds_obj, i) {

  # do not change this
  foi <- xds_obj$terms$FoI[[i]]

  # do not change this
  # attach the variables by name 
  with(get_XH_vars(y, xds_obj, i),{
    # do not change this
    # attach the parameters by name 
    with(xds_obj$XH_obj[[i]], {
      # Do not change this 
      dH <- Births(t, y, xds_obj, i) - Deaths(t, H, xds_obj, i)
      # Change this
      dX1 <- F_1(foi, X1, X2) - Deaths(t, X1, xds_obj, i)
      dX2 <- F_2(foi, X1, X2) - Deaths(t, X2, xds_obj, i)
      
      # Change this 
      derivs <- c(dH, dX1, dX2)
      
      # return the derivatives 
      return(derivs)
    })
  })
}

## -----------------------------------------------------------------------------
#' @title Compute the derivatives for parasite infection dynamics in human population strata 
#' @description Implements [Update_XHt] for the newXHname model
#' @inheritParams ramp.xds::Update_XHt
#' @return a [numeric] vector
#' @export
Update_XHt.newXHname <- function(t, y, xds_obj, i) {
  ar <- xds_obj$terms$AR[[i]]

  # do not change this
  # attach the variables by name 
  with(get_XH_vars(y, xds_obj, i),{
    # do not change this
    # attach the parameters by name 
    with(xds_obj$XH_obj[[i]], {
      # For a static-population example; add demographic updates if needed.
      Ht <- H
      # Change this. Use ar for infection events in a discrete-time model.
      X1t <- X1 + F_1(ar, X1, X2)
      X2t <- X2 + F_2(ar, X1, X2)
      
      # Change this 
      states <- c(Ht, X1t, X2t)
      
      # return the updated state variables
      return(states)
    })
  })
}

## -----------------------------------------------------------------------------
#' @title Make a `newXHname` **XH** model object
#' @param nStrata is the number of population strata
#' @param options a [list] that could override defaults
#' @param p1 the first parameter 
#' @param p2 the second parameter 
#' @param p3 the third parameter 
#' @return a [list]
#' @export
make_XH_obj_newXHname = function(nStrata, options=list(),
                         p1=1, p2=2, p3=3){
  with(options,{
    XH_obj = list()
    class(XH_obj) <- c("newXHname")

    # Change this
    XH_obj$p1 = checkIt(p1, nStrata)
    XH_obj$p2 = checkIt(p2, nStrata)
    XH_obj$p3 = checkIt(p3, nStrata)

    return(XH_obj)
  })}

## -----------------------------------------------------------------------------
#' @title Setup XH_obj.newXHname
#' @description Implements [setup_XH_obj] for the newXHname model
#' @inheritParams ramp.xds::setup_XH_obj
#' @return an **`xds`** object
#' @export
setup_XH_obj.newXHname = function(XHname, residence, xds_obj, options=list(), i=1){
  xds_obj$XHname = "newXHname"
  xds_obj$XH_obj[[i]] = make_XH_obj_newXHname(xds_obj$nStrata[i], options)
  xds_obj <- setup_XH_inits(xds_obj, options=options, i=i)
  xds_obj$XH_obj[[i]]$skill_set <- get_skillset_XH(xds_obj, i)
  xds_obj <- setup_blood_search_weights("setup", xds_obj, i=i)
  xds_obj <- setup_timespent("setup", xds_obj, list(residence=residence), i)
  xds_obj <- setup_mass_treatment(xds_obj, i=i)
  xds_obj <- setup_demog_obj(xds_obj, i)
  xds_obj <- setup_time_away("setup", xds_obj, i=i)
  xds_obj <- setup_travel_eir("no_travel", xds_obj, i=i)
  return(xds_obj)
}

## -----------------------------------------------------------------------------
#' @title Get parameters for `newXHname` (**XH**)
#' @inheritParams ramp.xds::get_XH_pars
#' @return a [list]
#' @export
get_XH_pars.newXHname <- function(xds_obj, i=1) {
  with(xds_obj$XH_obj[[i]], list(p1=p1, p2=p2, p3=p3))
}

#' @title Change parameters for `newXHname` (**XH**)
#' @inheritParams ramp.xds::change_XH_pars
#' @return an **`xds`** object
#' @export
change_XH_pars.newXHname <- function(xds_obj, options=list(), i=1) {
  with(xds_obj$XH_obj[[i]], with(options, {
    if (!is.null(p1)) xds_obj$XH_obj[[i]]$p1 <- p1
    if (!is.null(p2)) xds_obj$XH_obj[[i]]$p2 <- p2
    if (!is.null(p3)) xds_obj$XH_obj[[i]]$p3 <- p3
    xds_obj
  }))
}

## -----------------------------------------------------------------------------
#' @title Add indices for human population variables
#' @description Implements [setup_XH_ix] for the newXHname model.
#' @inheritParams ramp.xds::setup_XH_ix
#' @return none
#' @importFrom utils tail
#' @export
setup_XH_ix.newXHname <- function(xds_obj, i) {with(xds_obj,{

  H_ix <- seq(from = max_ix+1, length.out=xds_obj$nStrata[i])
  max_ix <- tail(H_ix, 1)
  
  X1_ix <- seq(from = max_ix+1, length.out=xds_obj$nStrata[i])
  max_ix <- tail(X1_ix, 1)

  X2_ix <- seq(from = max_ix+1, length.out=xds_obj$nStrata[i])
  max_ix <- tail(X2_ix, 1)
  
  xds_obj$max_ix = max_ix
  xds_obj$XH_obj[[i]]$ix = list(H_ix=H_ix, X1_ix=X1_ix, X2_ix=X2_ix)
  return(xds_obj)
})}

## -----------------------------------------------------------------------------
#' @title Return the variables as a list
#' @description This method dispatches on the type of `xds_obj$XH_obj`
#' @inheritParams ramp.xds::get_XH_vars
#' @return a [list]
#' @export
get_XH_vars.newXHname <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
      H = y[H_ix]
      X1 = y[X1_ix]
      X2 = y[X2_ix]
      X0 = H-X1-X2 
      return(list(H=H, X1=X1, X2=X2, X0=X0))})
}

## -----------------------------------------------------------------------------
#' @title parse the output of deSolve and return variables for the newXHname model
#' @description Implements [parse_XH_orbits] for the newXHname model
#' @inheritParams ramp.xds::parse_XH_orbits
#' @return a named list  
#' @export
parse_XH_orbits.newXHname <- function(outputs, xds_obj, i) {
  with(xds_obj$XH_obj[[i]]$ix,{
    H <- outputs[,H_ix]
    X1 <- outputs[,X1_ix]
    X2 <- outputs[,X2_ix]
    X0 <- H-X1-X2
    vars <- list(H=H, X1=X1, X2=X2, X0=X0)
    return(vars)
})}

## -----------------------------------------------------------------------------
#' @title Setup initial values for *newXHname*
#' 
#' @inheritParams ramp.xds::setup_XH_inits
#' 
#' @return a **`ramp.xds`** object 
#' 
#' @export
setup_XH_inits.newXHname = function(xds_obj, options=list(), i=1){
  xds_obj$XH_obj[[i]]$inits = make_XH_inits_newXHname(xds_obj$nStrata[i], options=options)
  return(xds_obj)
}

## -----------------------------------------------------------------------------
#' @title Make initial values for the XHname human model, with defaults
#' @param nStrata the number of strata in the model
#' @param options a [list] to override defaults
#' @param X1 the initial value for X1 
#' @param X2 the initial value for X2
#' @return a [list]
#' @export
make_XH_inits_newXHname = function(nStrata, options = list(), H=1000, X1=1, X2=1){with(options,{
  H = unname(as.vector(checkIt(H, nStrata)))
  X1 = unname(as.vector(checkIt(X1, nStrata)))
  X2 = unname(as.vector(checkIt(X2, nStrata)))
  return(list(H=H, X1=X1, X2=X2))
})}

## -----------------------------------------------------------------------------
#' @title Change initial values for `newXHname` (**XH**)
#' 
#' @inheritParams ramp.xds::change_XH_inits 
#' 
#' @return an **`xds`** object
#' @export
change_XH_inits.newXHname <- function(xds_obj, options=list(), i=1) {
  with(get_XH_inits(xds_obj, i), with(options,{
    xds_obj$XH_obj[[i]]$inits$H = H
    xds_obj$XH_obj[[i]]$inits$X1 = X1
    xds_obj$XH_obj[[i]]$inits$X2 = X2
    return(xds_obj)
  }))}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_I] for the newXHname model.
#' @inheritParams ramp.xds::F_I
#' @return a [numeric] vector of length `nStrata`
#' @export
F_I.newXHname <- function(t, y, xds_obj, i) {
  with(get_XH_vars(y, xds_obj, i), {
    with(xds_obj$XH_obj[[i]], {
      infectious_density <- p1*X1
      return(infectious_density)
    })
  })
}

## -----------------------------------------------------------------------------
#' @title Size of human population
#' @description Implements [F_H] for the newXHname model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.newXHname <- function(t, y, xds_obj, i){
  return(get_XH_vars(y, xds_obj, i)$H)
}

## -----------------------------------------------------------------------------
#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_infectivity] for the newXHname model.
#' @inheritParams ramp.xds::F_infectivity
#' @return a [numeric] vector of length `nStrata`
#' @export
F_infectivity.newXHname <- function(y, xds_obj, i) {
  with(xds_obj$XH_obj[[i]],{ 
    ########################
    # retrieve or compute it 
    ########################
    b = xds_obj$XH_obj[[i]]$p1
    ########################
    # return it 
    ########################
    return(b)
  })
}

## -----------------------------------------------------------------------------
#' @title Compute the true prevalence of infection / parasite rate
#' @description Implements [F_prevalence] for the newXHname model.
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @export
F_prevalence.newXHname <- function(vars, XH_obj) {
  with(vars, X1/H)
}

## -----------------------------------------------------------------------------
#' @title Compute net infectiousness for `newXHname` (**XH**)
#' @description Implements [F_ni] for the newXHname model.
#' @inheritParams ramp.xds::F_ni
#' @return a [numeric] vector of length `nStrata`
#' @export
F_ni.newXHname <- function(vars, XH_obj) {
  with(vars, with(XH_obj, p1*X1/H))
}

## -----------------------------------------------------------------------------
#' @title Compute human transmitting capacity for `newXHname` (**XH**)
#' @description Implements [get_HTC] for the newXHname model.
#' @inheritParams ramp.xds::get_HTC
#' @return a [numeric] vector
#' @export
get_HTC.newXHname <- function(xds_obj, i) {
  stop("Replace this with the model's human transmitting capacity formula")
}

## -----------------------------------------------------------------------------
#' @title Get the skill set for `newXHname` (**XH**)
#' @inheritParams ramp.xds::get_skillset_XH
#' @return the **XH** module skill set, as a list
#' @export
get_skillset_XH.newXHname <- function(xds_obj, i) {
  list(
    delay=FALSE,
    demography=TRUE,
    prevalence=TRUE,
    malaria=TRUE,
    diagnostics="linear"
  )
}

#' @title Check `newXHname` before solving
#' @inheritParams ramp.xds::check_XH
#' @return an **`xds`** object
#' @export
check_XH.newXHname <- function(xds_obj, i) {
  xds_obj
}

## -----------------------------------------------------------------------------
#' @title Compute prevalence by light microscopy
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pfpr_by_lm.newXHname <- function(vars, XH_obj) {
  with(vars, with(XH_obj, q_lm*X1/H))
}

#' @title Compute prevalence by RDT
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pfpr_by_rdt.newXHname <- function(vars, XH_obj) {
  with(vars, with(XH_obj, q_rdt*X1/H))
}

#' @title Compute prevalence by PCR
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pfpr_by_pcr.newXHname <- function(vars, XH_obj) {
  with(vars, with(XH_obj, q_pcr*X1/H))
}

## -----------------------------------------------------------------------------
#' @title Compute steady states for `newXHname` (**XH**)
#' @inheritParams ramp.xds::steady_state_X
#' @return a named list of steady-state variables
#' @export
steady_state_X.newXHname <- function(foi, H, xds_obj, i=1) {
  stop("Replace this with the model's steady-state calculation")
}

## -----------------------------------------------------------------------------
#' Add lines for the `newXHname` state variables
#'
#' @param time simulation time points
#' @param XH a list with the outputs of parse_XH_orbits.newXHname
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
add_lines_X_newXHname = function(time, XH, nStrata, clrs=c("darkblue","darkred"), llty=1){
  with(XH,{
    if(nStrata==1) {
      lines(time, X1, col=clrs[1], lty = llty[1])
      lines(time, X2, col=clrs[2], lty = llty[1])
    }
    if(nStrata>1){
      if (length(clrs)==2) clrs=matrix(clrs, 2, nStrata)
      if (length(llty)==1) llty=rep(llty, nStrata)

      for(i in 1:nStrata){
        lines(time, X1[,i], col=clrs[1,i], lty = llty[i])
        lines(time, X2[,i], col=clrs[2,i], lty = llty[i])
      }
    }
  })}

## -----------------------------------------------------------------------------
#' Plot the density of infected individuals for the newXHname model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @export
xds_plot_X.newXHname = function(xds_obj, i=1, clrs=c("darkblue","darkred"), llty=1, add=FALSE){
  XH <- xds_obj$outputs$orbits$XH[[i]]
  time <- xds_obj$outputs$time
  if (add==FALSE)
    plot(time, 0*time, type="n", ylim=c(0, max(XH$H)),
         ylab="State density", xlab="Time")
  add_lines_X_newXHname(time, XH, xds_obj$nStrata[i], clrs, llty)
}

