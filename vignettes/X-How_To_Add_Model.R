## -----------------------------------------------------------------------------
#' @title Compute the derivatives for `newXHname` (**XH**)
#' @description For the `newXHname` differential equation system:
#' \deqn{\frac{dH}{dt} = \mathrm{Births}(t,y) - \mathrm{Deaths}(t,H)}
#' \deqn{\frac{dX_1}{dt} = F_1(\mathrm{foi},X_1,X_2) - \mathrm{Deaths}(t,X_1)}
#' \deqn{\frac{dX_2}{dt} = F_2(\mathrm{foi},X_1,X_2) - \mathrm{Deaths}(t,X_2)}
#' Replace the placeholder functions and these equations together.
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
#' @title Update states for `newXHname` (**XH**)
#' @description For the `newXHname` discrete-time system, the updating function is:
#' \deqn{H_{t+1} = H_t}
#' \deqn{X_{1,t+1} = X_{1,t} + F_1(\mathrm{ar},X_{1,t},X_{2,t})}
#' \deqn{X_{2,t+1} = X_{2,t} + F_2(\mathrm{ar},X_{1,t},X_{2,t})}
#' Replace the placeholder functions and these equations together.
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
#' @title Get the state variables for `newXHname` (**XH**)
#' @description Extract the state variables for the `newXHname` module
#' from `y` and return them by name in a list for use during model solving.
#' Also compute the remaining class `X0 = H - X1 - X2`.
#' @note Variable indices are stored at `xds_obj$XH_obj[[i]]$ix`.
#' @inheritParams ramp.xds::get_XH_vars
#' @return a named [list] of state variables and module-specific derived quantities
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
#' @title Make a module object for `newXHname` (**XH**)
#' @description Create a `newXHname` module object and initialize its parameters.
#' Inspect parameters with [get_XH_pars.newXHname] and change supported values
#' with [change_XH_pars.newXHname]. Values in `options` override parameter arguments.
#' @seealso [get_XH_pars.newXHname], [change_XH_pars.newXHname]
#' @param nStrata is the number of population strata
#' @param options a [list] that could override defaults
#' @param p1 the first parameter 
#' @param p2 the second parameter 
#' @param p3 the third parameter 
#' @param q test positivity
#' @return a [list]
#' @export
make_XH_obj_newXHname = function(nStrata, options=list(),
                         p1=1, p2=2, p3=3, q=0.8){
  with(options,{
    XH_obj = list()
    class(XH_obj) <- c("newXHname")

    # Change this
    XH_obj$p1 = checkIt(p1, nStrata)
    XH_obj$p2 = checkIt(p2, nStrata)
    XH_obj$p3 = checkIt(p3, nStrata)
    XH_obj$q = checkIt(q, nStrata)

    return(XH_obj)
  })}

## -----------------------------------------------------------------------------
#' @title Set up the module object for `newXHname` (**XH**)
#' @description Set up a `newXHname` module object for the **XH** component
#' and the \eqn{i^{th}} host species.
#' + Set the module name to `newXHname`.
#' + Create the module object with [make_XH_obj_newXHname], passing `nStrata[i]` and `options`.
#' + Set initial values using [ramp.xds::setup_XH_inits], passing `options` to override defaults.
#' + Initialize blood search weights using [ramp.xds::setup_blood_search_weights.setup].
#' + Initialize the time-spent matrix using [ramp.xds::setup_timespent.setup].
#' + Initialize mass-treatment functions using [ramp.xds::setup_mass_treatment].
#' + Initialize birth and death rates using [ramp.xds::setup_demog_obj].
#' + Initialize time away using [ramp.xds::setup_time_away.setup].
#' + Initialize travel EIR using [ramp.xds::setup_travel_eir.no_travel].
#' + Return the updated **`xds`** model object.
#' @note The module object is stored at `xds_obj$XH_obj[[i]]`.
#' @inheritParams ramp.xds::setup_XH_obj
#' @return an **`xds`** model object
#' @export
setup_XH_obj.newXHname = function(XHname, residence, xds_obj, options=list(), i=1){
  xds_obj$XHname = "newXHname"
  xds_obj$XH_obj[[i]] = make_XH_obj_newXHname(xds_obj$nStrata[i], options)
  xds_obj <- setup_XH_inits(xds_obj, options=options, i=i)
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
#' @description Return the parameters for the `newXHname` module as a named list.
#' @seealso [make_XH_obj_newXHname], [change_XH_pars.newXHname]
#' @inheritParams ramp.xds::get_XH_pars
#' @return a named [list] of parameter values
#' @export
get_XH_pars.newXHname <- function(xds_obj, i=1) {
  with(xds_obj$XH_obj[[i]], list(p1=p1, p2=p2, p3=p3, q=q))
}

#' @title Change parameters for `newXHname` (**XH**)
#' @description Change parameter values named in `options` for the `newXHname`
#' module, leaving other values unchanged.
#' @seealso [make_XH_obj_newXHname], [get_XH_pars.newXHname]
#' @inheritParams ramp.xds::change_XH_pars
#' @return an **`xds`** model object
#' @export
change_XH_pars.newXHname <- function(xds_obj, options=list(), i=1) {
  with(xds_obj$XH_obj[[i]], with(options, {
    if (!is.null(p1)) xds_obj$XH_obj[[i]]$p1 <- p1
    if (!is.null(p2)) xds_obj$XH_obj[[i]]$p2 <- p2
    if (!is.null(p3)) xds_obj$XH_obj[[i]]$p3 <- p3
    xds_obj$XH_obj[[i]]$q <- q
    xds_obj
  }))
}

## -----------------------------------------------------------------------------
#' @title Make initial values for `newXHname` (**XH**)
#' @description Create a named list of initial values for the `newXHname` module.
#' Values in `options` override the defaults.
#' @details Initial values must be numeric scalars or have length `nStrata`.
#' Scalars are expanded to that length.
#' @param nStrata the number of strata in the model
#' @param options a [list] to override defaults
#' @param X1 the initial value for X1 
#' @param X2 the initial value for X2
#' @return a named [list] of initial values
#' @seealso [setup_XH_inits.newXHname], [change_XH_inits.newXHname]
#' @export
make_XH_inits_newXHname = function(nStrata, options = list(), H=1000, X1=1, X2=1){with(options,{
  H = unname(as.vector(checkIt(H, nStrata)))
  X1 = unname(as.vector(checkIt(X1, nStrata)))
  X2 = unname(as.vector(checkIt(X2, nStrata)))
  return(list(H=H, X1=X1, X2=X2))
})}

## -----------------------------------------------------------------------------
#' @title Set up initial values for `newXHname` (**XH**)
#' @description Set up initial values for the `newXHname` module using
#' [make_XH_inits_newXHname]. Values in `options` override the defaults.
#' @note Initial values are stored at `xds_obj$XH_obj[[i]]$inits`.
#' 
#' @inheritParams ramp.xds::setup_XH_inits
#' @seealso [make_XH_inits_newXHname], [change_XH_inits.newXHname]
#' 
#' @return an **`xds`** model object
#' 
#' @export
setup_XH_inits.newXHname = function(xds_obj, options=list(), i=1){
  xds_obj$XH_obj[[i]]$inits = make_XH_inits_newXHname(xds_obj$nStrata[i], options=options)
  return(xds_obj)
}

## -----------------------------------------------------------------------------
#' @title Change initial values for `newXHname` (**XH**)
#' @description Change initial values named in `options` for the `newXHname`
#' module, leaving other values unchanged.
#' @note Stored initial values for this module are at
#' `xds_obj$XH_obj[[i]]$inits`.
#' 
#' @inheritParams ramp.xds::change_XH_inits 
#' @seealso [make_XH_inits_newXHname], [setup_XH_inits.newXHname]
#' 
#' @return an **`xds`** model object
#' @export
change_XH_inits.newXHname <- function(xds_obj, options=list(), i=1) {
  with(get_XH_inits(xds_obj, i), with(options,{
    xds_obj$XH_obj[[i]]$inits$H = H
    xds_obj$XH_obj[[i]]$inits$X1 = X1
    xds_obj$XH_obj[[i]]$inits$X2 = X2
    return(xds_obj)
  }))}

## -----------------------------------------------------------------------------
#' @title Set up variable indices for `newXHname` (**XH**)
#' @description Assign consecutive variable indices for the `newXHname` module,
#' starting after `xds_obj$max_ix`, and update the maximum index.
#' Allocate blocks in the order `H_ix`, `X1_ix`, `X2_ix`.
#' Each block has length `nStrata[i]`.
#' @note Indices are stored at `xds_obj$XH_obj[[i]]$ix`.
#' The maximum index is stored at `xds_obj$max_ix`.
#' @inheritParams ramp.xds::setup_XH_ix
#' @return an **`xds`** model object
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
#' @title Parse outputs for `newXHname` (**XH**)
#' @description Extract the model state outputs for the `newXHname` module
#' from `outputs` and return them by name in a list.
#' Also compute `X0 = H-X1-X2`.
#' @note Variable indices are stored at `xds_obj$XH_obj[[i]]$ix`.
#' @inheritParams ramp.xds::parse_XH_orbits
#' @return a [list]
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
#' @title Compute observed prevalence
#' @description Return test positivity times true prevalence.
#' @inheritParams ramp.xds::F_pr_obs
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr_obs.newXHname <- function(vars, XH_obj) {
  XH_obj$q * ramp.xds::F_prevalence(vars, XH_obj)
}

## -----------------------------------------------------------------------------
#' @title Compute observed prevalence by light microscopy
#' @inheritParams ramp.xds::F_pfpr_by_lm
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pfpr_by_lm.newXHname <- function(vars, XH_obj) {
  with(vars, with(XH_obj, q_lm*X1/H))
}

#' @title Compute observed prevalence by RDT
#' @inheritParams ramp.xds::F_pfpr_by_rdt
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pfpr_by_rdt.newXHname <- function(vars, XH_obj) {
  with(vars, with(XH_obj, q_rdt*X1/H))
}

#' @title Compute observed prevalence by PCR
#' @inheritParams ramp.xds::F_pfpr_by_pcr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pfpr_by_pcr.newXHname <- function(vars, XH_obj) {
  with(vars, with(XH_obj, q_pcr*X1/H))
}

## -----------------------------------------------------------------------------
#' @title Check `newXHname` (**XH**)
#' @description Return the model unchanged; no checks are implemented in this template.
#' @inheritParams ramp.xds::check_XH
#' @return an **`xds`** model object
#' @export
check_XH.newXHname <- function(xds_obj, i) {
  xds_obj
}

#' @title Get the skill set for `newXHname` (**XH**)
#' @description Return the capabilities and compatibilities of the `newXHname` module as a list.
#' @inheritParams ramp.xds::get_skillset_XH
#' @return the **XH** module skill set, as a list
#' @export
get_skillset_XH.newXHname <- function(xds_obj, i) {
  list(
    dts=FALSE,
    xde=TRUE,
    delay=FALSE,
    demography=TRUE,
    malaria=TRUE,
    pr_obs=TRUE,
    pf_lm=TRUE,
    pf_rdt=TRUE,
    pf_pcr=TRUE,
    mda=FALSE,
    msat=FALSE
  )
}

## -----------------------------------------------------------------------------
#' @title Compute steady states for `newXHname` (**XH**)
#' @description Compute steady-state infection densities for the `newXHname`
#' module given daily FoI `foi` and fixed host population density `H`.
#' @note Document the equation system and assumptions supported by this method.
#' @inheritParams ramp.xds::steady_state_X
#' @return a named [list] of steady-state variables; replace this description
#' with the actual return type and state names
#' @export
steady_state_X.newXHname <- function(foi, H, xds_obj, i=1) {
  stop("Replace this with the model's steady-state calculation")
}

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

