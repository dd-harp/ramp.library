# TODO: Reconcile parameter coverage: make_L_obj_stages initializes w1, w2,
# and w3, and get_L_pars.stages returns them, but change_L_pars.stages does
# not update them. Decide whether to support changes or document the restriction.
# TODO: Reconcile the Emergence anchor formula psi_P*P with F_emerge.stages,
# which returns psi_P*P/2. Confirm the intended female-emergence convention.

#' @title The `stages` module for the **L** component
#' @description
#' A staged aquatic mosquito development model with four larval instars.
#' Population density is structured into eggs (\eqn{E}), four larval instars
#' (\eqn{L_1}–\eqn{L_4}), and pupae (\eqn{P}). Density-dependent regulation
#' operates through a weighted mean crowding index
#' \deqn{L = w_1 L_1 + w_2 L_2 + w_3 L_3 + L_4,}
#' allowing each instar to contribute differently to crowding. Each larval
#' instar has its own maturation rate, density-dependent maturation delay,
#' density-independent mortality, and density-dependent mortality.
#'
#' @section State Variables:
#' \describe{
#'   \item{`E`}{egg density per habitat}
#'   \item{`L1`}{first instar larval density per habitat}
#'   \item{`L2`}{second instar larval density per habitat}
#'   \item{`L3`}{third instar larval density per habitat}
#'   \item{`L4`}{fourth instar larval density per habitat}
#'   \item{`P`}{pupal density per habitat}
#' }
#'
#' @section Parameters:
#' \describe{
#'   \item{`psi_E`}{egg maturation (hatching) rate (\eqn{\psi_E})}
#'   \item{`phi_E`}{egg mortality rate (\eqn{\phi_E})}
#'   \item{`psi1`, `psi2`, `psi3`, `psi4`}{larval instar maturation rates (\eqn{\psi_i})}
#'   \item{`xi1`, `xi2`, `xi3`, `xi4`}{delayed maturation responses to mean crowding (\eqn{\xi_i})}
#'   \item{`phi1`, `phi2`, `phi3`, `phi4`}{density-independent larval mortality rates (\eqn{\phi_i})}
#'   \item{`theta1`, `theta2`, `theta3`, `theta4`}{density-dependent larval mortality slopes (\eqn{\theta_i})}
#'   \item{`phi_P`}{pupal mortality rate (\eqn{\phi_P})}
#'   \item{`psi_P`}{pupal maturation (emergence) rate (\eqn{\psi_P})}
#'   \item{`w1`, `w2`, `w3`}{contributions of \eqn{L_1}, \eqn{L_2}, and \eqn{L_3} to mean crowding (\eqn{w_i})}
#' }
#'
#' @section Egg Laying:
#' \describe{
#'   \item{`eta`}{the egg laying rate (\eqn{\eta})}
#' }
#'
#' @section Mean Crowding:
#'
#' \deqn{L = w_1 L_1 + w_2 L_2 + w_3 L_3 + L_4}
#'
#' The fourth instar \eqn{L_4} acts as the reference stage (\eqn{w_4 = 1}).
#' Setting all \eqn{w_i = 1} recovers unweighted total larval density.
#'
#' @section Dynamics:
#'
#' \deqn{dE/dt = \eta - (\psi_E + \phi_E) E}
#' \deqn{dL_1/dt = \psi_E E - \psi_1 L_1 e^{-\xi_1 L} - (\phi_1 + \theta_1 L) L_1}
#' \deqn{dL_2/dt = \psi_1 L_1 e^{-\xi_1 L} - \psi_2 L_2 e^{-\xi_2 L} - (\phi_2 + \theta_2 L) L_2}
#' \deqn{dL_3/dt = \psi_2 L_2 e^{-\xi_2 L} - \psi_3 L_3 e^{-\xi_3 L} - (\phi_3 + \theta_3 L) L_3}
#' \deqn{dL_4/dt = \psi_3 L_3 e^{-\xi_3 L} - \psi_4 L_4 e^{-\xi_4 L} - (\phi_4 + \theta_4 L) L_4}
#' \deqn{dP/dt = \psi_4 L_4 e^{-\xi_4 L} - (\phi_P + \psi_P) P}
#'
#' @section Emergence:
#'
#' The emergence rate of adult, female mosquitoes from each habitat is:
#' \deqn{\alpha = \psi_P P}
#'
#' @name stages
#' @rdname stages
NULL

#' @title Compute the derivatives for `stages` (**L**)
#'
#' @description Compute the derivatives for the `stages` module:
#'
#' \deqn{dE/dt = \eta - (\psi_E + \phi_E) E}
#' \deqn{dL_1/dt = \psi_E E - \psi_1 L_1 e^{-\xi_1 L} - (\phi_1 + \theta_1 L) L_1}
#' \deqn{dL_2/dt = \psi_1 L_1 e^{-\xi_1 L} - \psi_2 L_2 e^{-\xi_2 L} - (\phi_2 + \theta_2 L) L_2}
#' \deqn{dL_3/dt = \psi_2 L_2 e^{-\xi_2 L} - \psi_3 L_3 e^{-\xi_3 L} - (\phi_3 + \theta_3 L) L_3}
#' \deqn{dL_4/dt = \psi_3 L_3 e^{-\xi_3 L} - \psi_4 L_4 e^{-\xi_4 L} - (\phi_4 + \theta_4 L) L_4}
#' \deqn{dP/dt = \psi_4 L_4 e^{-\xi_4 L} - (\phi_P + \psi_P) P}
#'
#' @inheritSection stages State Variables
#' @inheritSection stages Egg Laying
#' @inheritSection stages Parameters
#' @inheritSection stages Mean Crowding
#' @inheritSection stages Dynamics
#' @inheritParams ramp.xds::dLdt
#' @note Derivatives are returned in blocks for `E`, `L1`, `L2`, `L3`,
#' `L4`, and `P`, each in habitat order.
#' The egg input is read from `xds_obj$terms$eta[[s]]`.
#' @return a [numeric] vector of derivatives in the module's state-variable
#' order, with length `6 * nHabitats`
#' @seealso [stages], [make_L_obj_stages], [ramp.xds::EggLaying]
#' @keywords internal
#' @export
dLdt.stages <- function(t, y, xds_obj, s) {
  eta <- as.vector(xds_obj$terms$eta[[s]])
  with(get_L_vars(y, xds_obj, s), {
    with(xds_obj$L_obj[[s]], {
      L  <- w1*L1 + w2*L2 + w3*L3 + L4
      m1 <- psi1 * L1 * exp(-xi1 * L)
      m2 <- psi2 * L2 * exp(-xi2 * L)
      m3 <- psi3 * L3 * exp(-xi3 * L)
      m4 <- psi4 * L4 * exp(-xi4 * L)
      dE  <- eta        - (psi_E + phi_E) * E
      dL1 <- psi_E * E  - m1 - (phi1 + theta1 * L) * L1
      dL2 <- m1         - m2 - (phi2 + theta2 * L) * L2
      dL3 <- m2         - m3 - (phi3 + theta3 * L) * L3
      dL4 <- m3         - m4 - (phi4 + theta4 * L) * L4
      dP  <- m4         - (phi_P + psi_P) * P
      return(c(dE, dL1, dL2, dL3, dL4, dP))
    })
  })
}

#' @title Get the state variables for `stages` (**L**)
#' @description Extract the state variables for the `stages` module
#' from `y` and return them by name in a list for use during model solving.
#' @inheritParams ramp.xds::get_L_vars
#' @return a named [list] of state variables and module-specific derived quantities
#' @keywords internal
#' @export
#' @note Variable indices are stored at `xds_obj$L_obj[[s]]$ix`.
get_L_vars.stages <- function(y, xds_obj, s) {
  with(xds_obj$L_obj[[s]]$ix, {
    list(
      E  = y[E_ix],
      L1 = y[L1_ix],
      L2 = y[L2_ix],
      L3 = y[L3_ix],
      L4 = y[L4_ix],
      P  = y[P_ix]
    )
  })
}

#' @title Mosquito bionomics for `stages` (**L**)
#'
#' @description Leave directly stored parameters unchanged. This module
#' does not configure bionomic parameter ports.
#' @inheritSection stages Parameters
#' @inheritParams ramp.xds::LBionomics
#' @seealso [stages], [LEffectSizes.stages]
#'
#' @return an **`xds`** model object
#' @keywords internal
#' @export
LBionomics.stages <- function(t, y, xds_obj, s) {
  return(xds_obj)
}


#' @title Apply effect sizes for `stages` (**L**)
#' @description Leave directly stored parameters unchanged. This module
#' does not configure bionomic parameter ports or apply effect sizes.
#' @inheritSection stages Parameters
#' @inheritParams ramp.xds::LEffectSizes
#' @seealso [stages], [LBionomics.stages]
#' @return an **`xds`** model object
#' @keywords internal
#' @export
LEffectSizes.stages <- function(t, y, xds_obj, s) {
  return(xds_obj)
}

#' @title Compute adult mosquito emergence for `stages` (**L**)
#' @description Return adult female mosquito emergence per habitat.
#' @inheritSection stages Emergence
#' @note The implementation returns `psi_P*P/2`, whereas the inherited
#' Emergence section gives `psi_P*P`; the factor of two requires review.
#' @inheritParams ramp.xds::F_emerge
#' @return a [numeric] vector of length `nHabitats`
#' @seealso [stages], [dLdt.stages], [ramp.xds::Emergence]
#' @keywords internal
#' @export
F_emerge.stages <- function(t, y, xds_obj, s) {
  P <- y[xds_obj$L_obj[[s]]$ix$P_ix]
  with(xds_obj$L_obj[[s]], {
    return(psi_P * P/2)
  })
}


#' @title Make a module object for `stages` (**L**)
#' @description Create a `stages` module object and initialize its parameters.
#' Inspect parameters with [get_L_pars.stages] and change supported values
#' with [change_L_pars.stages]. Values in `options` override parameter arguments.
#'
#' The following parameters will be set to the values in
#' `options`. If they are not found, default values will be used.
#'
#' - \eqn{\psi_E} or `psi_E`: egg maturation (hatching) rate
#' - \eqn{\phi_E} or `phi_E`: egg mortality rate
#' - \eqn{\psi_i} or `psi1`, `psi2`, `psi3`, `psi4`: larval instar maturation rates
#' - \eqn{\xi_i} or `xi1`, `xi2`, `xi3`, `xi4`: delayed maturation responses to crowding
#' - \eqn{\phi_i} or `phi1`, `phi2`, `phi3`, `phi4`: density-independent larval mortality rates
#' - \eqn{\theta_i} or `theta1`, `theta2`, `theta3`, `theta4`: density-dependent larval mortality slopes
#' - \eqn{\phi_P} or `phi_P`: pupal mortality rate
#' - \eqn{\psi_P} or `psi_P`: pupal maturation (emergence) rate
#' - \eqn{w_1, w_2, w_3} or `w1`, `w2`, `w3`: instar crowding weights
#'
#' @param nHabitats the number of habitats in the model
#' @param options a named [list]
#'
#' @param psi_E egg maturation (hatching) rate
#' @param phi_E egg mortality rate
#' @param psi1 first instar maturation rate
#' @param xi1 first instar delayed maturation response to crowding
#' @param phi1 first instar density-independent mortality rate
#' @param theta1 first instar density-dependent mortality slope
#' @param psi2 second instar maturation rate
#' @param xi2 second instar delayed maturation response to crowding
#' @param phi2 second instar density-independent mortality rate
#' @param theta2 second instar density-dependent mortality slope
#' @param psi3 third instar maturation rate
#' @param xi3 third instar delayed maturation response to crowding
#' @param phi3 third instar density-independent mortality rate
#' @param theta3 third instar density-dependent mortality slope
#' @param psi4 fourth instar maturation rate
#' @param xi4 fourth instar delayed maturation response to crowding
#' @param phi4 fourth instar density-independent mortality rate
#' @param theta4 fourth instar density-dependent mortality slope
#' @param phi_P pupal mortality rate
#' @param psi_P pupal maturation (emergence) rate
#' @param w1 contribution of \eqn{L_1} to mean crowding
#' @param w2 contribution of \eqn{L_2} to mean crowding
#' @param w3 contribution of \eqn{L_3} to mean crowding
#'
#' @seealso Called by: [setup_L_obj.stages]. Related: [dLdt.stages]
#' @seealso [get_L_pars.stages], [change_L_pars.stages]
#' @return **`L_obj`** an **L** component object
#' @keywords internal
#' @export
make_L_obj_stages = function(nHabitats, options = list(),
    psi_E  = 1/2,  phi_E  = 1/10,
    psi1   = 1/2,  xi1    = 0,    phi1   = 1/2, theta1 = 1/100,
    psi2   = 1/2,  xi2    = 0,    phi2   = 1/2, theta2 = 1/100,
    psi3   = 1/2,  xi3    = 0,    phi3   = 1/2, theta3 = 1/100,
    psi4   = 1/2,  xi4    = 0,    phi4   = 1/2, theta4 = 1/100,
    phi_P  = 1/2,  psi_P  = 1/3,
    w1 = .3, w2 = .5, w3 = .8) {
  with(options, {
    L_obj <- list()
    class(L_obj) <- "stages"
    L_obj$nHabitats <- nHabitats

    L_obj$psi_E  <- checkIt(psi_E,  nHabitats)
    L_obj$phi_E  <- checkIt(phi_E,  nHabitats)
    L_obj$psi1   <- checkIt(psi1,   nHabitats)
    L_obj$xi1    <- checkIt(xi1,    nHabitats)
    L_obj$phi1   <- checkIt(phi1,   nHabitats)
    L_obj$theta1 <- checkIt(theta1, nHabitats)
    L_obj$psi2   <- checkIt(psi2,   nHabitats)
    L_obj$xi2    <- checkIt(xi2,    nHabitats)
    L_obj$phi2   <- checkIt(phi2,   nHabitats)
    L_obj$theta2 <- checkIt(theta2, nHabitats)
    L_obj$psi3   <- checkIt(psi3,   nHabitats)
    L_obj$xi3    <- checkIt(xi3,    nHabitats)
    L_obj$phi3   <- checkIt(phi3,   nHabitats)
    L_obj$theta3 <- checkIt(theta3, nHabitats)
    L_obj$psi4   <- checkIt(psi4,   nHabitats)
    L_obj$xi4    <- checkIt(xi4,    nHabitats)
    L_obj$phi4   <- checkIt(phi4,   nHabitats)
    L_obj$theta4 <- checkIt(theta4, nHabitats)
    L_obj$phi_P  <- checkIt(phi_P,  nHabitats)
    L_obj$psi_P  <- checkIt(psi_P,  nHabitats)
    L_obj$w1     <- checkIt(w1,     nHabitats)
    L_obj$w2     <- checkIt(w2,     nHabitats)
    L_obj$w3     <- checkIt(w3,     nHabitats)

    return(L_obj)
  })
}

#' @title Set up the module object for `stages` (**L**)
#' @description
#'
#' Set up a `stages` module object for the **L** component and the \eqn{s^{th}} vector species.
#'
#' + Create the module object with [make_L_obj_stages], passing `nHabitats` and `options`.
#' + Set the skill set using [ramp.xds::get_skillset_L].
#' + Set initial values using [ramp.xds::setup_L_inits], passing `options` to override defaults.
#' + Initialize habitat search weights using [ramp.xds::setup_habitat_search_weights.setup].
#' + Return the updated **`xds`** model object.
#'
#' @seealso [stages], [make_L_obj_stages]
#'
#' @note The module object is stored at `xds_obj$L_obj[[s]]`
#' Its classes are `stages` and `stages_` followed by `xds_obj$xds`.
#'
#' @inheritParams ramp.xds::setup_L_obj
#' @return an **`xds`** model object
#' @keywords internal
#' @export
setup_L_obj.stages = function(Lname, membership, xds_obj, options=list(), s=1){
  L_obj <- make_L_obj_stages(xds_obj$nHabitats, options)
  class(L_obj) <- c("stages", paste("stages_", xds_obj$xds, sep=""))
  xds_obj$L_obj[[s]] = L_obj
  xds_obj$L_obj[[s]]$skill_set <- get_skillset_L(xds_obj, s)
  xds_obj <- setup_L_inits(xds_obj, options=options, s=s)
  xds_obj <- setup_habitat_search_weights("setup", xds_obj, list(membership=membership), s=s)
  return(xds_obj)
}

#' @title Get parameters for `stages` (**L**)
#' @seealso [make_L_obj_stages], [change_L_pars.stages]
#' @description Return the parameters for the `stages` module as a named list.
#' @inheritSection stages Parameters
#' @inheritParams ramp.xds::get_L_pars
#' @return a named [list] of parameter values
#' @seealso [dLdt.stages] or [change_L_pars.stages]
#' @keywords internal
#' @export
get_L_pars.stages <- function(xds_obj, s = 1) {
  with(xds_obj$L_obj[[s]], list(
    psi_E  = psi_E,  phi_E  = phi_E,
    psi1   = psi1,   xi1    = xi1,    phi1   = phi1,   theta1 = theta1,
    psi2   = psi2,   xi2    = xi2,    phi2   = phi2,   theta2 = theta2,
    psi3   = psi3,   xi3    = xi3,    phi3   = phi3,   theta3 = theta3,
    psi4   = psi4,   xi4    = xi4,    phi4   = phi4,   theta4 = theta4,
    phi_P  = phi_P,  psi_P  = psi_P,
    w1     = w1,     w2     = w2,     w3     = w3
  ))
}

#' @title Change parameters for `stages` (**L**)
#' @seealso [make_L_obj_stages], [get_L_pars.stages]
#' @description Change parameter values named in `options` for the `stages`
#' module, leaving other values unchanged.
#'
#' Named values passed in `options` update the corresponding
#' parameter fields directly; this module has no bionomic parameter ports.
#' Call [ramp.xds::LBionomics] at time zero after updating the parameter fields,
#' without calling [ramp.xds::LEffectSizes].
#'
#' Recognised names: `psi_E`, `phi_E`, `psi1`–`psi4`, `xi1`–`xi4`,
#' `phi1`–`phi4`, `theta1`–`theta4`, `phi_P`, `psi_P`.
#' @inheritSection stages Parameters
#' @inheritParams ramp.xds::change_L_pars
#' @seealso [dLdt.stages] or [make_L_obj_stages]
#' @return an **`xds`** model object
#' @keywords internal
#' @export
change_L_pars.stages <- function(xds_obj, options = list(), s = 1) {
  nHabitats <- xds_obj$nHabitats
  with(xds_obj$L_obj[[s]], with(options, {
    xds_obj$L_obj[[s]]$psi_E  <- checkIt(psi_E,  nHabitats)
    xds_obj$L_obj[[s]]$phi_E  <- checkIt(phi_E,  nHabitats)
    xds_obj$L_obj[[s]]$psi1   <- checkIt(psi1,   nHabitats)
    xds_obj$L_obj[[s]]$xi1    <- checkIt(xi1,    nHabitats)
    xds_obj$L_obj[[s]]$phi1   <- checkIt(phi1,   nHabitats)
    xds_obj$L_obj[[s]]$theta1 <- checkIt(theta1, nHabitats)
    xds_obj$L_obj[[s]]$psi2   <- checkIt(psi2,   nHabitats)
    xds_obj$L_obj[[s]]$xi2    <- checkIt(xi2,    nHabitats)
    xds_obj$L_obj[[s]]$phi2   <- checkIt(phi2,   nHabitats)
    xds_obj$L_obj[[s]]$theta2 <- checkIt(theta2, nHabitats)
    xds_obj$L_obj[[s]]$psi3   <- checkIt(psi3,   nHabitats)
    xds_obj$L_obj[[s]]$xi3    <- checkIt(xi3,    nHabitats)
    xds_obj$L_obj[[s]]$phi3   <- checkIt(phi3,   nHabitats)
    xds_obj$L_obj[[s]]$theta3 <- checkIt(theta3, nHabitats)
    xds_obj$L_obj[[s]]$psi4   <- checkIt(psi4,   nHabitats)
    xds_obj$L_obj[[s]]$xi4    <- checkIt(xi4,    nHabitats)
    xds_obj$L_obj[[s]]$phi4   <- checkIt(phi4,   nHabitats)
    xds_obj$L_obj[[s]]$theta4 <- checkIt(theta4, nHabitats)
    xds_obj$L_obj[[s]]$phi_P  <- checkIt(phi_P,  nHabitats)
    xds_obj$L_obj[[s]]$psi_P  <- checkIt(psi_P,  nHabitats)
    y <- get_inits(xds_obj, flatten=TRUE)
    xds_obj <- LBionomics(0, y, xds_obj, s)
    return(xds_obj)
  }))
}

#' @title Make initial values for `stages` (**L**)
#' @description Create a named list of initial values for the `stages` module.
#' Values in `options` override the defaults.
#' Initial values of \eqn{E, L_1, L_2, L_3, L_4, P} can be
#' set by passing the corresponding names in `options`.
#' @details Initial values must be numeric scalars or have length `nHabitats`.
#' Scalars are expanded to that length.
#' @param nHabitats the number of habitats in the model
#' @param options a [list] that overrides default values
#' @param E initial egg density
#' @param L1 initial first instar density
#' @param L2 initial second instar density
#' @param L3 initial third instar density
#' @param L4 initial fourth instar density
#' @param P initial pupal density
#' @return a named [list] of initial values
#' @keywords internal
#' @export
#' @seealso [setup_L_inits.stages], [change_L_inits.stages]
make_L_inits_stages = function(nHabitats, options = list(),
    E = 1, L1 = 1, L2 = 1, L3 = 1, L4 = 1, P = 1) {
  with(options, {
    list(
      E  = checkIt(E,  nHabitats),
      L1 = checkIt(L1, nHabitats),
      L2 = checkIt(L2, nHabitats),
      L3 = checkIt(L3, nHabitats),
      L4 = checkIt(L4, nHabitats),
      P  = checkIt(P,  nHabitats)
    )
  })
}

#' @title Set up initial values for `stages` (**L**)
#' @description Set up initial values for the `stages` module using
#' [make_L_inits_stages]. Values in `options` override the defaults.
#'
#' Initial values of \eqn{E, L_1, L_2, L_3, L_4, P} can be passed in `options`
#' by name.
#' @inheritParams ramp.xds::setup_L_inits
#' @seealso [make_L_inits_stages]
#' @return an **`xds`** model object
#' @keywords internal
#' @export
#' @seealso [change_L_inits.stages]
#' @note Initial values are stored at `xds_obj$L_obj[[s]]$inits`.
#' @seealso [ramp.xds::get_L_inits]
setup_L_inits.stages = function(xds_obj, options = list(), s=1){
  xds_obj$L_obj[[s]]$inits = make_L_inits_stages(xds_obj$nHabitats, options)
  return(xds_obj)
}

#' @title Change initial values for `stages` (**L**)
#' @description Change initial values named in `options` for the `stages`
#' module, leaving other values unchanged.
#'
#' Initial values of \eqn{E, L_1, L_2, L_3, L_4, P} are reset
#' if they are passed as named components of `options`.
#' @inheritParams ramp.xds::change_L_inits
#' @return an **`xds`** model object
#' @keywords internal
#' @export
#' @seealso [make_L_inits_stages], [setup_L_inits.stages]
#' @note Stored initial values for this module are at
#' `xds_obj$L_obj[[s]]$inits`.
change_L_inits.stages <- function(xds_obj, options = list(), s = 1) {
  with(xds_obj$L_obj[[s]]$inits, with(options, {
    xds_obj$L_obj[[s]]$inits$E  <- E
    xds_obj$L_obj[[s]]$inits$L1 <- L1
    xds_obj$L_obj[[s]]$inits$L2 <- L2
    xds_obj$L_obj[[s]]$inits$L3 <- L3
    xds_obj$L_obj[[s]]$inits$L4 <- L4
    xds_obj$L_obj[[s]]$inits$P  <- P
    return(xds_obj)
  }))
}

#' @title Set up variable indices for `stages` (**L**)
#' @description Assign consecutive variable indices for the `stages` module,
#' starting after `xds_obj$max_ix`, and update the maximum index.
#' Allocate blocks in the order `E_ix`, `L1_ix`, `L2_ix`, `L3_ix`, `L4_ix`, `P_ix`.
#' Each block has length `nHabitats`.
#' @inheritParams ramp.xds::setup_L_ix
#' @return an **`xds`** model object
#' @importFrom utils tail
#' @keywords internal
#' @export
#' @note Indices are stored at `xds_obj$L_obj[[s]]$ix`.
#' The maximum index is stored at `xds_obj$max_ix`.
setup_L_ix.stages <- function(xds_obj, s) {
  with(xds_obj, {
    E_ix  <- seq(from = max_ix + 1,          length.out = nHabitats)
    L1_ix <- seq(from = tail(E_ix,  1) + 1,  length.out = nHabitats)
    L2_ix <- seq(from = tail(L1_ix, 1) + 1,  length.out = nHabitats)
    L3_ix <- seq(from = tail(L2_ix, 1) + 1,  length.out = nHabitats)
    L4_ix <- seq(from = tail(L3_ix, 1) + 1,  length.out = nHabitats)
    P_ix  <- seq(from = tail(L4_ix, 1) + 1,  length.out = nHabitats)

    xds_obj$max_ix <- tail(P_ix, 1)
    xds_obj$L_obj[[s]]$ix <- list(
      E_ix  = E_ix,
      L1_ix = L1_ix,
      L2_ix = L2_ix,
      L3_ix = L3_ix,
      L4_ix = L4_ix,
      P_ix  = P_ix
    )
    return(xds_obj)
  })
}

#' @title Parse outputs for `stages` (**L**)
#' @description Extract the model state outputs for the `stages` module
#' from `outputs` and return them by name in a list.
#' Return `E`, `L1`, `L2`, `L3`, `L4`, and `P`.
#' @note Variable indices are stored at `xds_obj$L_obj[[s]]$ix`.
#' @inheritParams ramp.xds::parse_L_orbits
#' @return a [list]
#' @keywords internal
#' @export
parse_L_orbits.stages <- function(outputs, xds_obj, s) {
  with(xds_obj$L_obj[[s]]$ix, {
    list(
      E  = outputs[, E_ix],
      L1 = outputs[, L1_ix],
      L2 = outputs[, L2_ix],
      L3 = outputs[, L3_ix],
      L4 = outputs[, L4_ix],
      P  = outputs[, P_ix]
    )
  })
}

#' @title Check `stages` (**L**)
#' @description Run no consistency checks
#'
#' @inheritParams ramp.xds::check_L
#'
#' @return an **`xds`** model object
#' @keywords internal
#' @export
check_L.stages = function(xds_obj, s){
  return(xds_obj)
}

# the aquatic mosquito `stages` model

#' @title Get the skill set for `stages` (**L**)
#'
#' @description Return the capabilities and compatibilities of the `stages` module as a list.
#'
#' @inheritParams ramp.xds::get_skillset_L
#'
#' @return the **L** module skill set, as a list
#'
#' @keywords internal
#' @export
get_skillset_L.stages = function(xds_obj, s){
  skills = list(
    psi_is_port = FALSE,
    phi_is_port = FALSE,
    xi_is_port = FALSE,
    theta_is_port = FALSE,
    dts = FALSE,
    delay = FALSE
  )
  return(skills)
}

#' @title Compute steady states for `stages` (**L**)
#' @description Compute steady-state aquatic mosquito densities for the `stages`
#' module given the daily egg input `eta` in each habitat for the
#' differential-equation system in [dLdt.stages].
#' Uses [stats::nlm] to minimise the sum of squared derivatives.
#' @note This function does not use deSolve
#' @inheritParams ramp.xds::steady_state_L
#' @return a named [list] with elements `E`, `L1`, `L2`, `L3`, `L4`, `P`
#' @importFrom stats nlm
#' @keywords internal
#' @export
steady_state_L.stages = function(eta, xds_obj, s = 1) {
  with(xds_obj$L_obj[[s]], {
    n <- length(eta)
    obj_fn <- function(y, eta, L_obj) {
      with(L_obj, {
        E  <- y[1:n]
        L1 <- y[(n + 1):(2 * n)]
        L2 <- y[(2 * n + 1):(3 * n)]
        L3 <- y[(3 * n + 1):(4 * n)]
        L4 <- y[(4 * n + 1):(5 * n)]
        P  <- y[(5 * n + 1):(6 * n)]
        Lc <- w1 * L1 + w2 * L2 + w3 * L3 + L4
        m1 <- psi1 * L1 * exp(-xi1 * Lc)
        m2 <- psi2 * L2 * exp(-xi2 * Lc)
        m3 <- psi3 * L3 * exp(-xi3 * Lc)
        m4 <- psi4 * L4 * exp(-xi4 * Lc)
        dE  <- eta       - (psi_E + phi_E) * E
        dL1 <- psi_E * E - m1 - (phi1 + theta1 * Lc) * L1
        dL2 <- m1        - m2 - (phi2 + theta2 * Lc) * L2
        dL3 <- m2        - m3 - (phi3 + theta3 * Lc) * L3
        dL4 <- m3        - m4 - (phi4 + theta4 * Lc) * L4
        dP  <- m4        - (phi_P + psi_P) * P
        sum(c(dE, dL1, dL2, dL3, dL4, dP)^2)
      })
    }
    L_obj <- xds_obj$L_obj[[s]]
    y0    <- rep(eta, 6)
    y     <- abs(nlm(obj_fn, y0, L_obj = L_obj, eta = eta)$estimate)
    list(
      E  = y[1:n],
      L1 = y[(n + 1):(2 * n)],
      L2 = y[(2 * n + 1):(3 * n)],
      L3 = y[(3 * n + 1):(4 * n)],
      L4 = y[(4 * n + 1):(5 * n)],
      P  = y[(5 * n + 1):(6 * n)]
    )
  })
}
