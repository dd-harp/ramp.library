
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

# the aquatic mosquito `stages` model

#' @title The **L** module skill set
#'
#' @description The **L** skill set is a list of
#' a module's capabilities
#'
#' @param Lname the name of the **L** module
#'
#' @return *L* module skill set, as a list
#'
#' @keywords internal
#' @export
skill_set_L.stages = function(Lname = "stages"){
  list(trivial=FALSE)
}

#' @title Check the `stages` module
#' @description Run no consistency checks
#'
#' @param xds_obj an **`xds`** model object
#' @param s the vector species index
#'
#' @return an **`xds`** object
#' @keywords internal
#' @export
check_L.stages = function(xds_obj, s){
  return(xds_obj)
}

#' @title Compute derivatives for `stages` (**L**)
#'
#' @description
#' Implements the staged aquatic development model with four larval instars.
#' Density-dependent effects act through the weighted mean crowding index
#' \eqn{L = w_1 L_1 + w_2 L_2 + w_3 L_3 + L_4.}
#'
#' **Variables:**
#'
#' - \eqn{E}: egg density per habitat
#' - \eqn{L_1, L_2, L_3, L_4}: larval instar densities per habitat
#' - \eqn{P}: pupal density per habitat
#'
#' **Input Term:**
#'
#' - \eqn{\eta} or `eta`: egg deposition rate (from [F_eggs])
#'
#' **Dynamical System:**
#'
#' \deqn{dE/dt = \eta - (\psi_E + \phi_E) E}
#' \deqn{dL_1/dt = \psi_E E - \psi_1 L_1 e^{-\xi_1 L} - (\phi_1 + \theta_1 L) L_1}
#' \deqn{dL_2/dt = \psi_1 L_1 e^{-\xi_1 L} - \psi_2 L_2 e^{-\xi_2 L} - (\phi_2 + \theta_2 L) L_2}
#' \deqn{dL_3/dt = \psi_2 L_2 e^{-\xi_2 L} - \psi_3 L_3 e^{-\xi_3 L} - (\phi_3 + \theta_3 L) L_3}
#' \deqn{dL_4/dt = \psi_3 L_3 e^{-\xi_3 L} - \psi_4 L_4 e^{-\xi_4 L} - (\phi_4 + \theta_4 L) L_4}
#' \deqn{dP/dt = \psi_4 L_4 e^{-\xi_4 L} - (\phi_P + \psi_P) P}
#'
#' **Output Term:**
#'
#' - The function [F_emerge] computes the net emergence rate (\eqn{\alpha}):
#'
#' \deqn{\alpha = \psi_P P}
#'
#' @inheritParams ramp.xds::dLdt
#' @return a [numeric] vector
#' @seealso [make_L_obj_stages]
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

#' @title Compute emergent adults for `stages` (**L** component)
#' @description The number of female adults emerging from the habitats,
#' per day, is:
#' \deqn{\psi_P P/2 }
#' @inheritParams ramp.xds::F_emerge
#' @return a [numeric] vector of length `nHabitats`
#' @seealso [dLdt.stages]
#' @keywords internal
#' @export
F_emerge.stages <- function(t, y, xds_obj, s) {
  P <- y[xds_obj$L_obj[[s]]$ix$P_ix]
  with(xds_obj$L_obj[[s]], {
    return(psi_P * P/2)
  })
}

#' @title Mosquito bionomics for `stages` (**L**)
#'
#' @description Resets all effect sizes
#' (`es_psi_E`, `es_phi_E`, `es_psi1`–`es_psi4`, `es_xi1`–`es_xi4`,
#' `es_phi1`–`es_phi4`, `es_theta1`–`es_theta4`, `es_phi_P`, `es_psi_P`)
#' to 1. The time-varying baseline values (`_t` fields) are set by
#' [change_L_pars.stages] and persist until explicitly updated.
#' @inheritParams ramp.xds::LBionomics
#'
#' @return an **`xds`** object
#' @keywords internal
#' @export
LBionomics.stages <- function(t, y, xds_obj, s) {
  with(xds_obj$L_obj[[s]], {
    xds_obj$L_obj[[s]]$es_psi_E  <- rep(1, nHabitats)
    xds_obj$L_obj[[s]]$es_phi_E  <- rep(1, nHabitats)
    xds_obj$L_obj[[s]]$es_psi1   <- rep(1, nHabitats)
    xds_obj$L_obj[[s]]$es_xi1    <- rep(1, nHabitats)
    xds_obj$L_obj[[s]]$es_phi1   <- rep(1, nHabitats)
    xds_obj$L_obj[[s]]$es_theta1 <- rep(1, nHabitats)
    xds_obj$L_obj[[s]]$es_psi2   <- rep(1, nHabitats)
    xds_obj$L_obj[[s]]$es_xi2    <- rep(1, nHabitats)
    xds_obj$L_obj[[s]]$es_phi2   <- rep(1, nHabitats)
    xds_obj$L_obj[[s]]$es_theta2 <- rep(1, nHabitats)
    xds_obj$L_obj[[s]]$es_psi3   <- rep(1, nHabitats)
    xds_obj$L_obj[[s]]$es_xi3    <- rep(1, nHabitats)
    xds_obj$L_obj[[s]]$es_phi3   <- rep(1, nHabitats)
    xds_obj$L_obj[[s]]$es_theta3 <- rep(1, nHabitats)
    xds_obj$L_obj[[s]]$es_psi4   <- rep(1, nHabitats)
    xds_obj$L_obj[[s]]$es_xi4    <- rep(1, nHabitats)
    xds_obj$L_obj[[s]]$es_phi4   <- rep(1, nHabitats)
    xds_obj$L_obj[[s]]$es_theta4 <- rep(1, nHabitats)
    xds_obj$L_obj[[s]]$es_phi_P  <- rep(1, nHabitats)
    xds_obj$L_obj[[s]]$es_psi_P  <- rep(1, nHabitats)
    return(xds_obj)
  })
}

#' @title Apply effect sizes for `stages` (**L**)
#' @description Implements [LEffectSizes] for the `stages` model.
#' Computes effective parameter values as the product of each time-varying
#' baseline (`_t`) and its corresponding effect size (`es_`).
#' @inheritParams ramp.xds::LEffectSizes
#' @return an **`xds`** object
#' @keywords internal
#' @export
LEffectSizes.stages <- function(t, y, xds_obj, s) {
  with(xds_obj$L_obj[[s]], {
    xds_obj$L_obj[[s]]$psi_E  <- psi_E_t  * es_psi_E
    xds_obj$L_obj[[s]]$phi_E  <- phi_E_t  * es_phi_E
    xds_obj$L_obj[[s]]$psi1   <- psi1_t   * es_psi1
    xds_obj$L_obj[[s]]$xi1    <- xi1_t    * es_xi1
    xds_obj$L_obj[[s]]$phi1   <- phi1_t   * es_phi1
    xds_obj$L_obj[[s]]$theta1 <- theta1_t * es_theta1
    xds_obj$L_obj[[s]]$psi2   <- psi2_t   * es_psi2
    xds_obj$L_obj[[s]]$xi2    <- xi2_t    * es_xi2
    xds_obj$L_obj[[s]]$phi2   <- phi2_t   * es_phi2
    xds_obj$L_obj[[s]]$theta2 <- theta2_t * es_theta2
    xds_obj$L_obj[[s]]$psi3   <- psi3_t   * es_psi3
    xds_obj$L_obj[[s]]$xi3    <- xi3_t    * es_xi3
    xds_obj$L_obj[[s]]$phi3   <- phi3_t   * es_phi3
    xds_obj$L_obj[[s]]$theta3 <- theta3_t * es_theta3
    xds_obj$L_obj[[s]]$psi4   <- psi4_t   * es_psi4
    xds_obj$L_obj[[s]]$xi4    <- xi4_t    * es_xi4
    xds_obj$L_obj[[s]]$phi4   <- phi4_t   * es_phi4
    xds_obj$L_obj[[s]]$theta4 <- theta4_t * es_theta4
    xds_obj$L_obj[[s]]$phi_P  <- phi_P_t  * es_phi_P
    xds_obj$L_obj[[s]]$psi_P  <- psi_P_t  * es_psi_P
    return(xds_obj)
  })
}

#' @title Set up `stages` (**L**)
#' @description The function sets up `L_obj` for the \eqn{s^{th}} species
#' by calling [make_L_obj_stages]
#' @inheritParams ramp.xds::setup_L_obj
#' @return an **`xds`** object
#' @seealso [make_L_obj_stages]
#' @keywords internal
#' @export
setup_L_obj.stages = function(Lname, xds_obj, s, options=list()){
  L_obj <- make_L_obj_stages(xds_obj$nHabitats, options)
  class(L_obj) <- c("stages", paste("stages_", xds_obj$xds, sep=""))
  xds_obj$L_obj[[s]] = L_obj
  xds_obj <- LBionomics(0, 0, xds_obj, s)
  return(xds_obj)
}

#' @title Make `L_obj` for `stages` (**L** component)
#' @description The following parameters will be set to the values in
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

    # effective parameter values (updated each step by LEffectSizes)
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

    # time-varying baselines (updated by change_L_pars)
    L_obj$psi_E_t  <- L_obj$psi_E
    L_obj$phi_E_t  <- L_obj$phi_E
    L_obj$psi1_t   <- L_obj$psi1
    L_obj$xi1_t    <- L_obj$xi1
    L_obj$phi1_t   <- L_obj$phi1
    L_obj$theta1_t <- L_obj$theta1
    L_obj$psi2_t   <- L_obj$psi2
    L_obj$xi2_t    <- L_obj$xi2
    L_obj$phi2_t   <- L_obj$phi2
    L_obj$theta2_t <- L_obj$theta2
    L_obj$psi3_t   <- L_obj$psi3
    L_obj$xi3_t    <- L_obj$xi3
    L_obj$phi3_t   <- L_obj$phi3
    L_obj$theta3_t <- L_obj$theta3
    L_obj$psi4_t   <- L_obj$psi4
    L_obj$xi4_t    <- L_obj$xi4
    L_obj$phi4_t   <- L_obj$phi4
    L_obj$theta4_t <- L_obj$theta4
    L_obj$phi_P_t  <- L_obj$phi_P
    L_obj$psi_P_t  <- L_obj$psi_P

    # effect sizes (reset to 1 by LBionomics each step)
    L_obj$es_psi_E  <- rep(1, nHabitats)
    L_obj$es_phi_E  <- rep(1, nHabitats)
    L_obj$es_psi1   <- rep(1, nHabitats)
    L_obj$es_xi1    <- rep(1, nHabitats)
    L_obj$es_phi1   <- rep(1, nHabitats)
    L_obj$es_theta1 <- rep(1, nHabitats)
    L_obj$es_psi2   <- rep(1, nHabitats)
    L_obj$es_xi2    <- rep(1, nHabitats)
    L_obj$es_phi2   <- rep(1, nHabitats)
    L_obj$es_theta2 <- rep(1, nHabitats)
    L_obj$es_psi3   <- rep(1, nHabitats)
    L_obj$es_xi3    <- rep(1, nHabitats)
    L_obj$es_phi3   <- rep(1, nHabitats)
    L_obj$es_theta3 <- rep(1, nHabitats)
    L_obj$es_psi4   <- rep(1, nHabitats)
    L_obj$es_xi4    <- rep(1, nHabitats)
    L_obj$es_phi4   <- rep(1, nHabitats)
    L_obj$es_theta4 <- rep(1, nHabitats)
    L_obj$es_phi_P  <- rep(1, nHabitats)
    L_obj$es_psi_P  <- rep(1, nHabitats)

    return(L_obj)
  })
}

#' @title Get parameters for `stages` (**L**)
#' @description Get the **L** component parameters
#' @param xds_obj an **`xds`** model object
#' @param s the vector species index
#' @return a [list]
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
#' @description Set the values of **L** component parameters.
#' Named values passed in `options` update the corresponding
#' time-varying baseline (`_t`) fields, which take effect at the
#' next call to [LEffectSizes.stages].
#'
#' Recognised names: `psi_E`, `phi_E`, `psi1`–`psi4`, `xi1`–`xi4`,
#' `phi1`–`phi4`, `theta1`–`theta4`, `phi_P`, `psi_P`.
#' @inheritParams ramp.xds::change_L_pars
#' @seealso [dLdt.stages] or [make_L_obj_stages]
#' @return an **`xds`** object
#' @keywords internal
#' @export
change_L_pars.stages <- function(xds_obj, s = 1, options = list()) {
  nHabitats <- xds_obj$nHabitats
  with(xds_obj$L_obj[[s]], with(options, {
    xds_obj$L_obj[[s]]$psi_E_t  <- checkIt(psi_E,  nHabitats)
    xds_obj$L_obj[[s]]$phi_E_t  <- checkIt(phi_E,  nHabitats)
    xds_obj$L_obj[[s]]$psi1_t   <- checkIt(psi1,   nHabitats)
    xds_obj$L_obj[[s]]$xi1_t    <- checkIt(xi1,    nHabitats)
    xds_obj$L_obj[[s]]$phi1_t   <- checkIt(phi1,   nHabitats)
    xds_obj$L_obj[[s]]$theta1_t <- checkIt(theta1, nHabitats)
    xds_obj$L_obj[[s]]$psi2_t   <- checkIt(psi2,   nHabitats)
    xds_obj$L_obj[[s]]$xi2_t    <- checkIt(xi2,    nHabitats)
    xds_obj$L_obj[[s]]$phi2_t   <- checkIt(phi2,   nHabitats)
    xds_obj$L_obj[[s]]$theta2_t <- checkIt(theta2, nHabitats)
    xds_obj$L_obj[[s]]$psi3_t   <- checkIt(psi3,   nHabitats)
    xds_obj$L_obj[[s]]$xi3_t    <- checkIt(xi3,    nHabitats)
    xds_obj$L_obj[[s]]$phi3_t   <- checkIt(phi3,   nHabitats)
    xds_obj$L_obj[[s]]$theta3_t <- checkIt(theta3, nHabitats)
    xds_obj$L_obj[[s]]$psi4_t   <- checkIt(psi4,   nHabitats)
    xds_obj$L_obj[[s]]$xi4_t    <- checkIt(xi4,    nHabitats)
    xds_obj$L_obj[[s]]$phi4_t   <- checkIt(phi4,   nHabitats)
    xds_obj$L_obj[[s]]$theta4_t <- checkIt(theta4, nHabitats)
    xds_obj$L_obj[[s]]$phi_P_t  <- checkIt(phi_P,  nHabitats)
    xds_obj$L_obj[[s]]$psi_P_t  <- checkIt(psi_P,  nHabitats)
    return(xds_obj)
  }))
}

#' @title Setup initial values for `stages` (**L**)
#' @description This sets initial values of the variables \eqn{E, L_1, L_2,
#' L_3, L_4, P} by calling [make_L_inits_stages]. Default values are used
#' unless other values are passed in `options` by name.
#' @inheritParams ramp.xds::setup_L_inits
#' @seealso [make_L_inits_stages]
#' @return an **`xds`** object
#' @keywords internal
#' @export
setup_L_inits.stages = function(xds_obj, s, options = list()){
  xds_obj$L_obj[[s]]$inits = make_L_inits_stages(xds_obj$nHabitats, options)
  return(xds_obj)
}

#' @title Make initial values for `stages` (**L** component)
#' @description Initial values of \eqn{E, L_1, L_2, L_3, L_4, P} can be
#' set by passing the corresponding names in `options`.
#' @param nHabitats the number of habitats in the model
#' @param options a [list] that overwrites default values
#' @param E initial egg density
#' @param L1 initial first instar density
#' @param L2 initial second instar density
#' @param L3 initial third instar density
#' @param L4 initial fourth instar density
#' @param P initial pupal density
#' @return a [list] with initial conditions
#' @keywords internal
#' @export
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

#' @title List variables for `stages` (**L**)
#' @description Extract the **L** component variables from the
#' vector of state variables (`y`) and return them as a named list.
#' @inheritParams ramp.xds::get_L_vars
#' @return a named [list]
#' @keywords internal
#' @export
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

#' @title Change initial values for `stages` (**L**)
#' @description Initial values of \eqn{E, L_1, L_2, L_3, L_4, P} are reset
#' if they are passed as named components of `options`.
#' @inheritParams ramp.xds::change_L_inits
#' @return an **`xds`** object
#' @keywords internal
#' @export
change_L_inits.stages <- function(xds_obj, s = 1, options = list()) {
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

#' @title Setup variable indices for `stages` (**L** component)
#' @description Set the values of the indices for the **L** component
#' variables for the `stages` module. Indices are allocated in the order
#' \eqn{E, L_1, L_2, L_3, L_4, P}, each block of length `nHabitats`.
#' @inheritParams ramp.xds::setup_L_ix
#' @return an **`xds`** object
#' @importFrom utils tail
#' @keywords internal
#' @export
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
#' @description Returns the columns representing \eqn{E, L_1, L_2, L_3,
#' L_4, P} from a matrix where each column is a state variable. The
#' variables are returned as a named list.
#' @inheritParams ramp.xds::parse_L_orbits
#' @return a named [list]
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

#' @title Compute the steady state of `dLdt.stages` (**L** component)
#' @description Given an egg deposition rate `eta`,
#' return a steady state value for the equations in [dLdt.stages].
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
