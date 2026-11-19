# specialized methods for the adult mosquito BRQS model

#' @title The `BRQS` module for the **MY** component
#' @description
#' Implements the **MY** component using a BRQS (Blood-feeding / Resting /
#' gravid Queue / Sugar-feeding) model of adult mosquito ecology and infection
#' dynamics. Mosquitoes cycle through sugar-feeding (**S**), blood-feeding
#' (**B**), resting (**R**), and gravid (**Q**) states. Each behavioral state
#' is tracked as uninfected (**u**), exposed (**y**), or infectious (**z**).
#'
#' @section State Variables:
#' \describe{
#'   \item{`Su`, `Bu`, `Ru`, `Qu`}{uninfected mosquitoes in the sugar-feeding, blood-feeding, resting, and gravid states}
#'   \item{`Sy`, `By`, `Ry`, `Qy`}{exposed mosquitoes in each behavioral state}
#'   \item{`Sz`, `Bz`, `Rz`, `Qz`}{infectious mosquitoes in each behavioral state}
#'   \item{`S`, `B`, `R`, `Q`}{total mosquitoes in each behavioral state, derived by [get_MY_vars.BRQS]}
#'   \item{`M`}{total adult mosquito density, derived during parsing}
#'   \item{`Y`}{total exposed and infectious mosquito density, derived during parsing}
#'   \item{`Z`}{total infectious mosquito density, derived during parsing}
#' }
#'
#' @section Dynamical Terms:
#' \describe{
#'   \item{`Lambda` or \eqn{\Lambda}}{adult mosquito emergence rate, supplied by [ramp.xds::Emergence]}
#'   \item{`kappa` or \eqn{\kappa}}{net infectiousness, supplied by [ramp.xds::Transmission]}
#' }
#'
#' @section Parameters:
#' \describe{
#'   \item{`f`}{blood feeding rate --- use [ramp.xds::setup_f] or [ramp.xds::setup_blood_feeding_rates]}
#'   \item{`q`}{human fraction --- use [ramp.xds::setup_q] or [ramp.xds::setup_human_fractions]}
#'   \item{`g_b`, `g_q`, `g_s`, `g_r`}{mortality rates in the blood-feeding, gravid, sugar-feeding, and resting states}
#'   \item{`rho`}{rate of leaving the sugar-feeding state to blood-feed}
#'   \item{`zeta`}{rate at which blood-feeding mosquitoes take a sugar meal}
#'   \item{`xi`}{rate of leaving the resting state}
#'   \item{`eta`}{rate of transition from resting to gravid}
#'   \item{`nu`, `theta`}{parameters governing the egg-laying / gravid state}
#'   \item{`delta`}{fraction of gravid mosquitoes returning to blood feeding}
#'   \item{`omega`}{fraction of emerging mosquitoes entering the sugar-feeding state}
#'   \item{`tau`}{rate of progression from exposed to infectious}
#'   \item{`eip`}{extrinsic incubation period; \eqn{\tau} is its reciprocal}
#'   \item{`sigma_b`, `sigma_q`, `sigma_s`}{emigration rates in the blood-feeding, gravid, and sugar-feeding states}
#'   \item{`mu_b`, `mu_q`, `mu_s`}{fraction of emigrants lost in each dispersal state}
#'   \item{`K_b`, `K_q`, `K_s`}{dispersal matrices for the blood-feeding, gravid, and sugar-feeding states}
#'   \item{`eggsPerBatch`}{number of eggs per batch}
#' }
#'
#' After basic setup, [change_MY_pars.BRQS] can be used to change parameter
#' values. The model routes emergence between the blood-feeding and
#' sugar-feeding states using \eqn{\omega}; behavioral transition rates
#' \eqn{\rho}, \eqn{\zeta}, \eqn{\xi}, \eqn{\eta}, \eqn{\nu}, and \eqn{\theta}
#' govern movement among the four behavioral states.
#'
#' @section Survival and Dispersal:
#' Blood-feeding, gravid, and sugar-feeding mosquitoes have separate
#' demographic matrices, \eqn{\Omega_b}, \eqn{\Omega_q}, and \eqn{\Omega_s}.
#' Each combines its state-specific mortality, emigration rate, emigration-loss
#' fraction, and dispersal matrix.
#'
#' @section Blood Feeding:
#' Blood feeding occurs from the blood-feeding state at rate \eqn{f}. A
#' human blood meal transmits infection to a susceptible mosquito with
#' probability \eqn{\kappa}; the force of infection is \eqn{fq\kappa}.
#'
#' @section Dynamics:
#' The derivative method tracks the twelve state variables listed above.
#' Mosquitoes move among the four behavioral states, acquire infection during
#' human blood feeding, and progress from exposed to infectious at rate
#' \eqn{\tau}. The implementation is [dMYdt.BRQS].
#'
#' @section Infective Biting:
#' The patch-level human biting density, `fqM`, is based on blood-feeding
#' mosquitoes; infective biting density, `fqZ`, is based on infectious
#' blood-feeding mosquitoes. [F_fqM.BRQS] and [F_fqZ.BRQS] compute these
#' quantities. [ramp.xds::compute_EIR] combines `fqZ` with the biting
#' distribution matrix and resident-bite fraction to compute stratum-level
#' `EIR`; [ramp.xds::F_hbr] uses `fqM` to compute stratum-level `HBR`.
#'
#' @section Egg Laying:
#' Gravid mosquitoes lay eggs at rate \eqn{\nu}; egg production is the number
#' of gravid mosquitoes multiplied by the egg batch rate and eggs per batch:
#' \deqn{G=\nu(Q_u+Q_y+Q_z)\Gamma.}
#' [F_eggs.BRQS] computes \eqn{G}, which is passed to [ramp.xds::EggLaying]
#' for distribution among aquatic habitats.
#'
#' @seealso
#' + Basic setup: [make_MY_obj_BRQS], [setup_MY_obj.BRQS]
#' + Dynamics: [dMYdt.BRQS]
#' + Biting and eggs: [F_fqZ.BRQS], [F_fqM.BRQS], [F_eggs.BRQS]
#' + Parameters: [get_MY_pars.BRQS], [change_MY_pars.BRQS]
#' + Variables: [make_MY_inits_BRQS], [setup_MY_inits.BRQS], [parse_MY_orbits.BRQS]
#'
#' @name BRQS
#' @rdname BRQS
NULL

#' @title The skill set
#'
#' @inheritParams ramp.xds::get_skillset_MY
#'
#' @return the **MY** module skill set, as a list
#'
#' @keywords internal
#' @export
get_skillset_MY.BRQS = function(xds_obj, s){
  skills = list(
    delay = FALSE,
    not_implemented = TRUE
  )
  return(skills)
}

#' @title Compute derivatives for `BRQS` (**MY**)
#' @description Implements [dMYdt] for the BRQS ODE model.
#' @inheritSection BRQS State Variables
#' @inheritSection BRQS Dynamical Terms
#' @inheritSection BRQS Parameters
#' @inheritSection BRQS Dynamics
#' @inheritParams ramp.xds::dMYdt
#' @seealso [BRQS], [F_fqZ.BRQS], [F_eggs.BRQS]
#' @return a [numeric] vector
#' @keywords internal
#' @export
dMYdt.BRQS <- function(t, y, xds_obj, s){

  Lambda = as.vector(xds_obj$Lambda[[s]])
  kappa = as.vector(xds_obj$kappa[[s]])

  with(get_MY_vars(y, xds_obj, s),{
    with(xds_obj$MY_obj[[s]],{

      dBudt <- (1-omega)*Lambda + rho*Su + (nu+theta)*delta*Qu + xi*R - (f+zeta)*Bu - Omega_b %*% Bu
      dRudt <- f*(1-q*kappa)*Bu - (xi + eta + g)*Ru
      dQudt <- eta*Ru - (nu + theta)*Qu - Omega_q %*% Qu
      dSudt <- omega*Lambda + zeta*Bu + (nu+theta)*(1-delta)*Qu - rho*Su - Omega_s %*% Su

      dBydt <- rho*Sy + (nu + theta)*delta*Qy + xi*R - (f+zeta)*By - Omega_b %*% By - tau*By
      dRydt <- f*q*kappa*Bu + f*By - (xi + eta + g)*Ry - tau*Ry
      dQydt <- eta*Ry - (nu + theta)*Qy - Omega_q %*% Qy - tau*Qy
      dSydt <- zeta*By + (nu+theta)*(1-delta)*Qy - rho*Sy - Omega_s%*%Sy - tau*Sy

      dBzdt <- rho*Sz + (nu + theta)*delta*Qz + xi*R - (f+zeta)*Bz - Omega_b %*% Bz + tau*By
      dRzdt <- f*Bz - (xi + eta + g)*Rz + tau*Ry
      dQzdt <- eta*Rz - (nu + theta)*Qz - Omega_q %*% Qz + tau*Qy
      dSzdt <- zeta*Bz + (nu+theta)*(1-delta)*Qz -rho*Sz -Omega_s %*% Sz + tau*Sy

      return(c(dSudt, dBudt, dRudt, dQudt,
               dSydt, dBydt, dRydt, dQydt,
               dSzdt, dBzdt, dRzdt, dQzdt))
    })
  })
}

#' @title Return the variables as a list
#' @description This method dispatches on the type of `xds_obj$MY_obj[[s]]`
#' @inheritSection BRQS State Variables
#' @inheritParams ramp.xds::get_MY_vars
#' @return a [list]
#' @seealso [BRQS], [setup_MY_obj.BRQS]
#' @keywords internal
#' @export
get_MY_vars.BRQS<- function(y, xds_obj, s){
  with(xds_obj$MY_obj[[s]]$ix, return(list(
    Su <- y[Su_ix],
    Bu <- y[Bu_ix],
    Ru <- y[Ru_ix],
    Qu <- y[Qu_ix],
    Sy <- y[Sy_ix],
    By <- y[By_ix],
    Ry <- y[Ry_ix],
    Qy <- y[Qy_ix],
    Sz <- y[Sz_ix],
    Bz <- y[Bz_ix],
    Rz <- y[Rz_ix],
    Qz <- y[Qz_ix],
    B = Bu+By+Bz,
    R = Ru+Ry+Rz,
    Q = Qu+Qy+Qz,
    S = Su+Sy+Sz,
    M = B+R+Q+S
)))}

#' @title Mosquito bionomics for `BRQS` (**MY**)
#' @description Implements [MBionomics] for the BRQS model
#' @inheritSection BRQS Parameters
#' @inheritSection BRQS Survival and Dispersal
#' @inheritSection BRQS Blood Feeding
#' @inheritParams ramp.xds::MBionomics
#' @seealso [BRQS], [MEffectSizes.BRQS]
#' @return an **`xds`** model object
#' @keywords internal
#' @export
MBionomics.BRQS <- function(t, y, xds_obj, s) {with(xds_obj$MY_obj[[s]],{
  xds_obj$MY_obj[[s]]$es_f       <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_q       <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_nu      <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_g_s     <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_g_b     <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_g_r     <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_g_q     <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_sigma_b <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_sigma_q <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_sigma_s <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_mu_b    <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_mu_q    <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_mu_s    <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_theta   <- rep(1, nPatches)
  return(pars)
})}


#' @title Apply effect sizes for `BRQS` (**MY**)
#' @description Implements [MEffectSizes] for the BRQS model
#' @inheritSection BRQS Parameters
#' @inheritSection BRQS Survival and Dispersal
#' @inheritParams ramp.xds::MEffectSizes
#' @seealso [BRQS], [MBionomics.BRQS]
#' @return an **`xds`** model object
#' @keywords internal
#' @export
MEffectSizes.BRQS <- function(t, y, xds_obj, s) {with(xds_obj$MY_obj[[s]],{
  xds_obj$MY_obj[[s]]$f <- es_f*f_t
  xds_obj$MY_obj[[s]]$q <- es_q*q_t
  g <- es_g*g_t
  sigma_b <- es_sigma_b*sigma_b_t
  sigma_q <- es_sigma_q*sigma_q_t
  xds_obj$MY_obj[[s]]$g <- g
  xds_obj$MY_obj[[s]]$sigma_b <- sigma_b
  xds_obj$MY_obj[[s]]$sigma_q <- sigma_q
  xds_obj$MY_obj[[s]]$Omega_b = compute_Omega_xde(g, sigma_b, mu, calKb)
  xds_obj$MY_obj[[s]]$Omega_q = compute_Omega_xde(g, sigma_q, mu, calKq)
  return(pars)
})}

#' @title Infective mosquito biting rate for `BRQS`
#' @description Implements [F_fqZ] for the BRQS model.
#' @inheritSection BRQS Infective Biting
#' @inheritParams ramp.xds::F_fqZ
#' @seealso [BRQS], [F_fqM.BRQS]
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
F_fqZ.BRQS <- function(t, y, xds_obj, s) {
  f = get_f(xds_obj, s)
  q = get_f(xds_obj, s)
  fqZ = f*q*y[xds_obj$ix$MY[[s]]$Z_b_ix]
  return(fqZ)
}

#' @title Mosquito biting rate for `BRQS`
#' @description Returns `fqM`, the daily mosquito biting density on humans
#' for each patch in the BRQS model.
#' @inheritSection BRQS Infective Biting
#' @inheritParams ramp.xds::F_fqM
#' @seealso [BRQS], [ramp.xds::F_hbr]
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
F_fqM.BRQS <- function(t, y, xds_obj, s) {
  f = get_f(xds_obj, s)
  q = get_f(xds_obj, s)
  M = with(xds_obj$ix$MY[[s]], y[Bu_ix] + y[By_ix] + y[Bz_ix])
  fqM = f*q*M
  return(fqM)
}

#' @title Egg production rate for `BRQS`
#' @description Implements [F_eggs] for the BRQS model.
#' @inheritSection BRQS Egg Laying
#' @inheritParams ramp.xds::F_eggs
#' @seealso [BRQS], [ramp.xds::EggLaying]
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
F_eggs.BRQS <- function(t, y, xds_obj, s) {

  G <- with(xds_obj$ix$MY[[s]], y[Qu_ix] + y[Qy_ix] + y[Qz_ix])
  eggs = with(xds_obj$MY_obj[[s]], {
    return(G*nu*eggsPerBatch)
  })
  return(eggs)
}


#' @title Set up `BRQS` (**MY**)
#' @description
#'
#' This sets up a `BRQS` model object for the **MY** module for the \eqn{s^{th}} vector species
#' + the object is created by [make_MY_obj_BRQS]
#' + the initial conditions are set by [setup_MY_inits]
#' + `options` is passed to both to override the defaults
#'
#' @inheritSection BRQS State Variables
#' @inheritSection BRQS Parameters
#' @inheritSection BRQS Dynamics
#' @seealso [BRQS], [make_MY_obj_BRQS], [setup_MY_inits.BRQS]
#'
#' @note The module object is stored at `xds_obj$MY_obj[[s]]`
#'
#' @inheritParams ramp.xds::setup_MY_obj
#' @return a [list] vector
#' @keywords internal
#' @export
setup_MY_obj.BRQS = function(MYname, xds_obj, options=list(), s=1){
  xds_obj$MYname <- "BRQS"
  xds_obj$MY_obj[[s]] = make_MY_obj_BRQS(xds_obj$nPatches, options)
  xds_obj <- setup_F_circadian("setup", xds_obj, s=s)
  xds_obj$MY_obj[[s]]$skill_set <- get_skillset_MY(xds_obj, s)
  xds_obj <- setup_K_matrix("zero", xds_obj, options=list(which_K = "Kb"), s=s)
  xds_obj <- setup_K_matrix("zero", xds_obj, options=list(which_K = "Kq"), s=s)
  xds_obj <- setup_K_matrix("zero", xds_obj, options=list(which_K = "Ks"), s=s)
  xds_obj <- setup_MY_inits(xds_obj, options=options, s=s)
  return(xds_obj)
}

#' @title Make parameters for `BRQS` (**MY**)
#' @param nPatches is the number of patches, an integer
#' @param options a [list] of values that overrides the defaults
#' @param eip extrinsic incubation period
#' @param g mosquito mortality rate
#' @param sigma_b emigration rate while blood feeding
#' @param sigma_q emigration rate while egg laying
#' @param mu emigration loss
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @return a [list]
#' @keywords internal
#' @export
make_MY_obj_BRQS = function(nPatches, options=list(), eip=12,
                          g=1/12, sigma_b=1/8, sigma_q=1/8, mu=0, f=0.5, q=0.95,
                          nu=1, eggsPerBatch=60){

  with(options,{
    MY_obj <- list()
    class(MY_obj) <- "BRQS"

    MY_obj$nPatches <- nPatches

    eip_par <- list()
    class(eip_par) <- 'static'
    MY_obj$eip_par <- eip_par
    MY_obj$eip     <- eip

    MY_obj$f_t          <- checkIt(f, nPatches)
    MY_obj$es_f         <- rep(1, nPatches)

    MY_obj$nu_t         <- checkIt(nu, nPatches)
    MY_obj$es_nu        <- rep(1, nPatches)

    MY_obj$theta_t        <- checkIt(theta, nPatches)
    MY_obj$es_theta        <- rep(1, nPatches)

    MY_obj$g_b_t         <- checkIt(g_b, nPatches)
    MY_obj$es_g_b        <- rep(1, nPatches)

    MY_obj$g_q_t          <- checkIt(g_q, nPatches)
    MY_obj$es_g_q         <- rep(1, nPatches)

    MY_obj$g_s_t          <- checkIt(g_s, nPatches)
    MY_obj$es_g_s         <- rep(1, nPatches)

    MY_obj$g_r_t          <- checkIt(g_r, nPatches)
    MY_obj$es_g_r         <- rep(1, nPatches)

    MY_obj$sigma_b_t    <- checkIt(sigma_b, nPatches)
    MY_obj$es_sigma_b   <- rep(1, nPatches)

    MY_obj$sigma_q_t    <- checkIt(sigma_q, nPatches)
    MY_obj$es_sigma_q   <- rep(1, nPatches)

    MY_obj$sigma_s_t    <- checkIt(sigma_q, nPatches)
    MY_obj$es_sigma_s   <- rep(1, nPatches)

    MY_obj$mu_b_t    <- checkIt(mu_b, nPatches)
    MY_obj$es_mu_b   <- rep(1, nPatches)

    MY_obj$mu_q_t    <- checkIt(mu_q, nPatches)
    MY_obj$es_mu_q   <- rep(1, nPatches)

    MY_obj$mu_s_t    <- checkIt(mu_q, nPatches)
    MY_obj$es_mu_s   <- rep(1, nPatches)

    MY_obj$xi_t           <- checkIt(xi, nPatches)
    MY_obj$es_xi         <- rep(1, nPatches)

    MY_obj$eta_t           <- checkIt(eta, nPatches)
    MY_obj$es_eta         <- rep(1, nPatches)

    MY_obj$q_t          <- checkIt(q, nPatches)
    MY_obj$es_q         <- rep(1, nPatches)

    MY_obj$eggsPerBatch <- eggsPerBatch

    MY_obj$phi <- 1/MY_obj$eip


    Omega_par <- list()
    class(Omega_par) <- "static"
    MY_obj$Omega_par <- Omega_par
    MY_obj$Omega_b <- with(MY_obj, compute_Omega_xde(g, sigma_b, mu, calK))
    MY_obj$Omega_q <- with(MY_obj, compute_Omega_xde(g, sigma_q, mu, calK))
    MY_obj$Omega_s <- with(MY_obj, compute_Omega_xde(g, sigma_s, mu, calK))
    base <- 'BRQS'
    class(base) <- 'BRQS'
    MY_obj$baseline <- base

    return(MY_obj)
})}

#' @title Return the parameters as a list
#' @description This method dispatches on the type of `xds_obj$MY_obj[[s]]`.
#' @inheritSection BRQS Parameters
#' @param xds_obj an **`xds`** model object
#' @param s the vector species index
#' @return a [list]
#' @seealso [BRQS], [change_MY_pars.BRQS]
#' @keywords internal
#' @export
get_MY_pars.BRQS <- function(xds_obj, s=1) {
  with(xds_obj$MY_obj[[s]], list(
    f=f_t, q=q_t, g=g_t, sigma_b=sigma_b_t, sigma_q = sigma_q_t, eip=eip, mu=mu_t,
    nu=nu_t, eggsPerBatch=eggsPerBatch, calK=calK
  ))
}

#' @title Change parameters for `BRQS` (**MY**)
#' @description For the `BRQS` module (**MY**) for the \eqn{s^{th}} species,
#' change parameter values named in `options` leaving other values unchanged.
#' @inheritSection BRQS Parameters
#' @inheritParams ramp.xds::change_MY_pars
#' @return an **`xds`** model object
#' @keywords internal
#' @export
change_MY_pars.BRQS <- function(xds_obj, options=list(), s=1) {
  nHabitats <- xds_obj$nHabitats
  with(xds_obj$MY_obj[[s]], with(options,{

    xds_obj$MY_obj[[s]]$f_obj$f = f
    xds_obj$MY_obj[[s]]$q_obj$q = q
    xds_obj$MY_obj[[s]]$g_obj$g = g
    xds_obj$MY_obj[[s]]$sigma_b_obj$sigma_b = sigma_b
    xds_obj$MY_obj[[s]]$sigma_q_obj$sigma_q = sigma_q
    xds_obj$MY_obj[[s]]$eip_obj$eip = eip
    xds_obj$MY_obj[[s]]$mu_obj$mu = mu
    xds_obj$MY_obj[[s]]$nu_obj$nu = nu
    xds_obj$MY_obj[[s]]$eggsPerBatch = eggsPerBatch

    y <- get_inits(xds_obj, flatten=TRUE)
    xds_obj <- MBionomics(0, y, xds_obj, s)
    xds_obj <- MEffectSizes(0, y, xds_obj, s)


    xds_obj$MY_obj[[s]]$f_t = f
    xds_obj$MY_obj[[s]]$q_t = q
    xds_obj$MY_obj[[s]]$nu_t = nu
    xds_obj$MY_obj[[s]]$g_s_t = g_s
    xds_obj$MY_obj[[s]]$g_b_t = g_b
    xds_obj$MY_obj[[s]]$g_r_t = g_r
    xds_obj$MY_obj[[s]]$g_q_t = g_q
    xds_obj$MY_obj[[s]]$mu_s_t = mu_s
    xds_obj$MY_obj[[s]]$mu_b_t = mu_b
    xds_obj$MY_obj[[s]]$mu_q_t = mu_q
    xds_obj$MY_obj[[s]]$sigma_b_t = sigma_b
    xds_obj$MY_obj[[s]]$sigma_q_t = sigma_q
    xds_obj$MY_obj[[s]]$sigma_s_t = sigma_s
    xds_obj$MY_obj[[s]]$eip = eip
    xds_obj$MY_obj[[s]]$eggsPerBatch = eggsPerBatch
    return(pars)
  }))}

#' @title Setup initial values for `BRQS` (**MY**)
#' @description Set initial values for all the variables in the `BRQS` module
#' (**MY**) for the \eqn{s^{th}} species.
#' + [make_MY_inits_BRQS] supplies default values
#' + The values in `options` override the defaults
#' + To inspect variable names, use [get_MY_inits] or
#' see the inline module documentation.
#' @inheritSection BRQS State Variables
#' @inheritParams ramp.xds::setup_MY_inits
#' @return a [list]
#' @keywords internal
#' @export
setup_MY_inits.BRQS = function(xds_obj, options=list(), s=1){
  xds_obj$MYinits[[s]] = with(xds_obj$MY_obj[[s]], make_MY_inits_BRQS(nPatches, options))
  return(xds_obj)
}

#' @title Make initial values for `BRQS` (**MY**)
#' @inheritSection BRQS State Variables
#' @seealso [BRQS], [setup_MY_inits.BRQS]
#' @param nPatches the number of patches
#' @param options a [list] of values that overrides the defaults
#' @param Su uninfected, sugar feeding mosquitoes
#' @param Bu uninfected, blood feeding mosquitoes
#' @param Ru uninfected, resting mosquitoes
#' @param Qu uninfected, gravid mosqutioes
#' @param Sy infected, sugar feeding mosquitoes
#' @param By infected, blood feeding mosquitoes
#' @param Ry infected, resting mosquitoes
#' @param Qy infected, gravid mosqutioes
#' @param Sz infectious, sugar feeding mosquitoes
#' @param Bz infectious, blood feeding mosquitoes
#' @param Rz infectious, resting mosquitoes
#' @param Qz infectious, gravid mosqutioes
#' @return a [list]
#' @keywords internal
#' @export
make_MY_inits_BRQS = function(nPatches, options = list(),
                              Su=100, Bu=100, Ru=1, Qu=10,
                              Sy=0, By=0, Ry=0, Qy=0,
                              Sz=0, Bz=0, Rz=0, Qz=0){

  with(options,{
    Su = checkit(Su, nPatches)
    Bu = checkit(Bu, nPatches)
    Ru = checkit(Ru, nPatches)
    Qu = checkit(Qu, nPatches)
    Sy = checkit(Sy, nPatches)
    By = checkit(By, nPatches)
    Ry = checkit(Ry, nPatches)
    Qy = checkit(Qy, nPatches)
    Sz = checkit(Sz, nPatches)
    Bz = checkit(Bz, nPatches)
    Rz = checkit(Rz, nPatches)
    Qz = checkit(Qz, nPatches)

    return(list(Su=Su, Bu=Bu, Ru=Ru, Qu=Qu,
                Sy=Sy, By=By, Ry=Ry, Qy=Qy,
                Sz=Sz, Bz=Bz, Rz=Rz, Qz=Qz))
  })}


#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [setup_MY_ix] for the BRQS model.
#' @inheritSection BRQS State Variables
#' @inheritParams ramp.xds::setup_MY_ix
#' @return none
#' @importFrom utils tail
#' @keywords internal
#' @export
setup_MY_ix.BRQS <- function(xds_obj, s) {with(xds_obj,{

  Su_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Su_ix, 1)

  Bu_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Bu_ix, 1)

  Ru_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Ru_ix, 1)

  Qu_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Qu_ix, 1)

  Sy_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Sy_ix, 1)

  By_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(By_ix, 1)

  Ry_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Ry_ix, 1)

  Qy_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Qy_ix, 1)

  Sz_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Sz_ix, 1)

  Bz_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Bz_ix, 1)

  Rz_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Rz_ix, 1)

  Qz_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Qz_ix, 1)

  xds_obj$max_ix = max_ix
  xds_obj$MY[[s]]$ix = list(
    Su_ix = Su_ix, Bu_ix = Bu_ix, Ru_ix=Ru_ix, Qu_ix = Qu_ix,
    Sy_ix = Sy_ix, By_ix = By_ix, Ry_ix=Ry_ix, Qy_ix = Qy_ix,
    Sz_ix = Sz_ix, Bz_ix = Bz_ix, Rz_ix=Rz_ix, Qz_ix = Qz_ix)


  return(pars)
})}


#' @title Parse outputs for `BRQS` (**MY**)
#' @description Implements [parse_MY_orbits] for the BRQS model
#' @inheritSection BRQS State Variables
#' @inheritSection BRQS Infective Biting
#' @inheritParams ramp.xds::parse_MY_orbits
#' @seealso [BRQS], [get_MY_vars.BRQS]
#' @return none
#' @keywords internal
#' @export
parse_MY_orbits.BRQS <- function(outputs, xds_obj, s) {with(xds_obj$ix$MY[[s]],{
  Su = outputs[,Su_ix]
  Bu = outputs[,Bu_ix]
  Ru = outputs[,Ru_ix]
  Qu = outputs[,Qu_ix]
  Sy = outputs[,Sy_ix]
  By = outputs[,By_ix]
  Ry = outputs[,Ry_ix]
  Qy = outputs[,Qy_ix]
  Sz = outputs[,Sz_ix]
  Bz = outputs[,Bz_ix]
  Rz = outputs[,Rz_ix]
  Qz = outputs[,Qz_ix]
  M =  Su+Bu+Ru+Qu+Sy+By+Ry+Qy+Sz+Bz+Rz+Qz
  Y =  Sy+By+Ry+Qy+Sz+Bz+Rz+Qz
  Z =  Sz+Bz+Rz+Qz
  y = Y/M
  z = Z/M
  return(list(
    Su=Su, Bu=Bu, Ru=Ru, Qu=Qu,
    Sy=Sy, By=By, Ry=Ry, Qy=Qy,
    Sz=Sz, Bz=Bz, Rz=Rz, Qz=Qz, M=M, Y=Y, Z=Z, y=y, z=z))
})}

#' @title Change initial values for `BRQS` (**MY**)
#' @description For the `BRQS` module (**MY**) for the \eqn{s^{th}} species,
#' change initial values named in `options` leaving other values unchanged.
#' @inheritSection BRQS State Variables
#' @seealso [BRQS], [setup_MY_inits.BRQS]
#' @inheritParams ramp.xds::change_MY_inits
#' @return none
#' @keywords internal
#' @export
change_MY_inits.BRQS <- function(xds_obj, options = list(), s=1) {
  inits = with(get_MY_inits(xds_obj, s), with(options,
    list(Su=Su, Bu=Bu, Ru=Ru, Qu=Qu,
         Sy=Sy, By=By, Ry=Ry, Qy=Qy,
         Sz=Sz, Bz=Bz, Rz=Rz, Qz=Qz)))

  xds_obj$MY_obj[[s]]$inits=inits
  return(xds_obj)
}


#' @title Get the feeding rate
#' @param xds_obj an **`xds`** model object
#' @param s the vector species index
#' @return a [numeric] vector
#' @keywords internal
#' @export
get_f.BRQS = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], f_t*es_f)
}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** model object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @keywords internal
#' @export
get_q.BRQS = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], q_t*es_q)
}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** model object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @keywords internal
#' @export
get_g.BRQS = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], g_t*es_g)
}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** model object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @keywords internal
#' @export
get_sigma.BRQS = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], sigma_b_t*es_sigma_b)
}
