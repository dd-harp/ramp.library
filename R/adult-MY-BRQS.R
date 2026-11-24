# specialized methods for the adult mosquito BRQS model

# TODO: Reconcile constructor inputs with the state-specific parameter fields;
# several referenced parameters are not defined by its arguments or defaults.
# Confirm the intended sigma_s and mu_s inputs: the constructor currently
# initializes them from sigma_q and mu_q.
# TODO: Align get_MY_pars.BRQS and constructor demographic matrices with the
# state-specific parameter fields; both still reference legacy generic inputs.
# TODO: Review F_fqZ.BRQS and F_fqM.BRQS: both obtain q using get_f, not get_q.
# TODO: Review dMYdt.BRQS: all three blood-feeding equations use total R in
# xi*R rather than the corresponding Ru, Ry, or Rz. Confirm the intended flows.
# TODO: Reconcile dMYdt.BRQS's generic g in resting-state losses with the
# state-specific mortality parameters, including g_r.
# TODO: Align dMYdt.BRQS's top-level Lambda and kappa reads with the shared
# dynamical terms at xds_obj$terms$Lambda and xds_obj$terms$kappa.

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
#' of gravid mosquitoes multiplied by the oviposition rate and eggs per batch:
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

#' @title Compute the derivatives for `BRQS` (**MY**)
#' @description Compute the derivatives for the `BRQS` module:
#' \deqn{\begin{array}{rl}
#' dS_u/dt &= \omega\Lambda+\zeta B_u+(\nu+\theta)(1-\delta)Q_u-\rho S_u-\Omega_s\cdot S_u \\
#' dB_u/dt &= (1-\omega)\Lambda+\rho S_u+(\nu+\theta)\delta Q_u+\xi R-(f+\zeta)B_u-\Omega_b\cdot B_u \\
#' dR_u/dt &= f(1-q\kappa)B_u-(\xi+\eta+g)R_u \\
#' dQ_u/dt &= \eta R_u-(\nu+\theta)Q_u-\Omega_q\cdot Q_u \\
#' dS_y/dt &= \zeta B_y+(\nu+\theta)(1-\delta)Q_y-\rho S_y-\Omega_s\cdot S_y-\tau S_y \\
#' dB_y/dt &= \rho S_y+(\nu+\theta)\delta Q_y+\xi R-(f+\zeta)B_y-\Omega_b\cdot B_y-\tau B_y \\
#' dR_y/dt &= fq\kappa B_u+fB_y-(\xi+\eta+g)R_y-\tau R_y \\
#' dQ_y/dt &= \eta R_y-(\nu+\theta)Q_y-\Omega_q\cdot Q_y-\tau Q_y \\
#' dS_z/dt &= \zeta B_z+(\nu+\theta)(1-\delta)Q_z-\rho S_z-\Omega_s\cdot S_z+\tau S_y \\
#' dB_z/dt &= \rho S_z+(\nu+\theta)\delta Q_z+\xi R-(f+\zeta)B_z-\Omega_b\cdot B_z+\tau B_y \\
#' dR_z/dt &= fB_z-(\xi+\eta+g)R_z+\tau R_y \\
#' dQ_z/dt &= \eta R_z-(\nu+\theta)Q_z-\Omega_q\cdot Q_z+\tau Q_y
#' \end{array}}
#' @inheritSection BRQS State Variables
#' @inheritSection BRQS Dynamical Terms
#' @inheritSection BRQS Parameters
#' @inheritSection BRQS Dynamics
#' @inheritParams ramp.xds::dMYdt
#' @seealso [BRQS], [F_fqZ.BRQS], [F_eggs.BRQS]
#' @note Derivatives are returned in blocks for `Su`, `Bu`, `Ru`, `Qu`,
#' `Sy`, `By`, `Ry`, `Qy`, `Sz`, `Bz`, `Rz`, and `Qz`, each in patch order.
#' @return a [numeric] vector of derivatives in the module's state-variable order
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

#' @title Get the state variables for `BRQS` (**MY**)
#' @description Extract the state variables for the `BRQS` module
#' from `y` and return them by name in a list for use during model solving.
#' Also compute totals `B`, `R`, `Q`, `S`, and `M` across behavioral states.
#' @inheritSection BRQS State Variables
#' @inheritParams ramp.xds::get_MY_vars
#' @return a named [list] of state variables and module-specific derived quantities
#' @seealso [BRQS], [setup_MY_obj.BRQS]
#' @keywords internal
#' @export
#' @note Variable indices are stored at `xds_obj$MY_obj[[s]]$ix`.
get_MY_vars.BRQS<- function(y, xds_obj, s){
  with(xds_obj$MY_obj[[s]]$ix, {
    Su <- y[Su_ix]
    Bu <- y[Bu_ix]
    Ru <- y[Ru_ix]
    Qu <- y[Qu_ix]
    Sy <- y[Sy_ix]
    By <- y[By_ix]
    Ry <- y[Ry_ix]
    Qy <- y[Qy_ix]
    Sz <- y[Sz_ix]
    Bz <- y[Bz_ix]
    Rz <- y[Rz_ix]
    Qz <- y[Qz_ix]
    B <- Bu+By+Bz
    R <- Ru+Ry+Rz
    Q <- Qu+Qy+Qz
    S <- Su+Sy+Sz
    return(list(
      Su=Su, Bu=Bu, Ru=Ru, Qu=Qu,
      Sy=Sy, By=By, Ry=Ry, Qy=Qy,
      Sz=Sz, Bz=Bz, Rz=Rz, Qz=Qz,
      B=B, R=R, Q=Q, S=S, M=B+R+Q+S
    ))
  })
}

#' @title Mosquito bionomics for `BRQS` (**MY**)
#' @description Leave directly stored parameters unchanged. This module
#' does not configure bionomic parameter ports.
#' @inheritSection BRQS Parameters
#' @inheritSection BRQS Survival and Dispersal
#' @inheritSection BRQS Blood Feeding
#' @inheritParams ramp.xds::MBionomics
#' @seealso [BRQS], [MEffectSizes.BRQS]
#' @return an **`xds`** model object
#' @keywords internal
#' @export
MBionomics.BRQS <- function(t, y, xds_obj, s) {
  return(xds_obj)
}


#' @title Apply effect sizes for `BRQS` (**MY**)
#' @description Leave directly stored parameters unchanged. This module
#' does not configure bionomic parameter ports or apply effect sizes.
#' @inheritSection BRQS Parameters
#' @inheritSection BRQS Survival and Dispersal
#' @inheritParams ramp.xds::MEffectSizes
#' @seealso [BRQS], [MBionomics.BRQS]
#' @return an **`xds`** model object
#' @keywords internal
#' @export
MEffectSizes.BRQS <- function(t, y, xds_obj, s) {
  return(xds_obj)
}

#' @title Compute infective mosquito biting density for `BRQS` (**MY**)
#' @description Return the daily infective mosquito biting density on humans,
#' `fqZ`, for each patch in the `BRQS` module.
#' Use only the blood-feeding behavioral states: `fqZ = f*q*Bz`.
#' @inheritSection BRQS Infective Biting
#' @inheritParams ramp.xds::F_fqZ
#' @seealso [BRQS], [F_fqM.BRQS]
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
#' @note The current implementation obtains `q` using `get_f`
#' rather than `get_q`; this discrepancy requires repair.
F_fqZ.BRQS <- function(t, y, xds_obj, s) {
  f = get_f(xds_obj, s)
  q = get_f(xds_obj, s)
  fqZ = f*q*y[xds_obj$MY_obj[[s]]$ix$Bz_ix]
  return(fqZ)
}

#' @title Compute mosquito biting density for `BRQS` (**MY**)
#' @description Return the daily mosquito biting density on humans,
#' `fqM`, for each patch in the `BRQS` module.
#' Use only the blood-feeding behavioral states: `fqM = f*q*(Bu+By+Bz)`.
#' @inheritSection BRQS Infective Biting
#' @inheritParams ramp.xds::F_fqM
#' @seealso [BRQS], [ramp.xds::F_hbr], [F_fqZ.BRQS]
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
#' @note The current implementation obtains `q` using `get_f`
#' rather than `get_q`; this discrepancy requires repair.
F_fqM.BRQS <- function(t, y, xds_obj, s) {
  f = get_f(xds_obj, s)
  q = get_f(xds_obj, s)
  M = with(xds_obj$MY_obj[[s]]$ix, y[Bu_ix] + y[By_ix] + y[Bz_ix])
  fqM = f*q*M
  return(fqM)
}

#' @title Compute egg production for `BRQS` (**MY**)
#' @description Return the number of eggs laid by the adult mosquito population, per patch
#' @inheritSection BRQS Egg Laying
#' @inheritParams ramp.xds::F_eggs
#' @seealso [BRQS], [ramp.xds::EggLaying]
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
F_eggs.BRQS <- function(t, y, xds_obj, s) {

  G <- with(xds_obj$MY_obj[[s]]$ix, y[Qu_ix] + y[Qy_ix] + y[Qz_ix])
  eggs = with(xds_obj$MY_obj[[s]], {
    return(G*nu*eggsPerBatch)
  })
  return(eggs)
}


#' @title Get the feeding rate for `BRQS` (**MY**)
#' @description Return the stored realized feeding rate per patch.
#' @inheritSection BRQS Parameters
#' @inheritParams ramp.xds::get_f
#' @seealso [BRQS]
#' @details Read `f` without evaluating ports or applying effect sizes.
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
get_f.BRQS = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], f)
}

#' @title Get the human blood fraction for `BRQS` (**MY**)
#' @description Return the stored realized human blood fraction per patch.
#' @inheritSection BRQS Parameters
#' @inheritParams ramp.xds::get_q
#' @seealso [BRQS]
#' @details Read `q` without evaluating ports or applying effect sizes.
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
get_q.BRQS = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], q)
}

#' @title Get the adult mosquito mortality rate for `BRQS` (**MY**)
#' @description Return the stored realized adult mosquito mortality rate per patch.
#' @inheritSection BRQS Parameters
#' @inheritParams ramp.xds::get_g
#' @seealso [BRQS]
#' @details Read `g` without evaluating ports or applying effect sizes.
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
get_g.BRQS = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], g)
}

#' @title Get the patch emigration rate for `BRQS` (**MY**)
#' @description Return the stored realized patch emigration rate per patch.
#' @inheritSection BRQS Parameters
#' @inheritParams ramp.xds::get_sigma
#' @seealso [BRQS]
#' @details Read `sigma_b` without evaluating ports or applying effect sizes.
#' This is the blood-feeding emigration rate, not the egg-laying emigration rate.
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
get_sigma.BRQS = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], sigma_b)
}


#' @title Make a module object for `BRQS` (**MY**)
#' @description Create a `BRQS` module object and initialize its parameters.
#' Inspect parameters with [get_MY_pars.BRQS] and change supported values
#' with [change_MY_pars.BRQS]. Values in `options` override parameter arguments.
#' @seealso [get_MY_pars.BRQS], [change_MY_pars.BRQS]
#' @param nPatches is the number of patches, an integer
#' @param options a [list] of values that overrides the defaults
#' @param eip extrinsic incubation period
#' @param g mosquito mortality rate
#' @param sigma_b emigration rate while blood feeding
#' @param sigma_q emigration rate while egg laying
#' @param mu emigration loss
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu the oviposition rate: the number of egg batches laid in each patch,
#' per mosquito, per day
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

    MY_obj$f          <- checkIt(f, nPatches)

    MY_obj$nu         <- checkIt(nu, nPatches)

    MY_obj$theta        <- checkIt(theta, nPatches)

    MY_obj$g_b         <- checkIt(g_b, nPatches)

    MY_obj$g_q          <- checkIt(g_q, nPatches)

    MY_obj$g_s          <- checkIt(g_s, nPatches)

    MY_obj$g_r          <- checkIt(g_r, nPatches)

    MY_obj$sigma_b    <- checkIt(sigma_b, nPatches)

    MY_obj$sigma_q    <- checkIt(sigma_q, nPatches)

    MY_obj$sigma_s    <- checkIt(sigma_q, nPatches)

    MY_obj$mu_b    <- checkIt(mu_b, nPatches)

    MY_obj$mu_q    <- checkIt(mu_q, nPatches)

    MY_obj$mu_s    <- checkIt(mu_q, nPatches)

    MY_obj$xi           <- checkIt(xi, nPatches)

    MY_obj$eta           <- checkIt(eta, nPatches)

    MY_obj$q          <- checkIt(q, nPatches)

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


#' @title Set up the module object for `BRQS` (**MY**)
#' @description
#'
#' Set up a `BRQS` module object for the **MY** component and the \eqn{s^{th}} vector species.
#'
#' + Set the module name to `BRQS`.
#' + Create the module object with [make_MY_obj_BRQS], passing `nPatches` and `options`.
#' + Initialize the circadian function using [ramp.xds::setup_F_circadian.setup].
#' + Set the skill set using [ramp.xds::get_skillset_MY].
#' + Initialize the mosquito dispersal matrices using [ramp.xds::setup_K_matrix.zero],
#'   with `which_K = "Kb"`, `which_K = "Kq"`, and `which_K = "Ks"`.
#' + Set initial values using [ramp.xds::setup_MY_inits], passing `options` to override defaults.
#' + Return the updated **`xds`** model object.
#'
#' @inheritSection BRQS State Variables
#' @inheritSection BRQS Parameters
#' @inheritSection BRQS Dynamics
#' @seealso [BRQS], [make_MY_obj_BRQS], [setup_MY_inits.BRQS]
#'
#' @note The module object is stored at `xds_obj$MY_obj[[s]]`
#'
#' @inheritParams ramp.xds::setup_MY_obj
#' @return an **`xds`** model object
#' @keywords internal
#' @export
setup_MY_obj.BRQS = function(MYname, xds_obj, options=list(), s=1){
  xds_obj$MYname <- "BRQS"
  xds_obj$MY_obj[[s]] = make_MY_obj_BRQS(xds_obj$nPatches, options)
  xds_obj <- setup_F_circadian("setup", xds_obj, s=s)
  xds_obj <- setup_K_matrix("zero", xds_obj, options=list(which_K = "Kb"), s=s)
  xds_obj <- setup_K_matrix("zero", xds_obj, options=list(which_K = "Kq"), s=s)
  xds_obj <- setup_K_matrix("zero", xds_obj, options=list(which_K = "Ks"), s=s)
  xds_obj <- setup_MY_inits(xds_obj, options=options, s=s)
  return(xds_obj)
}

#' @title Get parameters for `BRQS` (**MY**)
#' @description Return the parameters for the `BRQS` module as a named list.
#' @inheritSection BRQS Parameters
#' @inheritParams ramp.xds::get_MY_pars
#' @return a named [list] of parameter values
#' @seealso [BRQS], [change_MY_pars.BRQS]
#' @seealso [make_MY_obj_BRQS]
#' @keywords internal
#' @export
get_MY_pars.BRQS <- function(xds_obj, s=1) {
  with(xds_obj$MY_obj[[s]], list(
    f=f, q=q, g=g, sigma_b=sigma_b, sigma_q = sigma_q, eip=eip, mu=mu,
    nu=nu, eggsPerBatch=eggsPerBatch, calK=calK
  ))
}

#' @title Change parameters for `BRQS` (**MY**)
#' @seealso [make_MY_obj_BRQS], [get_MY_pars.BRQS]
#' @description Change parameter values named in `options` for the `BRQS`
#' module, leaving other values unchanged.
#' Update the direct parameter fields for `f`, `q`, `nu`, `g_s`, `g_b`, `g_r`, `g_q`,
#' `mu_s`, `mu_b`, `mu_q`, `sigma_b`, `sigma_q`, `sigma_s`, `theta`, `xi`, and `eta`
#' when supplied, as well as `eip` and `eggsPerBatch`.
#' Call [ramp.xds::MBionomics] at time zero after the updates,
#' without calling [ramp.xds::MEffectSizes].
#' @inheritSection BRQS Parameters
#' @inheritParams ramp.xds::change_MY_pars
#' @return an **`xds`** model object
#' @keywords internal
#' @export
change_MY_pars.BRQS <- function(xds_obj, options=list(), s=1) {
  parameters <- c(
    "f", "q", "nu", "g_s", "g_b", "g_r", "g_q",
    "mu_s", "mu_b", "mu_q", "sigma_b", "sigma_q", "sigma_s",
    "theta", "xi", "eta"
  )
  for (par in intersect(parameters, names(options))) {
    xds_obj$MY_obj[[s]][[par]] <- checkIt(options[[par]], xds_obj$nPatches)
  }
  for (par in intersect(c("eip", "eggsPerBatch"), names(options))) {
    xds_obj$MY_obj[[s]][[par]] <- options[[par]]
  }
  y <- get_inits(xds_obj, flatten=TRUE)
  xds_obj <- MBionomics(0, y, xds_obj, s)
  return(xds_obj)
}

#' @title Make initial values for `BRQS` (**MY**)
#' @description Create a named list of initial values for the `BRQS` module.
#' Values in `options` override the defaults.
#' @details Initial values must be numeric scalars or have length `nPatches`.
#' Scalars are expanded to that length.
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
#' @return a named [list] of initial values
#' @keywords internal
#' @export
#' @seealso [change_MY_inits.BRQS]
make_MY_inits_BRQS = function(nPatches, options = list(),
                              Su=100, Bu=100, Ru=1, Qu=10,
                              Sy=0, By=0, Ry=0, Qy=0,
                              Sz=0, Bz=0, Rz=0, Qz=0){

  with(options,{
    Su = checkIt(Su, nPatches)
    Bu = checkIt(Bu, nPatches)
    Ru = checkIt(Ru, nPatches)
    Qu = checkIt(Qu, nPatches)
    Sy = checkIt(Sy, nPatches)
    By = checkIt(By, nPatches)
    Ry = checkIt(Ry, nPatches)
    Qy = checkIt(Qy, nPatches)
    Sz = checkIt(Sz, nPatches)
    Bz = checkIt(Bz, nPatches)
    Rz = checkIt(Rz, nPatches)
    Qz = checkIt(Qz, nPatches)

    return(list(Su=Su, Bu=Bu, Ru=Ru, Qu=Qu,
                Sy=Sy, By=By, Ry=Ry, Qy=Qy,
                Sz=Sz, Bz=Bz, Rz=Rz, Qz=Qz))
  })}

#' @title Set up initial values for `BRQS` (**MY**)
#' @description Set up initial values for the `BRQS` module using
#' [make_MY_inits_BRQS]. Values in `options` override the defaults.
#' @inheritSection BRQS State Variables
#' @inheritParams ramp.xds::setup_MY_inits
#' @return an **`xds`** model object
#' @keywords internal
#' @export
#' @seealso [change_MY_inits.BRQS]
#' @note Initial values are stored at `xds_obj$MY_obj[[s]]$inits`.
#' @seealso [ramp.xds::get_MY_inits]
setup_MY_inits.BRQS = function(xds_obj, options=list(), s=1){
  xds_obj$MY_obj[[s]]$inits = make_MY_inits_BRQS(xds_obj$nPatches, options)
  return(xds_obj)
}

#' @title Change initial values for `BRQS` (**MY**)
#' @description Change initial values named in `options` for the `BRQS`
#' module, leaving other values unchanged.
#' @inheritSection BRQS State Variables
#' @seealso [BRQS], [setup_MY_inits.BRQS]
#' @inheritParams ramp.xds::change_MY_inits
#' @return an **`xds`** model object
#' @keywords internal
#' @export
#' @seealso [make_MY_inits_BRQS]
#' @note Stored initial values for this module are at
#' `xds_obj$MY_obj[[s]]$inits`.
change_MY_inits.BRQS <- function(xds_obj, options = list(), s=1) {
  inits = with(get_MY_inits(xds_obj, s), with(options,
    list(Su=Su, Bu=Bu, Ru=Ru, Qu=Qu,
         Sy=Sy, By=By, Ry=Ry, Qy=Qy,
         Sz=Sz, Bz=Bz, Rz=Rz, Qz=Qz)))

  xds_obj$MY_obj[[s]]$inits=inits
  return(xds_obj)
}


#' @title Set up variable indices for `BRQS` (**MY**)
#' @description Assign consecutive variable indices for the `BRQS` module,
#' starting after `xds_obj$max_ix`, and update the maximum index.
#' Allocate blocks in the order `Su_ix`, `Bu_ix`, `Ru_ix`, `Qu_ix`, `Sy_ix`, `By_ix`, `Ry_ix`, `Qy_ix`, `Sz_ix`, `Bz_ix`, `Rz_ix`, `Qz_ix`.
#' Each block has length `nPatches`.
#' @inheritSection BRQS State Variables
#' @inheritParams ramp.xds::setup_MY_ix
#' @return an **`xds`** model object
#' @importFrom utils tail
#' @keywords internal
#' @export
#' @note Indices are stored at `xds_obj$MY_obj[[s]]$ix`.
#' The maximum index is stored at `xds_obj$max_ix`.
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
  xds_obj$MY_obj[[s]]$ix = list(
    Su_ix = Su_ix, Bu_ix = Bu_ix, Ru_ix=Ru_ix, Qu_ix = Qu_ix,
    Sy_ix = Sy_ix, By_ix = By_ix, Ry_ix=Ry_ix, Qy_ix = Qy_ix,
    Sz_ix = Sz_ix, Bz_ix = Bz_ix, Rz_ix=Rz_ix, Qz_ix = Qz_ix)


  return(xds_obj)
})}


#' @title Parse outputs for `BRQS` (**MY**)
#' @description Extract the model state outputs for the `BRQS` module
#' from `outputs` and return them by name in a list.
#' Sum behavioral states to compute total density `M`, infected density `Y`
#' (exposed and infectious), and infectious density `Z`.
#' Also compute `y = Y/M` and `z = Z/M`.
#' @note Variable indices are stored at `xds_obj$MY_obj[[s]]$ix`.
#' @inheritSection BRQS State Variables
#' @inheritSection BRQS Infective Biting
#' @inheritParams ramp.xds::parse_MY_orbits
#' @seealso [BRQS], [get_MY_vars.BRQS]
#' @return a [list]
#' @keywords internal
#' @export
parse_MY_orbits.BRQS <- function(outputs, xds_obj, s) {with(xds_obj$MY_obj[[s]]$ix,{
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

#' @title Get the skill set for `BRQS` (**MY**)
#'
#' @description Return the capabilities and compatibilities of the `BRQS` module as a list.
#'
#' @inheritParams ramp.xds::get_skillset_MY
#'
#' @return the **MY** module skill set, as a list
#'
#' @keywords internal
#' @export
get_skillset_MY.BRQS = function(xds_obj, s){
  skills = list(
    f_is_port = FALSE,
    q_is_port = FALSE,
    g_is_port = FALSE,
    sigma_is_port = FALSE,
    sigma_b_is_port = FALSE,
    sigma_q_is_port = FALSE,
    mu_is_port = FALSE,
    nu_is_port = FALSE,
    eip_is_port = FALSE,
    K_matrix_is_port = FALSE,
    dts = FALSE,
    delay = FALSE
  )
  return(skills)
}
