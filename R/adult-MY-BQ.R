# specialized methods for the adult mosquito BQ model
# TODO: Repair F_fqZ.BQ and F_fqM.BQ: obtain q using get_q, not get_f.
# TODO: Review change_MY_inits.BQ repair: read fallback values from inits.
# Confirm partial updates preserve unspecified values and other species.
# TODO: Reconcile emigration parameters: get_MY_pars.BQ reads sigma_t, but
# make_MY_obj_BQ initializes separate sigma_b and sigma_q parameter objects.
# Align the getter with the constructor and change_MY_pars.BQ, and test coverage.

#' @title The `BQ` module for the **MY** component
#' @description
#' Implements the **MY** component using a BQ (Blood-feeding / egg-laying Queue)
#' model of adult mosquito ecology and infection dynamics. Mosquitoes alternate
#' between blood-feeding (**B**) and gravid / egg-laying (**Q**) states.
#' Infection is tracked as uninfected (**u**), exposed (**y**), or infectious
#' (**z**).
#'
#' @section State Variables:
#' \describe{
#'   \item{`Bu`}{density of uninfected, blood-feeding mosquitoes}
#'   \item{`Qu`}{density of uninfected, gravid (egg-laying) mosquitoes}
#'   \item{`By`}{density of infected (exposed), blood-feeding mosquitoes}
#'   \item{`Qy`}{density of infected (exposed), gravid mosquitoes}
#'   \item{`Bz`}{density of infectious, blood-feeding mosquitoes}
#'   \item{`Qz`}{density of infectious, gravid mosquitoes}
#'   \item{`M`}{total adult mosquito density, derived during parsing}
#'   \item{`Y`}{total infected and infectious mosquito density, derived during parsing}
#'   \item{`Z`}{total infectious mosquito density, derived during parsing}
#' }
#'
#' @section Dynamical Terms:
#' \describe{
#'   \item{`Lambda` or \eqn{\Lambda}}{adult mosquito emergence rate, computed by [ramp.xds::Emergence]}
#'   \item{`kappa` or \eqn{\kappa}}{net infectiousness, computed by [ramp.xds::Transmission]}
#' }
#'
#' @section Parameters:
#' \describe{
#'   \item{`f`}{blood feeding rate --- use [ramp.xds::setup_f] or [ramp.xds::setup_blood_feeding_rates]}
#'   \item{`q`}{human fraction --- use [ramp.xds::setup_q] or [ramp.xds::setup_human_fractions]}
#'   \item{`g`}{mosquito mortality rate --- use [ramp.xds::setup_g] or [ramp.xds::setup_mosquito_mortality]}
#'   \item{`sigma_b`}{emigration rate during blood feeding}
#'   \item{`sigma_q`}{emigration rate during egg laying}
#'   \item{`mu`}{fraction of emigrating mosquitoes lost}
#'   \item{`eip` or \eqn{\tau}}{extrinsic incubation period; \eqn{\phi=1/\tau}}
#'   \item{`nu` or \eqn{\nu}}{oviposition rate: the number of egg batches laid in each patch, per mosquito, per day --- use [ramp.xds::setup_nu] or [ramp.xds::setup_egg_batch_rates]}
#'   \item{`eggsPerBatch` or \eqn{\Gamma}}{number of eggs per batch}
#' }
#'
#' After basic setup, [change_MY_pars.BQ] can be used to change parameter
#' values. Bionomic parameters are configured as static baselines by default.
#' @section Survival and Dispersal:
#' Blood-feeding and gravid mosquitoes have separate demographic matrices,
#' \eqn{\Omega_b} and \eqn{\Omega_q}. Each combines mortality \eqn{g}, the
#' corresponding emigration rate (\eqn{\sigma_b} or \eqn{\sigma_q}), the
#' emigration-loss fraction \eqn{\mu}, and its dispersal matrix (\eqn{K_b} or
#' \eqn{K_q}). The columns of each dispersal matrix sum to zero; its
#' off-diagonal elements describe the fractions of surviving emigrants that
#' arrive in other patches.
#'
#' @section Blood Feeding:
#' Human blood feeding occurs from the blood-feeding state at rate \eqn{f}.
#' A blood meal infects a susceptible mosquito with probability \eqn{\kappa};
#' consequently, \eqn{fq\kappa} is the force of infection on susceptible
#' blood-feeding mosquitoes.
#'
#' @section Dynamics:
#' Let \eqn{\phi=1/\tau}, where \eqn{\tau} is the extrinsic incubation period.
#' Mosquitoes move between blood-feeding and gravid states at rate \eqn{f} and
#' oviposit at rate \eqn{\nu}. Infection is acquired by uninfected
#' blood-feeding mosquitoes; exposed mosquitoes progress to the infectious
#' class at rate \eqn{\phi}. The compartment equations are:
#' \deqn{
#' \begin{array}{rl}
#' dB_u/dt &= \Lambda+\nu Q_u-fB_u-\Omega_b\cdot B_u \\\\
#' dQ_u/dt &= f(1-q\kappa)B_u-\nu Q_u-\Omega_q\cdot Q_u \\\\
#' dB_y/dt &= \nu Q_y-(f+\phi)B_y-\Omega_b\cdot B_y \\\\
#' dQ_y/dt &= fq\kappa B_u+fB_y-(\nu+\phi)Q_y-\Omega_q\cdot Q_y \\\\
#' dB_z/dt &= \phi B_y+\nu Q_z-fB_z-\Omega_b\cdot B_z \\\\
#' dQ_z/dt &= \phi Q_y+fB_z-\nu Q_z-\Omega_q\cdot Q_z
#' \end{array}
#' }
#'
#' @section Infective Biting:
#' The patch-level human biting density, `fqM`, is based on the total
#' blood-feeding mosquito density \eqn{B_u+B_y+B_z}. The infective biting
#' density, `fqZ`, is based on infectious blood-feeding mosquitoes \eqn{B_z}.
#' [F_fqM.BQ] and [F_fqZ.BQ] compute these quantities. [ramp.xds::compute_EIR]
#' combines `fqZ` with the biting distribution matrix and resident-bite
#' fraction to compute stratum-level `EIR`; [ramp.xds::F_hbr] uses `fqM` to
#' compute stratum-level `HBR`.
#'
#' @section Egg Laying:
#' Gravid mosquitoes lay eggs at rate \eqn{\nu}; the oviposition rate and
#' number of eggs per batch determine patch-level egg production:
#' \deqn{G=\nu(Q_u+Q_y+Q_z)\Gamma.}
#' [F_eggs.BQ] computes \eqn{G}, which is passed to [ramp.xds::EggLaying] for
#' distribution among aquatic habitats.
#'
#' @seealso
#' + Basic setup: [make_MY_obj_BQ], [setup_MY_obj.BQ]
#' + Dynamics: [dMYdt.BQ]
#' + Biting and eggs: [F_fqZ.BQ], [F_fqM.BQ], [F_eggs.BQ]
#' + Parameters: [get_MY_pars.BQ], [change_MY_pars.BQ]
#' + Variables: [make_MY_inits_BQ], [setup_MY_inits.BQ], [parse_MY_orbits.BQ]
#'
#' @name BQ
#' @rdname BQ
NULL

#' @title Compute the derivatives for `BQ` (**MY**)
#' @description Compute the derivatives for the `BQ` module:
#' \deqn{\begin{array}{rl}
#' dB_u/dt &= \Lambda+\nu Q_u-fB_u-\Omega_b\cdot B_u \\
#' dQ_u/dt &= f(1-q\kappa)B_u-\nu Q_u-\Omega_q\cdot Q_u \\
#' dB_y/dt &= \nu Q_y-(f+\phi)B_y-\Omega_b\cdot B_y \\
#' dQ_y/dt &= fq\kappa B_u+fB_y-(\nu+\phi)Q_y-\Omega_q\cdot Q_y \\
#' dB_z/dt &= \phi B_y+\nu Q_z-fB_z-\Omega_b\cdot B_z \\
#' dQ_z/dt &= \phi Q_y+fB_z-\nu Q_z-\Omega_q\cdot Q_z
#' \end{array}}
#' @inheritSection BQ State Variables
#' @inheritSection BQ Dynamical Terms
#' @inheritSection BQ Parameters
#' @inheritSection BQ Dynamics
#' @inheritParams ramp.xds::dMYdt
#' @seealso [BQ], [F_fqZ.BQ], [F_eggs.BQ]
#' @note Derivatives are returned in blocks for `Bu`, `Qu`, `By`, `Qy`,
#' `Bz`, and `Qz`, each in patch order.
#' @return a [numeric] vector of derivatives in the module's state-variable order
#' @keywords internal
#' @export
dMYdt.BQ <- function(t, y, xds_obj, s){

  Lambda = xds_obj$terms$Lambda[[s]]
  kappa = xds_obj$terms$kappa[[s]]

  with(get_MY_vars(y, xds_obj, s),{
    with(xds_obj$MY_obj[[s]],{

      dBu <- Lambda + nu*Qu - f*Bu - Omega_b %*% Bu
      dQu <- f*(1-q*kappa)*Bu - nu*Qu - Omega_q %*% Qu

      dBy <- nu*Qy - (f+phi)*By - Omega_b %*% By
      dQy <- f*q*kappa*Bu + f*By - (nu+phi)*Qy - Omega_q %*% Qy

      dBz <- phi*By + nu*Qz - f*Bz - Omega_b %*% Bz
      dQz <- phi*Qy + f*Bz - nu*Qz - Omega_q %*% Qz

      return(c(dBu, dQu, dBy, dQy, dBz, dQz))
    })
  })
}

#' @title Get the state variables for `BQ` (**MY**)
#' @description Extract the state variables for the `BQ` module
#' from `y` and return them by name in a list for use during model solving.
#' @inheritSection BQ State Variables
#' @inheritParams ramp.xds::get_MY_vars
#' @return a named [list] of state variables and module-specific derived quantities
#' @keywords internal
#' @export
#' @note Variable indices are stored at `xds_obj$MY_obj[[s]]$ix`.
get_MY_vars.BQ <- function(y, xds_obj, s){
  with(xds_obj$MY_obj[[s]]$ix,
       return(list(
         Bu = y[Bu_ix],
         Qu = y[Qu_ix],
         By = y[By_ix],
         Qy = y[Qy_ix],
         Bz = y[Bz_ix],
         Qz = y[Qz_ix]
)))}

#' @title Mosquito bionomics for `BQ` (**MY**)
#' @description Implements [MBionomics] for the BQ model.
#' @inheritSection BQ Parameters
#' @inheritSection BQ Survival and Dispersal
#' @inheritParams ramp.xds::MBionomics
#' @seealso [BQ], [MEffectSizes.BQ]
#' @return an **`xds`** model object
#' @keywords internal
#' @export
MBionomics.BQ <- function(t, y, xds_obj, s) {
  xds_obj$MY_obj[[s]]$f_t <- F_f(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$q_t <- F_q(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$g_t <- F_g(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$sigma_b_t <- F_sigma_b(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$sigma_q_t <- F_sigma_q(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$mu_t <- F_mu(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$nu_t <- F_nu(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$eip_t <- F_eip(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$es_f <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_q <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_g <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_sigma_b <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_sigma_q <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_mu <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_nu <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_eip <- rep(1, xds_obj$nPatches)
  return(xds_obj)
}


#' @title Apply effect sizes for `BQ` (**MY**)
#' @description Implements [MEffectSizes] for the BQ model.
#' @inheritSection BQ Parameters
#' @inheritSection BQ Survival and Dispersal
#' @inheritParams ramp.xds::MEffectSizes
#' @seealso [BQ], [MBionomics.BQ]
#' @return an **`xds`** model object
#' @keywords internal
#' @export
MEffectSizes.BQ <- function(t, y, xds_obj, s) {
  with(xds_obj$MY_obj[[s]], {
    xds_obj$MY_obj[[s]]$f <- f_t * es_f
    xds_obj$MY_obj[[s]]$q <- q_t * es_q
    xds_obj$MY_obj[[s]]$g <- g_t * es_g
    xds_obj$MY_obj[[s]]$sigma_b <- sigma_b_t * es_sigma_b
    xds_obj$MY_obj[[s]]$sigma_q <- sigma_q_t * es_sigma_q
    xds_obj$MY_obj[[s]]$mu <- mu_t * es_mu
    xds_obj$MY_obj[[s]]$nu <- nu_t * es_nu
    xds_obj$MY_obj[[s]]$eip <- eip_t * es_eip
    xds_obj$MY_obj[[s]]$phi <- 1 / xds_obj$MY_obj[[s]]$eip
    xds_obj <- change_Omega_b(xds_obj, s)
    xds_obj <- change_Omega_q(xds_obj, s)
    return(xds_obj)
  })
}

#' @title Compute infective mosquito biting density for `BQ` (**MY**)
#' @description Return the daily infective mosquito biting density on humans,
#' `fqZ`, for each patch in the `BQ` module.
#' Use only the blood-feeding behavioral states: `fqZ = f*q*Bz`.
#' @inheritSection BQ Infective Biting
#' @inheritParams ramp.xds::F_fqZ
#' @seealso [BQ], [F_fqM.BQ]
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
#' @note The current implementation obtains `q` using `get_f`
#' rather than `get_q`; this discrepancy requires repair.
F_fqZ.BQ <- function(t, y, xds_obj, s) {
  f = get_f(xds_obj, s)
  q = get_f(xds_obj, s)
  fqZ = f*q*y[xds_obj$MY_obj[[s]]$ix$Bz_ix]
  return(fqZ)
}

#' @title Compute mosquito biting density for `BQ` (**MY**)
#' @description Return the daily mosquito biting density on humans,
#' `fqM`, for each patch in the `BQ` module.
#' Use only the blood-feeding behavioral states: `fqM = f*q*(Bu+By+Bz)`.
#' @inheritSection BQ Infective Biting
#' @inheritParams ramp.xds::F_fqM
#' @seealso [BQ], [ramp.xds::F_hbr], [F_fqZ.BQ]
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
#' @note The current implementation obtains `q` using `get_f`
#' rather than `get_q`; this discrepancy requires repair.
F_fqM.BQ <- function(t, y, xds_obj, s) {
  f = get_f(xds_obj, s)
  q = get_f(xds_obj, s)
  M = with(xds_obj$MY_obj[[s]]$ix, y[Bu_ix] + y[By_ix] + y[Bz_ix])
  fqM = f*q*M
  return(fqM)
}

#' @title Compute egg production for `BQ` (**MY**)
#' @description Return the number of eggs laid by the adult mosquito population, per patch.
#' @inheritSection BQ Egg Laying
#' @inheritParams ramp.xds::F_eggs
#' @seealso [BQ], [ramp.xds::EggLaying]
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
F_eggs.BQ <- function(t, y, xds_obj, s) {

  G <- with(get_MY_vars(y, xds_obj, s), Qu + Qy + Qz)
  eggs = with(xds_obj$MY_obj[[s]], {
    return(G*nu*eggsPerBatch)
  })
  return(eggs)
}

#' @title Get the feeding rate for `BQ` (**MY**)
#' @description Return the stored realized feeding rate per patch.
#' @inheritSection BQ Parameters
#' @inheritParams ramp.xds::get_f
#' @seealso [BQ]
#' @details Read `f` without evaluating ports or applying effect sizes.
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
get_f.BQ = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], f)
}

#' @title Get the human blood fraction for `BQ` (**MY**)
#' @description Return the stored realized human blood fraction per patch.
#' @inheritSection BQ Parameters
#' @inheritParams ramp.xds::get_q
#' @seealso [BQ]
#' @details Read `q` without evaluating ports or applying effect sizes.
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
get_q.BQ = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], q)
}

#' @title Get the adult mosquito mortality rate for `BQ` (**MY**)
#' @description Return the stored realized adult mosquito mortality rate per patch.
#' @inheritSection BQ Parameters
#' @inheritParams ramp.xds::get_g
#' @seealso [BQ]
#' @details Read `g` without evaluating ports or applying effect sizes.
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
get_g.BQ = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], g)
}

#' @title Get the patch emigration rate for `BQ` (**MY**)
#' @description Return the stored realized patch emigration rate per patch.
#' @inheritSection BQ Parameters
#' @inheritParams ramp.xds::get_sigma
#' @seealso [BQ]
#' @details Read `sigma_b` without evaluating ports or applying effect sizes.
#' This is the blood-feeding emigration rate, not the egg-laying emigration rate.
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
get_sigma.BQ = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], sigma_b)
}


#' @title Make a module object for `BQ` (**MY**)
#' @description Create a `BQ` module object and initialize its parameters.
#' Inspect parameters with [get_MY_pars.BQ] and change supported values
#' with [change_MY_pars.BQ]. Values in `options` override parameter arguments.
#' @seealso [get_MY_pars.BQ], [change_MY_pars.BQ]
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
#' @seealso [BQ], [setup_MY_obj.BQ]
#' @importFrom ramp.xds checkIt
#' @keywords internal
#' @export
make_MY_obj_BQ = function(nPatches, options=list(), eip=12,
                          g=1/12, sigma_b=1/8, sigma_q=1/8, mu=0, f=0.5, q=0.95,
                          nu=1, eggsPerBatch=60){

  with(options,{
    MY_obj <- list()
    class(MY_obj) <- "BQ"

    MY_obj$nPatches <- nPatches

    eip_par <- list()
    class(eip_par) <- 'static'
    MY_obj$eip_par <- eip_par
    MY_obj$eip     <- eip

    MY_obj <- setup_eip_obj(checkIt(eip, nPatches), MY_obj)
    MY_obj <- setup_f_obj(checkIt(f, nPatches), MY_obj)
    MY_obj <- setup_q_obj(checkIt(q, nPatches), MY_obj)
    MY_obj <- setup_g_obj(checkIt(g, nPatches), MY_obj)
    MY_obj <- setup_mu_obj(checkIt(mu, nPatches), MY_obj)
    MY_obj <- setup_nu_obj(checkIt(nu, nPatches), MY_obj)
    MY_obj <- setup_sigma_q_obj(checkIt(sigma_q, nPatches), MY_obj)
    MY_obj <- setup_sigma_b_obj(checkIt(sigma_b, nPatches), MY_obj)
    MY_obj$eggsPerBatch <- eggsPerBatch
    MY_obj$phi <- 1/MY_obj$eip

    return(MY_obj)
})}


#' @title Set up the module object for `BQ` (**MY**)
#' @description
#'
#' Set up a `BQ` module object for the **MY** component and the \eqn{s^{th}} vector species.
#'
#' + Set the module name to `BQ`.
#' + Create the module object with [make_MY_obj_BQ], passing `nPatches` and `options`.
#' + Initialize the circadian function using [ramp.xds::setup_F_circadian.setup].
#' + Set the skill set using [ramp.xds::get_skillset_MY].
#' + Initialize the mosquito dispersal matrices using [ramp.xds::setup_K_matrix.zero],
#'   with `which_K = "Kb"` and `which_K = "Kq"`.
#' + Set initial values using [ramp.xds::setup_MY_inits], passing `options` to override defaults.
#' + Return the updated **`xds`** model object.
#'
#' @inheritSection BQ State Variables
#' @inheritSection BQ Parameters
#' @seealso [BQ], [make_MY_obj_BQ], [setup_MY_inits.BQ]
#'
#' @note The module object is stored at `xds_obj$MY_obj[[s]]`.
#'
#' @inheritParams ramp.xds::setup_MY_obj
#' @return an **`xds`** model object
#' @keywords internal
#' @export
setup_MY_obj.BQ = function(MYname, xds_obj, options=list(), s=1){
  xds_obj$MYname <- "BQ"
  xds_obj$MY_obj[[s]] = make_MY_obj_BQ(xds_obj$nPatches, options)
  xds_obj <- setup_F_circadian("setup", xds_obj, s=s)
  xds_obj$MY_obj[[s]]$skill_set <- get_skillset_MY(xds_obj, s)
  xds_obj <- setup_K_matrix("zero", xds_obj, options=list(which_K = "Kb"), s=s)
  xds_obj <- setup_K_matrix("zero", xds_obj, options=list(which_K = "Kq"), s=s)
  xds_obj <- setup_MY_inits(xds_obj, options=options, s=s)
  return(xds_obj)
}


#' @title Get parameters for `BQ` (**MY**)
#'
#' @description Return the parameters for the `BQ` module as a named list.
#' @inheritSection BQ Parameters
#'
#' @inheritParams ramp.xds::get_MY_pars
#'
#' @return a named [list] of parameter values
#' @seealso [BQ], [change_MY_pars.BQ]
#' @seealso [make_MY_obj_BQ]
#' @keywords internal
#' @export
get_MY_pars.BQ <- function(xds_obj, s=1) {
  with(xds_obj$MY_obj[[s]], list(
    f=f_t, q=q_t, g=g_t, sigma=sigma_t, eip=eip, mu=mu_t,
    nu=nu_t, eggsPerBatch=eggsPerBatch, calK=calK
  ))
}


#' @title Change parameters for `BQ` (**MY**)
#' @seealso [make_MY_obj_BQ], [get_MY_pars.BQ]
#' @description Change parameter values named in `options` for the `BQ`
#' module, leaving other values unchanged.
#' Refresh baseline bionomics at time zero using [ramp.xds::MBionomics],
#' without calling [ramp.xds::MEffectSizes].
#' @inheritSection BQ Parameters
#' @inheritParams ramp.xds::change_MY_pars
#' @return an **`xds`** model object
#' @keywords internal
#' @export
change_MY_pars.BQ <- function(xds_obj, options=list(), s=1) {
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

    return(xds_obj)
  }))}

#' @title Make initial values for `BQ` (**MY**)
#' @description Create a named list of initial values for the `BQ` module.
#' Values in `options` override the defaults.
#' @details Initial values must be numeric scalars or have length `nPatches`.
#' Scalars are expanded to that length.
#' @param nPatches the number of patches
#' @param options a [list] of values that overrides the defaults
#' @param Bu total uninfected, not gravid mosquito density at each patch
#' @param Qu total uninfected, gravid mosquito density at each patch
#' @param By infected, not gravid mosquito density at each patch
#' @param Qy infected, gravid mosquito density at each patch
#' @param Bz infectious, not gravid mosquito density at each patch
#' @param Qz infectious, gravid mosquito density at each patch
#' @inheritSection BQ State Variables
#' @seealso [BQ], [setup_MY_inits.BQ]
#' @return a named [list] of initial values
#' @keywords internal
#' @export
#' @seealso [change_MY_inits.BQ]
make_MY_inits_BQ = function(nPatches, options = list(),
                             Bu=5, Qu=1, By=1, Qy=1, Bz=1, Qz=1){
  with(options,{
    Bu = checkIt(Bu, nPatches)
    Qu = checkIt(Qu, nPatches)
    By = checkIt(By, nPatches)
    Qy = checkIt(Qy, nPatches)
    Bz = checkIt(Bz, nPatches)
    Qz = checkIt(Qz, nPatches)

    return(list(Bu=Bu, Qu=Qu, By=By, Qy=Qy, Bz=Bz, Qz=Qz))
  })}

#' @title Set up initial values for `BQ` (**MY**)
#' @description Set up initial values for the `BQ` module using
#' [make_MY_inits_BQ]. Values in `options` override the defaults.
#' @inheritSection BQ State Variables
#' @inheritParams ramp.xds::setup_MY_inits
#' @return an **`xds`** model object
#' @keywords internal
#' @export
#' @seealso [change_MY_inits.BQ]
#' @note Initial values are stored at `xds_obj$MY_obj[[s]]$inits`.
#' @seealso [ramp.xds::get_MY_inits]
setup_MY_inits.BQ = function(xds_obj, options=list(), s=1){
  xds_obj$MY_obj[[s]]$inits = make_MY_inits_BQ(xds_obj$nPatches, options)
  return(xds_obj)
}


#' @title Change initial values for `BQ` (**MY**)
#' @description Change initial values named in `options` for the `BQ`
#' module, leaving other values unchanged.
#' @inheritParams ramp.xds::change_MY_inits
#' @return an **`xds`** model object
#' @keywords internal
#' @export
#' @seealso [make_MY_inits_BQ], [setup_MY_inits.BQ]
#' @note Stored initial values for this module are at
#' `xds_obj$MY_obj[[s]]$inits`.
change_MY_inits.BQ <- function(xds_obj, options=list(), s=1) {
  with(xds_obj$MY_obj[[s]]$inits,
    with(options,{
      xds_obj$MY_obj[[s]]$inits$Bu = Bu
      xds_obj$MY_obj[[s]]$inits$Qu = Qu
      xds_obj$MY_obj[[s]]$inits$By = By
      xds_obj$MY_obj[[s]]$inits$Qy = Qy
      xds_obj$MY_obj[[s]]$inits$Bz = Bz
      xds_obj$MY_obj[[s]]$inits$Qz = Qz
      return(xds_obj)
  }))}

#' @title Set up variable indices for `BQ` (**MY**)
#' @description Assign consecutive variable indices for the `BQ` module,
#' starting after `xds_obj$max_ix`, and update the maximum index.
#' Allocate blocks in the order `Bu_ix`, `Qu_ix`, `By_ix`, `Qy_ix`, `Bz_ix`, `Qz_ix`.
#' Each block has length `nPatches`.
#' @inheritSection BQ State Variables
#' @inheritParams ramp.xds::setup_MY_ix
#' @return an **`xds`** model object
#' @importFrom utils tail
#' @keywords internal
#' @export
#' @note Indices are stored at `xds_obj$MY_obj[[s]]$ix`.
#' The maximum index is stored at `xds_obj$max_ix`.
setup_MY_ix.BQ <- function(xds_obj, s) {with(xds_obj,{

  Bu_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Bu_ix, 1)

  Qu_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Qu_ix, 1)

  By_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(By_ix, 1)

  Qy_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Qy_ix, 1)

  Bz_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Bz_ix, 1)

  Qz_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Qz_ix, 1)

  xds_obj$max_ix = max_ix
  xds_obj$MY_obj[[s]]$ix = list(Bu_ix=Bu_ix, Qu_ix=Qu_ix,
                          By_ix=By_ix, Qy_ix=Qy_ix,
                          Bz_ix=Bz_ix, Qz_ix=Qz_ix)
  return(xds_obj)
})}

#' @title Parse outputs for `BQ` (**MY**)
#' @description Extract the model state outputs for the `BQ` module
#' from `outputs` and return them by name in a list.
#' Sum behavioral states to compute total density `M`, infected density `Y`
#' (exposed and infectious), and infectious density `Z`.
#' Also compute `y = Y/M` and `z = Z/M`.
#' @note Variable indices are stored at `xds_obj$MY_obj[[s]]$ix`.
#' @inheritSection BQ State Variables
#' @inheritSection BQ Infective Biting
#' @inheritParams ramp.xds::parse_MY_orbits
#' @seealso [BQ], [get_MY_vars.BQ]
#' @return a [list]
#' @keywords internal
#' @export
parse_MY_orbits.BQ <- function(outputs, xds_obj, s) {
  with(xds_obj$MY_obj[[s]]$ix,{
    Bu = outputs[,Bu_ix]
    Qu = outputs[,Qu_ix]
    By = outputs[,By_ix]
    Qy = outputs[,Qy_ix]
    Bz = outputs[,Bz_ix]
    Qz = outputs[,Qz_ix]
    M =  Bu+By+Bz+Qu+Qy+Qz
    Y =  By+Qy+Bz+Qz
    Z =  Bz+Qz
    y = Y/M
    z = Z/M
  return(list(Bu=Bu, Qu=Qu,
              By=By, Qy=Qy,
              Bz=Bz, Qz=Qz,
              M=M, Y=Y, Z=Z, y=y, z=z))
})}

#' @title Check `BQ` (**MY**)
#'
#' @inheritParams ramp.xds::check_MY
#'
#' @return an **`xds`** model object
#' @keywords internal
#' @export
check_MY.BQ = function(xds_obj, s){
  return(xds_obj)
}

#' @title Get the skill set for `BQ` (**MY**)
#'
#' @description Return the capabilities and compatibilities of the `BQ` module as a list.
#'
#' @inheritParams ramp.xds::get_skillset_MY
#'
#' @return the **MY** module skill set, as a list
#'
#' @keywords internal
#' @export
get_skillset_MY.BQ = function(xds_obj, s){
  skills = list(
    f_is_port = TRUE,
    q_is_port = TRUE,
    g_is_port = TRUE,
    sigma_is_port = FALSE,
    sigma_b_is_port = TRUE,
    sigma_q_is_port = TRUE,
    mu_is_port = TRUE,
    nu_is_port = TRUE,
    eip_is_port = TRUE,
    K_matrix_is_port = FALSE,
    dts = FALSE,
    delay = FALSE
  )
  return(skills)
}
