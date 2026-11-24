# specialized methods for the adult mosquito BQm model
# TODO: Repair F_fqZ.BQm and F_fqM.BQm: obtain q using get_q, not get_f.
# TODO: Review change_MY_inits.BQm repair: read fallback values from inits.
# Confirm partial updates preserve unspecified values and other species.
# TODO: Reconcile emigration parameters: get_MY_pars.BQm reads sigma_t, but
# make_MY_obj_BQm initializes separate sigma_b and sigma_q parameter objects.
# Align the getter with the constructor and change_MY_pars.BQm, and test coverage.

#' @title The `BQm` module for the **MY** component
#' @description
#' Implements a mosquito blood-feeding / egg-laying queue model. Its
#' initialized and parsed state vector tracks blood-feeding (**B**) and gravid
#' (**Q**) mosquitoes in uninfected (**u**), exposed (**y**), and infectious
#' (**z**) classes.
#'
#' @section State Variables:
#' \describe{
#'   \item{`Bu`, `By`, `Bz`}{uninfected, exposed, and infectious blood-feeding mosquitoes}
#'   \item{`Qu`, `Qy`, `Qz`}{uninfected, exposed, and infectious gravid mosquitoes}
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
#'   \item{`eip` or \eqn{\tau}}{extrinsic incubation period}
#'   \item{`nu` or \eqn{\nu}}{oviposition rate: the number of egg batches laid in each patch, per mosquito, per day --- use [ramp.xds::setup_nu] or [ramp.xds::setup_egg_batch_rates]}
#'   \item{`eggsPerBatch` or \eqn{\Gamma}}{number of eggs per batch}
#' }
#' After basic setup, [change_MY_pars.BQm] can be used to change parameter
#' values. Bionomic parameters are configured as static baselines by default.
#'
#' @section Survival and Dispersal:
#' Blood-feeding and gravid mosquitoes have separate demographic matrices,
#' \eqn{\Omega_b} and \eqn{\Omega_q}. Each combines mortality \eqn{g}, the
#' corresponding emigration rate (\eqn{\sigma_b} or \eqn{\sigma_q}), the
#' emigration-loss fraction \eqn{\mu}, and its dispersal matrix (\eqn{K_b} or
#' \eqn{K_q}).
#'
#' @section Dynamics:
#' Initialization and output parsing expose six state variables: \eqn{B_u},
#' \eqn{Q_u}, \eqn{B_y}, \eqn{Q_y}, \eqn{B_z}, and \eqn{Q_z}. The current
#' [dMYdt.BQm] method instead refers to aggregate \eqn{B} and \eqn{Q}, which
#' are not returned by [get_MY_vars.BQm]. Therefore its equations do not
#' currently correspond to the initialized six-state vector.
#'
#' @section Infective Biting:
#' The patch-level human biting density, `fqM`, is based on
#' \eqn{B_u+B_y+B_z}; the infective biting density, `fqZ`, is based on
#' \eqn{B_z}. [F_fqM.BQm] and [F_fqZ.BQm] compute these quantities.
#' [ramp.xds::compute_EIR] combines `fqZ` with the biting distribution matrix
#' and resident-bite fraction to compute stratum-level `EIR`;
#' [ramp.xds::F_hbr] uses `fqM` to compute stratum-level `HBR`.
#'
#' @section Egg Laying:
#' Gravid mosquitoes lay eggs at rate \eqn{\nu}; the oviposition rate and
#' number of eggs per batch determine patch-level egg production:
#' \deqn{G=\nu(Q_u+Q_y+Q_z)\Gamma.}
#' [F_eggs.BQm] computes \eqn{G}, which is passed to [ramp.xds::EggLaying]
#' for distribution among aquatic habitats.
#'
#' @seealso
#' + Basic setup: [make_MY_obj_BQm], [setup_MY_obj.BQm]
#' + Dynamics: [dMYdt.BQm]
#' + Biting and eggs: [F_fqZ.BQm], [F_fqM.BQm], [F_eggs.BQm]
#' + Parameters: [get_MY_pars.BQm], [change_MY_pars.BQm]
#' + Variables: [make_MY_inits_BQm], [setup_MY_inits.BQm], [parse_MY_orbits.BQm]
#'
#' @name BQm
#' @rdname BQm
NULL

#' @title Compute the derivatives for `BQm` (**MY**)
#' @description Compute the derivatives for the `BQm` module:
#' \deqn{\begin{array}{rl}
#' dB/dt &= \Lambda+\nu Q-fB-\Omega_b\cdot B \\
#' dQ/dt &= fB-\nu Q-\Omega_q\cdot Q
#' \end{array}}
#' @inheritSection BQm State Variables
#' @inheritSection BQm Dynamical Terms
#' @inheritSection BQm Parameters
#' @inheritSection BQm Dynamics
#' @inheritParams ramp.xds::dMYdt
#' @seealso [BQm], [F_fqZ.BQm], [F_eggs.BQm]
#' @note Derivatives are returned in blocks for `B`, then `Q`, each in patch order.
#' @return a [numeric] vector of derivatives in the module's state-variable order
#' @keywords internal
#' @export
dMYdt.BQm <- function(t, y, xds_obj, s){

  Lambda = xds_obj$terms$Lambda[[s]]

  with(get_MY_vars(y, xds_obj, s),{
    with(xds_obj$MY_obj[[s]],{

      dB <- Lambda + nu*Q - f*B - Omega_b %*% B
      dQ <- f*B - nu*Q - Omega_q %*% Q

      return(c(dB, dQ))
    })
  })
}

#' @title Get the state variables for `BQm` (**MY**)
#' @description Extract the state variables for the `BQm` module
#' from `y` and return them by name in a list for use during model solving.
#' @inheritSection BQm State Variables
#' @inheritParams ramp.xds::get_MY_vars
#' @return a named [list] of state variables and module-specific derived quantities
#' @keywords internal
#' @export
#' @note Variable indices are stored at `xds_obj$MY_obj[[s]]$ix`.
get_MY_vars.BQm <- function(y, xds_obj, s){
  with(xds_obj$MY_obj[[s]]$ix,
       return(list(
         Bu = y[Bu_ix],
         Qu = y[Qu_ix],
         By = y[By_ix],
         Qy = y[Qy_ix],
         Bz = y[Bz_ix],
         Qz = y[Qz_ix]
)))}

#' @title Mosquito bionomics for `BQm` (**MY**)
#' @description Implements [MBionomics] for the BQm model
#' @inheritSection BQm Parameters
#' @inheritSection BQm Survival and Dispersal
#' @inheritParams ramp.xds::MBionomics
#' @seealso [BQm], [MEffectSizes.BQm]
#' @return an **`xds`** model object
#' @keywords internal
#' @export
MBionomics.BQm <- function(t, y, xds_obj, s) {
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


#' @title Apply effect sizes for `BQm` (**MY**)
#' @description Implements [MEffectSizes] for the BQm model
#' @inheritSection BQm Parameters
#' @inheritSection BQm Survival and Dispersal
#' @inheritParams ramp.xds::MEffectSizes
#' @seealso [BQm], [MBionomics.BQm]
#' @return an **`xds`** model object
#' @keywords internal
#' @export
MEffectSizes.BQm <- function(t, y, xds_obj, s) {
  with(xds_obj$MY_obj[[s]], {
    xds_obj$MY_obj[[s]]$f <- f_t * es_f
    xds_obj$MY_obj[[s]]$q <- q_t * es_q
    xds_obj$MY_obj[[s]]$g <- g_t * es_g
    xds_obj$MY_obj[[s]]$sigma_b <- sigma_b_t * es_sigma_b
    xds_obj$MY_obj[[s]]$sigma_q <- sigma_q_t * es_sigma_q
    xds_obj$MY_obj[[s]]$mu <- mu_t * es_mu
    xds_obj$MY_obj[[s]]$nu <- nu_t * es_nu
    xds_obj$MY_obj[[s]]$eip <- eip_t * es_eip
    xds_obj <- change_Omega_b(xds_obj, s)
    xds_obj <- change_Omega_q(xds_obj, s)
    return(xds_obj)
  })
}

#' @title Compute infective mosquito biting density for `BQm` (**MY**)
#' @description Return the daily infective mosquito biting density on humans,
#' `fqZ`, for each patch in the `BQm` module.
#' Use only the blood-feeding behavioral states: `fqZ = f*q*Bz`.
#' @inheritSection BQm Infective Biting
#' @inheritParams ramp.xds::F_fqZ
#' @seealso [BQm], [F_fqM.BQm]
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
#' @note The current implementation obtains `q` using `get_f`
#' rather than `get_q`; this discrepancy requires repair.
F_fqZ.BQm <- function(t, y, xds_obj, s) {
  f = get_f(xds_obj, s)
  q = get_f(xds_obj, s)
  fqZ = f*q*y[xds_obj$MY_obj[[s]]$ix$Bz_ix]
  return(fqZ)
}

#' @title Compute mosquito biting density for `BQm` (**MY**)
#' @description Return the daily mosquito biting density on humans,
#' `fqM`, for each patch in the `BQm` module.
#' Use only the blood-feeding behavioral states: `fqM = f*q*(Bu+By+Bz)`.
#' @inheritSection BQm Infective Biting
#' @inheritParams ramp.xds::F_fqM
#' @seealso [BQm], [ramp.xds::F_hbr], [F_fqZ.BQm]
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
#' @note The current implementation obtains `q` using `get_f`
#' rather than `get_q`; this discrepancy requires repair.
F_fqM.BQm <- function(t, y, xds_obj, s) {
  f = get_f(xds_obj, s)
  q = get_f(xds_obj, s)
  M = with(xds_obj$MY_obj[[s]]$ix, y[Bu_ix] + y[By_ix] + y[Bz_ix])
  fqM = f*q*M
  return(fqM)
}

#' @title Compute egg production for `BQm` (**MY**)
#' @description Return the number of eggs laid by the adult mosquito population, per patch
#' @inheritSection BQm Egg Laying
#' @inheritParams ramp.xds::F_eggs
#' @seealso [BQm], [ramp.xds::EggLaying]
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
F_eggs.BQm <- function(t, y, xds_obj, s) {

  G <- with(get_MY_vars(y, xds_obj, s), Qu + Qy + Qz)
  eggs = with(xds_obj$MY_obj[[s]], {
    return(G*nu*eggsPerBatch)
  })
  return(eggs)
}

#' @title Get the feeding rate for `BQm` (**MY**)
#' @description Return the stored realized feeding rate per patch.
#' @inheritSection BQm Parameters
#' @inheritParams ramp.xds::get_f
#' @seealso [BQm]
#' @details Read `f` without evaluating ports or applying effect sizes.
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
get_f.BQm = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], f)
}

#' @title Get the human blood fraction for `BQm` (**MY**)
#' @description Return the stored realized human blood fraction per patch.
#' @inheritSection BQm Parameters
#' @inheritParams ramp.xds::get_q
#' @seealso [BQm]
#' @details Read `q` without evaluating ports or applying effect sizes.
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
get_q.BQm = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], q)
}

#' @title Get the adult mosquito mortality rate for `BQm` (**MY**)
#' @description Return the stored realized adult mosquito mortality rate per patch.
#' @inheritSection BQm Parameters
#' @inheritParams ramp.xds::get_g
#' @seealso [BQm]
#' @details Read `g` without evaluating ports or applying effect sizes.
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
get_g.BQm = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], g)
}

#' @title Get the patch emigration rate for `BQm` (**MY**)
#' @description Return the stored realized patch emigration rate per patch.
#' @inheritSection BQm Parameters
#' @inheritParams ramp.xds::get_sigma
#' @seealso [BQm]
#' @details Read `sigma_b` without evaluating ports or applying effect sizes.
#' This is the blood-feeding emigration rate, not the egg-laying emigration rate.
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
get_sigma.BQm = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], sigma_b)
}


#' @title Make a module object for `BQm` (**MY**)
#' @description Create a `BQm` module object and initialize its parameters.
#' Inspect parameters with [get_MY_pars.BQm] and change supported values
#' with [change_MY_pars.BQm]. Values in `options` override parameter arguments.
#' @seealso [get_MY_pars.BQm], [change_MY_pars.BQm]
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
#' @seealso [BQm], [setup_MY_obj.BQm]
#' @keywords internal
#' @importFrom ramp.xds checkIt
#' @export
make_MY_obj_BQm = function(nPatches, options=list(), eip=12,
                          g=1/12, sigma_b=1/8, sigma_q=1/8, mu=0, f=0.5, q=0.95,
                          nu=1, eggsPerBatch=60){

  with(options,{
    MY_obj <- list()
    class(MY_obj) <- "BQm"

    MY_obj$nPatches <- nPatches

    MY_obj <- setup_eip_obj(checkIt(eip, nPatches), MY_obj)
    MY_obj <- setup_f_obj(checkIt(f, nPatches), MY_obj)
    MY_obj <- setup_q_obj(checkIt(q, nPatches), MY_obj)
    MY_obj <- setup_g_obj(checkIt(g, nPatches), MY_obj)
    MY_obj <- setup_mu_obj(checkIt(mu, nPatches), MY_obj)
    MY_obj <- setup_nu_obj(checkIt(nu, nPatches), MY_obj)
    MY_obj <- setup_sigma_q_obj(checkIt(sigma_q, nPatches), MY_obj)
    MY_obj <- setup_sigma_b_obj(checkIt(sigma_b, nPatches), MY_obj)

    MY_obj$eggsPerBatch <- eggsPerBatch


    return(MY_obj)
})}


#' @title Set up the module object for `BQm` (**MY**)
#' @description
#'
#' Set up a `BQm` module object for the **MY** component and the \eqn{s^{th}} vector species.
#'
#' + Set the module name to `BQm`.
#' + Create the module object with [make_MY_obj_BQm], passing `nPatches` and `options`.
#' + Initialize the circadian function using [ramp.xds::setup_F_circadian.setup].
#' + Set the skill set using [ramp.xds::get_skillset_MY].
#' + Initialize the mosquito dispersal matrices using [ramp.xds::setup_K_matrix.zero],
#'   with `which_K = "Kb"` and `which_K = "Kq"`.
#' + Set initial values using [ramp.xds::setup_MY_inits], passing `options` to override defaults.
#' + Return the updated **`xds`** model object.
#'
#' @inheritSection BQm State Variables
#' @inheritSection BQm Parameters
#' @seealso [BQm], [make_MY_obj_BQm], [setup_MY_inits.BQm]
#'
#' @note The module object is stored at `xds_obj$MY_obj[[s]]`
#'
#' @inheritParams ramp.xds::setup_MY_obj
#' @return an **`xds`** model object
#' @keywords internal
#' @export
setup_MY_obj.BQm = function(MYname, xds_obj, options=list(), s=1){
  xds_obj$MYname <- "BQm"
  xds_obj$MY_obj[[s]] = make_MY_obj_BQm(xds_obj$nPatches, options)
  xds_obj <- setup_F_circadian("setup", xds_obj, s=s)
  xds_obj <- setup_K_matrix("zero", xds_obj, options=list(which_K = "Kb"), s=s)
  xds_obj <- setup_K_matrix("zero", xds_obj, options=list(which_K = "Kq"), s=s)
  xds_obj <- setup_MY_inits(xds_obj, options=options, s=s)
  return(xds_obj)
}


#' @title Get parameters for `BQm` (**MY**)
#'
#' @description Return the parameters for the `BQm` module as a named list.
#' @inheritSection BQm Parameters
#'
#' @inheritParams ramp.xds::get_MY_pars
#'
#' @return a named [list] of parameter values
#' @seealso [BQm], [change_MY_pars.BQm]
#' @seealso [make_MY_obj_BQm]
#' @keywords internal
#' @export
get_MY_pars.BQm <- function(xds_obj, s=1) {
  with(xds_obj$MY_obj[[s]], list(
    f=f_t, q=q_t, g=g_t, sigma=sigma_t, eip=eip, mu=mu_t,
    nu=nu_t, eggsPerBatch=eggsPerBatch,
    Kb_matrix=Kb_matrix, Kq_matrix=Kq_matrix
  ))
}


#' @title Change parameters for `BQm` (**MY**)
#' @seealso [make_MY_obj_BQm], [get_MY_pars.BQm]
#' @description Change parameter values named in `options` for the `BQm`
#' module, leaving other values unchanged.
#' Refresh baseline bionomics at time zero using [ramp.xds::MBionomics],
#' without calling [ramp.xds::MEffectSizes].
#' @inheritSection BQm Parameters
#' @inheritParams ramp.xds::change_MY_pars
#' @return an **`xds`** model object
#' @keywords internal
#' @export
change_MY_pars.BQm <- function(xds_obj, options=list(), s=1) {
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

#' @title Make initial values for `BQm` (**MY**)
#' @description Create a named list of initial values for the `BQm` module.
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
#' @inheritSection BQm State Variables
#' @seealso [BQm], [setup_MY_inits.BQm]
#' @return a named [list] of initial values
#' @keywords internal
#' @export
#' @seealso [change_MY_inits.BQm]
make_MY_inits_BQm = function(nPatches, options = list(),
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

#' @title Set up initial values for `BQm` (**MY**)
#' @description Set up initial values for the `BQm` module using
#' [make_MY_inits_BQm]. Values in `options` override the defaults.
#' @inheritSection BQm State Variables
#' @inheritParams ramp.xds::setup_MY_inits
#' @return an **`xds`** model object
#' @keywords internal
#' @export
#' @seealso [change_MY_inits.BQm]
#' @note Initial values are stored at `xds_obj$MY_obj[[s]]$inits`.
#' @seealso [ramp.xds::get_MY_inits]
setup_MY_inits.BQm = function(xds_obj, options=list(), s=1){
  xds_obj$MY_obj[[s]]$inits = make_MY_inits_BQm(xds_obj$nPatches, options)
  return(xds_obj)
}


#' @title Change initial values for `BQm` (**MY**)
#' @description Change initial values named in `options` for the `BQm`
#' module, leaving other values unchanged.
#' @inheritParams ramp.xds::change_MY_inits
#' @return an **`xds`** model object
#' @keywords internal
#' @export
#' @seealso [make_MY_inits_BQm], [setup_MY_inits.BQm]
#' @note Stored initial values for this module are at
#' `xds_obj$MY_obj[[s]]$inits`.
change_MY_inits.BQm <- function(xds_obj, options=list(), s=1) {
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

#' @title Set up variable indices for `BQm` (**MY**)
#' @description Assign consecutive variable indices for the `BQm` module,
#' starting after `xds_obj$max_ix`, and update the maximum index.
#' Allocate blocks in the order `Bu_ix`, `Qu_ix`, `By_ix`, `Qy_ix`, `Bz_ix`, `Qz_ix`.
#' Each block has length `nPatches`.
#' @inheritSection BQm State Variables
#' @inheritParams ramp.xds::setup_MY_ix
#' @return an **`xds`** model object
#' @importFrom utils tail
#' @keywords internal
#' @export
#' @note Indices are stored at `xds_obj$MY_obj[[s]]$ix`.
#' The maximum index is stored at `xds_obj$max_ix`.
setup_MY_ix.BQm <- function(xds_obj, s) {with(xds_obj,{

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

#' @title Parse outputs for `BQm` (**MY**)
#' @description Extract the model state outputs for the `BQm` module
#' from `outputs` and return them by name in a list.
#' Sum behavioral states to compute total density `M`, infected density `Y`
#' (exposed and infectious), and infectious density `Z`.
#' Also compute `y = Y/M` and `z = Z/M`.
#' @note Variable indices are stored at `xds_obj$MY_obj[[s]]$ix`.
#' @inheritSection BQm State Variables
#' @inheritSection BQm Infective Biting
#' @inheritParams ramp.xds::parse_MY_orbits
#' @seealso [BQm], [get_MY_vars.BQm]
#' @return a [list]
#' @keywords internal
#' @export
parse_MY_orbits.BQm <- function(outputs, xds_obj, s) {
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

#' @title Check `BQm` (**MY**)
#'
#' @inheritParams ramp.xds::check_MY
#'
#' @return an **`xds`** model object
#' @keywords internal
#' @export
check_MY.BQm = function(xds_obj, s){
  return(xds_obj)
}

#' @title Get the skill set for `BQm` (**MY**)
#'
#' @description Return the capabilities and compatibilities of the `BQm` module as a list.
#'
#' @inheritParams ramp.xds::get_skillset_MY
#'
#' @return the **MY** module skill set, as a list
#'
#' @keywords internal
#' @export
get_skillset_MY.BQm = function(xds_obj, s){
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
