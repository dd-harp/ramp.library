# specialized methods for the adult mosquito SEI model
# TODO: Repair steady_state_MY.SEI: M_eq is referenced but never computed.
# Review returned state coverage and verify the equilibrium before using it.
# TODO: Reconcile Update_MYt.SEI's list return with the numeric-vector interface,
# including state coverage/order.
# TODO: Review Update_MYt.SEI's top-level Lambda/kappa reads, MYday and ccc
# dependencies, and use of Omega as a one-step matrix rather than a rate matrix.

#' @title The `SEI` module for the **MY** component
#' @description
#' Implements the **MY** component using a Susceptible-Exposed-Infectious
#' (SEI) compartmental model of adult mosquito infection dynamics. Mosquito
#' abundance, exposed mosquitoes, and infectious mosquitoes are tracked
#' separately across patches.
#'
#' @section State Variables:
#' \describe{
#'   \item{`M`}{density of adult mosquitoes}
#'   \item{`Y`}{density of infected (exposed) adult mosquitoes}
#'   \item{`Z`}{density of infectious adult mosquitoes}
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
#'   \item{`sigma` or \eqn{\sigma}}{emigration rate --- use [ramp.xds::setup_sigma] or [ramp.xds::setup_mosy_emigration]}
#'   \item{`mu` or \eqn{\mu}}{fraction of emigrating mosquitoes lost --- use [ramp.xds::setup_mu] or [ramp.xds::setup_emigration_loss]}
#'   \item{`K` or `K_matrix`}{mosquito dispersal matrix --- use [ramp.xds::setup_K_matrix]}
#'   \item{`eip` or \eqn{\tau}}{extrinsic incubation period --- use [ramp.xds::setup_eip]}
#'   \item{`nu` or \eqn{\nu}}{oviposition rate: the number of egg batches laid in each patch, per mosquito, per day --- use [ramp.xds::setup_nu] or [ramp.xds::setup_egg_batch_rates]}
#'   \item{`eggsPerBatch` or \eqn{\Gamma}}{number of eggs per batch}
#' }
#' After basic setup, [change_MY_pars.SEI] can be used to change parameter
#' values. Bionomic parameters are configured as static baselines by default.
#'
#' @section Survival and Dispersal:
#' Mosquito survival and dispersal are described by the demographic matrix
#' \eqn{\Omega}. The parameter \eqn{g} is the per-capita mortality rate while
#' mosquitoes are in a patch. Mosquitoes emigrate from patches at per-capita
#' rates \eqn{\sigma}, and a fraction \eqn{\mu} of emigrating mosquitoes is
#' lost. The total rate of loss from a patch is \deqn{g+\mu\sigma.}
#' The columns of the dispersal matrix \eqn{K} sum to zero. For non-trivial
#' dispersal matrices, \eqn{\mbox{diag}(K)=-1}; its off-diagonal elements
#' describe the fractions of surviving emigrants that arrive in other patches.
#' The demographic matrix is
#' \deqn{\Omega=\mbox{diag}(g+\sigma\mu)-K\cdot\mbox{diag}(\sigma(1-\mu)).}
#'
#' @section Blood Feeding:
#' The blood feeding rate \eqn{f} and human fraction \eqn{q} determine the
#' rate of human blood feeding, \eqn{fq}. Natural constraints are
#' \eqn{f=0} when no blood hosts are available, \eqn{q=0} when no humans are
#' available, and \eqn{q=1} when humans are the only available hosts.
#'
#' @section Dynamics:
#' Adult mosquitoes emerge from aquatic habitats at rate \eqn{\Lambda}.
#' Mosquitoes become exposed when they take an infectious human blood meal;
#' exposed mosquitoes progress to the infectious class at rate
#' \eqn{1/\tau}, where \eqn{\tau} is the EIP. The model equations are
#' \deqn{
#' \begin{array}{rl}
#' dM/dt &= \Lambda-\Omega\cdot M \\\\
#' dY/dt &= fq\kappa(M-Y)-\Omega\cdot Y \\\\
#' dZ/dt &= (Y-Z)/\tau-\Omega\cdot Z
#' \end{array}
#' }
#'
#' @section DTS:
#' [Update_MYt.SEI] updates `M` and `Z` using the stored `Omega` directly.
#' The updating equations use `D = xds_obj$MYday` and the module field `ccc`.
#' The method returns these two updated states as a list; it does not update `Y`.
#'
#' @section Infective Biting:
#' The patch-level human biting density, `fqM`, is \eqn{f q M}; the
#' infective biting density, `fqZ`, is \eqn{f q Z}. [F_fqM.SEI] and
#' [F_fqZ.SEI] compute these quantities from the total and infectious mosquito
#' densities, respectively. [ramp.xds::compute_EIR] combines `fqZ` with the
#' biting distribution matrix and resident-bite fraction to compute
#' stratum-level `EIR`; [ramp.xds::F_hbr] uses `fqM` to compute stratum-level
#' `HBR`.
#'
#' @section Egg Laying:
#' The oviposition rate \eqn{\nu} and the number of eggs per batch
#' \eqn{\Gamma} determine the egg production rate per patch:
#' \deqn{G=\nu\Gamma M.}
#' [F_eggs.SEI] computes \eqn{G}, which is passed to [ramp.xds::EggLaying]
#' for distribution among aquatic habitats.
#'
#' @seealso
#' + Basic setup: [make_MY_obj_SEI], [setup_MY_obj.SEI]
#' + Dynamics: [dMYdt.SEI], [Update_MYt.SEI]
#' + Biting and eggs: [F_fqZ.SEI], [F_fqM.SEI], [F_eggs.SEI]
#' + Parameters: [get_MY_pars.SEI], [change_MY_pars.SEI]
#' + Variables: [make_MY_inits_SEI], [setup_MY_inits.SEI], [parse_MY_orbits.SEI]
#'
#' @name SEI
#' @rdname SEI
NULL

#' @title Compute the derivatives for `SEI` (**MY**)
#' @description Compute the derivatives for the `SEI` module:
#' \deqn{\begin{array}{rl}
#' dM/dt &= \Lambda-\Omega\cdot M \\
#' dY/dt &= fq\kappa(M-Y)-\Omega\cdot Y \\
#' dZ/dt &= (Y-Z)/\tau-\Omega\cdot Z
#' \end{array}}
#' Here \eqn{\tau} is the extrinsic incubation period (`eip`).
#' @inheritSection SEI State Variables
#' @inheritSection SEI Dynamical Terms
#' @inheritSection SEI Parameters
#' @inheritSection SEI Dynamics
#'
#' @inheritParams ramp.xds::dMYdt
#' @seealso [SEI], [F_fqZ.SEI], [F_eggs.SEI]
#' @note Derivatives are returned in blocks for `M`, `Y`, and `Z`, each in patch order.
#' @return a [numeric] vector of derivatives in the module's state-variable order
#' @keywords internal
#' @export
dMYdt.SEI <- function(t, y, xds_obj, s) {
  Lambda = xds_obj$terms$Lambda[[s]]
  kappa = xds_obj$terms$kappa[[s]]

  with(get_MY_vars(y, xds_obj, s),{
    with(xds_obj$MY_obj[[s]],{

      dM <- Lambda - (Omega %*% M)
      dY <- f*q*kappa*(M-Y) - Omega %*% Y
      dZ <- (Y-Z)/eip - (Omega %*% Z)

      return(c(dM, dY, dZ))
    })
  })
}


#' @title Update states for `SEI` (**MY**)
#' @description Update the state variables for the `SEI` module:
#' \deqn{\begin{array}{rl}
#' M_t &= c\Lambda D+\Omega\cdot M_{t-1} \\\\
#' Z_t &= \Omega\cdot\left[(1-e^{-fq\kappa D})(M_{t-1}-Z_{t-1})\right]
#'       +\Omega\cdot Z_{t-1}
#'       +(1-e^{-c})(1-e^{-fq\kappa D})\Lambda cD
#' \end{array}}
#' Here \eqn{c} is `ccc` and \eqn{D} is `xds_obj$MYday`.
#' @inheritSection SEI State Variables
#' @inheritSection SEI Dynamical Terms
#' @inheritSection SEI Parameters
#' @inheritSection SEI DTS
#' @inheritParams ramp.xds::Update_MYt
#' @seealso [SEI], [dMYdt.SEI]
#' @return a [list]
#' @keywords internal
#' @export
Update_MYt.SEI <- function(t, y, xds_obj, s) {
  Lambda = xds_obj$Lambda[[s]]
  kappa = xds_obj$kappa[[s]]
  D = xds_obj$MYday

  with(get_MY_vars(y, xds_obj, s),{
    with(xds_obj$MY_obj[[s]],{
      Mt <- ccc*Lambda*D + Omega %*% M
      Zt <- Omega%*%((1-exp(-f*q*kappa*xds_obj$MYday))*(M-Z)) + (Omega %*% Z)
      Zt <- Zt + (1-exp(-ccc))*(1-exp(-f*q*kappa*xds_obj$MYday))*Lambda*ccc*D
      return(list(M=unname(Mt), Z=unname(Zt)))
    })
  })
}

#' @title Get the state variables for `SEI` (**MY**)
#' @description Extract the state variables for the `SEI` module
#' from `y` and return them by name in a list for use during model solving.
#' @inheritSection SEI State Variables
#' @inheritParams ramp.xds::get_MY_vars
#' @return a named [list] of state variables and module-specific derived quantities
#' @keywords internal
#' @export
#' @note Variable indices are stored at `xds_obj$MY_obj[[s]]$ix`.
get_MY_vars.SEI <- function(y, xds_obj, s){
  with(xds_obj$MY_obj[[s]]$ix,
       return(list(
         M = y[M_ix],
         Y = y[Y_ix],
         Z = y[Z_ix]
       )))
}


#' @title Mosquito bionomics for `SEI` (**MY**)
#' @description Implements [MBionomics] for the SEI model
#' @inheritSection SEI Parameters
#' @inheritSection SEI Survival and Dispersal
#' @inheritSection SEI Blood Feeding
#' @inheritParams ramp.xds::MBionomics
#' @seealso [SEI], [MEffectSizes.SEI]
#' @return an **`xds`** model object
#' @keywords internal
#' @export
MBionomics.SEI <- function(t, y, xds_obj, s) {
  xds_obj$MY_obj[[s]]$f_t <- F_f(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$q_t <- F_q(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$g_t <- F_g(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$sigma_t <- F_sigma(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$mu_t <- F_mu(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$nu_t <- F_nu(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$eip_t <- F_eip(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$K_matrix_t <- F_K_matrix(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$es_f <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_q <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_g <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_sigma <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_mu <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_nu <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_eip <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_K_matrix <- matrix(1, xds_obj$nPatches, xds_obj$nPatches)
  return(xds_obj)
}


#' @title Apply effect sizes for `SEI` (**MY**)
#' @description Implements [MEffectSizes] for the SEI model
#' @inheritSection SEI Parameters
#' @inheritSection SEI Survival and Dispersal
#' @inheritParams ramp.xds::MEffectSizes
#' @seealso [SEI], [MBionomics.SEI]
#' @return an **`xds`** model object
#' @keywords internal
#' @export
MEffectSizes.SEI <- function(t, y, xds_obj, s) {
  with(xds_obj$MY_obj[[s]], {
    xds_obj$MY_obj[[s]]$f <- f_t * es_f
    xds_obj$MY_obj[[s]]$q <- q_t * es_q
    xds_obj$MY_obj[[s]]$g <- g_t * es_g
    xds_obj$MY_obj[[s]]$sigma <- sigma_t * es_sigma
    xds_obj$MY_obj[[s]]$mu <- mu_t * es_mu
    xds_obj$MY_obj[[s]]$nu <- nu_t * es_nu
    xds_obj$MY_obj[[s]]$eip <- eip_t * es_eip
    xds_obj$MY_obj[[s]]$K_matrix <- K_matrix_t * es_K_matrix
    xds_obj <- change_Omega(xds_obj, s)
    return(xds_obj)
  })
}

#' @title Compute infective mosquito biting density for `SEI` (**MY**)
#' @description Return the daily infective mosquito biting density on humans,
#' `fqZ`, for each patch in the `SEI` module.
#' Compute `fqZ = f*q*Z`.
#' @inheritSection SEI Infective Biting
#' @inheritParams ramp.xds::F_fqZ
#' @seealso [SEI], [F_fqM.SEI]
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
F_fqZ.SEI <- function(t, y, xds_obj, s) {
  f = get_f(xds_obj, s)
  q = get_q(xds_obj, s)
  Z = y[xds_obj$MY_obj[[s]]$ix$Z_ix]
  return(f*q*Z)
}

#' @title Compute mosquito biting density for `SEI` (**MY**)
#' @description Return the daily mosquito biting density on humans,
#' `fqM`, for each patch in the `SEI` module.
#' Compute `fqM = f*q*M`.
#' @inheritSection SEI Infective Biting
#' @inheritParams ramp.xds::F_fqM
#' @seealso [SEI], [ramp.xds::F_hbr], [F_fqZ.SEI]
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
F_fqM.SEI <- function(t, y, xds_obj, s) {
  f = get_f(xds_obj, s)
  q = get_q(xds_obj, s)
  M = y[xds_obj$MY_obj[[s]]$ix$M_ix]
  return(f*q*M)
}

#' @title Compute egg production for `SEI` (**MY**)
#' @description Return the number of eggs laid by the adult mosquito population, per patch
#' @inheritSection SEI Egg Laying
#' @inheritParams ramp.xds::F_eggs
#' @seealso [SEI], [ramp.xds::EggLaying]
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
F_eggs.SEI <- function(t, y, xds_obj, s) {
  M <- y[xds_obj$MY_obj[[s]]$ix$M_ix]
  with(xds_obj$MY_obj[[s]],{
    return(M*nu*eggsPerBatch)
  })
}


#' @title Get the feeding rate for `SEI` (**MY**)
#' @description Return the stored realized feeding rate per patch.
#' @inheritSection SEI Parameters
#' @inheritParams ramp.xds::get_f
#' @seealso [SEI]
#' @details Read `f` without evaluating ports or applying effect sizes.
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
get_f.SEI = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], f)
}

#' @title Get the human blood fraction for `SEI` (**MY**)
#' @description Return the stored realized human blood fraction per patch.
#' @inheritSection SEI Parameters
#' @inheritParams ramp.xds::get_q
#' @seealso [SEI]
#' @details Read `q` without evaluating ports or applying effect sizes.
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
get_q.SEI = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], q)
}

#' @title Get the adult mosquito mortality rate for `SEI` (**MY**)
#' @description Return the stored realized adult mosquito mortality rate per patch.
#' @inheritSection SEI Parameters
#' @inheritParams ramp.xds::get_g
#' @seealso [SEI]
#' @details Read `g` without evaluating ports or applying effect sizes.
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
get_g.SEI = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], g)
}

#' @title Get the patch emigration rate for `SEI` (**MY**)
#' @description Return the stored realized patch emigration rate per patch.
#' @inheritSection SEI Parameters
#' @inheritParams ramp.xds::get_sigma
#' @seealso [SEI]
#' @details Read `sigma` without evaluating ports or applying effect sizes.
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
get_sigma.SEI = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], sigma)
}


#' @title Make a module object for `SEI` (**MY**)
#' @description Create an `SEI` module object and initialize its parameters.
#' Inspect parameters with [get_MY_pars.SEI] and change supported values
#' with [change_MY_pars.SEI]. Values in `options` override parameter arguments.
#' @seealso [get_MY_pars.SEI], [change_MY_pars.SEI]
#' @param nPatches is the number of patches, an integer
#' @param options a [list] of values that overrides the defaults
#' @param eip the extrinsic incubation period
#' @param g mosquito mortality rate
#' @param sigma emigration rate
#' @param mu emigration loss
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu the oviposition rate: the number of egg batches laid in each patch,
#' per mosquito, per day
#' @param eggsPerBatch eggs laid per oviposition
#' @return a [list]
#' @seealso [SEI], [setup_MY_obj.SEI]
#' @keywords internal
#'
#' @export
make_MY_obj_SEI = function(nPatches, options=list(), eip =12,
                             g=1/12,  sigma=1/8,  mu=0,
                             f=0.3,  q=0.95,
                             nu=1,  eggsPerBatch=60){
    with(options,{
      MY_obj <- list()
      MY_obj$nPatches <- nPatches

      eip_par <- 'static'
      class(eip_par) <- 'static'

      MY_obj <- setup_eip_obj(checkIt(eip, nPatches), MY_obj)
      MY_obj <- setup_f_obj(checkIt(f, nPatches), MY_obj)
      MY_obj <- setup_q_obj(checkIt(q, nPatches), MY_obj)
      MY_obj <- setup_g_obj(checkIt(g, nPatches), MY_obj)
      MY_obj <- setup_mu_obj(checkIt(mu, nPatches), MY_obj)
      MY_obj <- setup_nu_obj(checkIt(nu, nPatches), MY_obj)
      MY_obj <- setup_sigma_obj(checkIt(sigma, nPatches), MY_obj)
      MY_obj$eggsPerBatch <- eggsPerBatch

      class(MY_obj) = "SEI"

      return(MY_obj)
})}

#' @title Set up the module object for `SEI` (**MY**)
#' @description
#'
#' Set up an `SEI` module object for the **MY** component and the \eqn{s^{th}} vector species.
#'
#' + Set the module name to `SEI`.
#' + Create the module object with [make_MY_obj_SEI], passing `nPatches` and `options`.
#' + Set the skill set using [ramp.xds::get_skillset_MY].
#' + Set initial values using [ramp.xds::setup_MY_inits], passing `options` to override defaults.
#' + Initialize the circadian function using [ramp.xds::setup_F_circadian.setup].
#' + Initialize the mosquito demographic matrix port.
#' + Initialize the mosquito dispersal matrix using [ramp.xds::setup_K_matrix.zero].
#' + Return the updated **`xds`** model object.
#'
#' @inheritSection SEI State Variables
#' @inheritSection SEI Parameters
#' @seealso [SEI], [make_MY_obj_SEI], [setup_MY_inits.SEI]
#'
#' @note The module object is stored at `xds_obj$MY_obj[[s]]`
#' Its classes are `SEI` and `SEI_` followed by `xds_obj$xds`.
#'
#' @inheritParams ramp.xds::setup_MY_obj
#' @return an **`xds`** model object
#' @keywords internal
#' @export
setup_MY_obj.SEI = function(MYname, xds_obj, options=list(), s=1){
  xds_obj$MYname <- "SEI"
  MY_obj <- make_MY_obj_SEI(xds_obj$nPatches, options)
  class(MY_obj) <- c("SEI", paste("SEI_", xds_obj$xds, sep=""))
  xds_obj$MY_obj[[s]] <- MY_obj
  xds_obj <- setup_MY_inits(xds_obj, options=options, s=s)
  xds_obj <- setup_F_circadian("setup", xds_obj, s=s)
  xds_obj <- setup_Omega_obj("xde", "static", xds_obj, s=s)
  xds_obj <- setup_Upsilon_obj(xds_obj, s=s, FALSE)
  xds_obj <- setup_K_matrix("zero", xds_obj, s=s)
  return(xds_obj)
}

#' @title Get parameters for `SEI` (**MY**)
#' @seealso [make_MY_obj_SEI], [change_MY_pars.SEI]
#'
#' @description Return the parameters for the `SEI` module as a named list.
#' @inheritSection SEI Parameters
#'
#' @inheritParams ramp.xds::get_MY_pars
#'
#' @return a named [list] of parameter values
#' @keywords internal
#' @export
get_MY_pars.SEI <- function(xds_obj, s=1) {
  with(xds_obj$MY_obj[[s]], list(
    f=f_t, q=q_t, g=g_t, sigma=sigma_t, eip=eip, mu=mu_t,
    nu=nu_t, eggsPerBatch=eggsPerBatch, K_matrix=K_matrix
  ))
}


#' @title Change parameters for `SEI` (**MY**)
#' @seealso [make_MY_obj_SEI], [get_MY_pars.SEI]
#' @description Change parameter values named in `options` for the `SEI`
#' module, leaving other values unchanged.
#' Refresh baseline bionomics at time zero using [ramp.xds::MBionomics],
#' without calling [ramp.xds::MEffectSizes].
#' @inheritSection SEI Parameters
#' @inheritParams ramp.xds::change_MY_pars
#' @return an **`xds`** model object
#' @keywords internal
#' @export
change_MY_pars.SEI <- function(xds_obj, options=list(), s=1) {
  nHabitats <- xds_obj$nHabitats
  with(xds_obj$MY_obj[[s]], with(options,{
    xds_obj$MY_obj[[s]]$f_obj$f = f
    xds_obj$MY_obj[[s]]$q_obj$q = q
    xds_obj$MY_obj[[s]]$g_obj$g = g
    xds_obj$MY_obj[[s]]$sigma_obj$sigma = sigma
    xds_obj$MY_obj[[s]]$eip_obj$eip = eip
    xds_obj$MY_obj[[s]]$mu_obj$mu = mu
    xds_obj$MY_obj[[s]]$nu_obj$nu = nu
    xds_obj$MY_obj[[s]]$eggsPerBatch = eggsPerBatch

    y <- get_inits(xds_obj, flatten=TRUE)
    xds_obj <- MBionomics(0, y, xds_obj, s)

    return(xds_obj)
  }))}


#' @title Make initial values for `SEI` (**MY**)
#' @description Create a named list of initial values for the `SEI` module.
#' Values in `options` override the defaults.
#' @details Initial values must be numeric scalars or have length `nPatches`.
#' Scalars are expanded to that length.
#' @param nPatches the number of patches
#' @param options a [list] of values that overrides the defaults
#' @param M total mosquito density at each patch
#' @param Y total infected mosquito density at each patch
#' @param Z infectious mosquito density at each patch
#' @inheritSection SEI State Variables
#' @seealso [SEI], [setup_MY_inits.SEI]
#' @return a named [list] of initial values
#' @keywords internal
#' @export
#' @seealso [change_MY_inits.SEI]
make_MY_inits_SEI = function(nPatches, options = list(),
                               M=5, Y=1, Z=0){
  with(options,{
    M = checkIt(M, nPatches)
    Y = checkIt(Y, nPatches)
    Z = checkIt(Z, nPatches)
    return(list(M=M, Y=Y, Z=Z))
  })
}


#' @title Set up initial values for `SEI` (**MY**)
#' @description Set up initial values for the `SEI` module using
#' [make_MY_inits_SEI]. Values in `options` override the defaults.
#' @inheritSection SEI State Variables
#' @inheritParams ramp.xds::setup_MY_inits
#' @return an **`xds`** model object
#' @keywords internal
#' @export
#' @seealso [change_MY_inits.SEI]
#' @note Initial values are stored at `xds_obj$MY_obj[[s]]$inits`.
#' @seealso [ramp.xds::get_MY_inits]
setup_MY_inits.SEI = function(xds_obj, options=list(), s=1){
  xds_obj$MY_obj[[s]]$inits = make_MY_inits_SEI(xds_obj$nPatches, options)
  return(xds_obj)
}

#' @title Change initial values for `SEI` (**MY**)
#' @description Change initial values named in `options` for the `SEI`
#' module, leaving other values unchanged.
#' @inheritParams ramp.xds::change_MY_inits
#' @return an **`xds`** model object
#' @keywords internal
#' @export
#' @seealso [make_MY_inits_SEI], [setup_MY_inits.SEI]
#' @note Stored initial values for this module are at
#' `xds_obj$MY_obj[[s]]$inits`.
change_MY_inits.SEI <- function(xds_obj, options=list(), s=1) {
  with(xds_obj$MY_obj[[s]]$inits,
    with(options,{
      xds_obj$MY_obj[[s]]$inits$M = M
      xds_obj$MY_obj[[s]]$inits$Y = Y
      xds_obj$MY_obj[[s]]$inits$Z = Z
      return(xds_obj)
}))}



#' @title Set up variable indices for `SEI` (**MY**)
#' @description Assign consecutive variable indices for the `SEI` module,
#' starting after `xds_obj$max_ix`, and update the maximum index.
#' Allocate blocks in the order `M_ix`, `Y_ix`, `Z_ix`.
#' Each block has length `nPatches`.
#' @inheritSection SEI State Variables
#' @inheritParams ramp.xds::setup_MY_ix
#' @return an **`xds`** model object
#' @importFrom utils tail
#' @keywords internal
#' @export
#' @note Indices are stored at `xds_obj$MY_obj[[s]]$ix`.
#' The maximum index is stored at `xds_obj$max_ix`.
setup_MY_ix.SEI <- function(xds_obj, s) {with(xds_obj,{

  M_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(M_ix, 1)

  Y_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Y_ix, 1)

  Z_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Z_ix, 1)

  xds_obj$max_ix = max_ix
  xds_obj$MY_obj[[s]]$ix = list(M_ix=M_ix, Y_ix=Y_ix, Z_ix=Z_ix)
  return(xds_obj)
})}


#' @title Parse outputs for `SEI` (**MY**)
#' @description Extract the model state outputs for the `SEI` module
#' from `outputs` and return them by name in a list.
#' Also compute `y = Y/M`, `z = Z/M`, `fqZ = f*q*Z`, and `fqM = f*q*M`.
#' @note Variable indices are stored at `xds_obj$MY_obj[[s]]$ix`.
#' @inheritSection SEI State Variables
#' @inheritSection SEI Infective Biting
#' @inheritParams ramp.xds::parse_MY_orbits
#' @seealso [SEI], [get_MY_vars.SEI]
#' @return a [list]
#' @keywords internal
#' @export
parse_MY_orbits.SEI <- function(outputs, xds_obj, s) {with(xds_obj$MY_obj[[s]]$ix,{
  M = outputs[,M_ix]
  Y = outputs[,Y_ix]
  Z = outputs[,Z_ix]
  f = get_f(xds_obj,s)
  q = get_q(xds_obj,s)
  y = Y/M
  z = Z/M
  return(list(M=M, Z=Z, Y=Y, y=y, z=z, fqZ=f*q*Z, fqM=f*q*M))
})}

#' @title Check `SEI` (**MY**)
#'
#' @inheritParams ramp.xds::check_MY
#'
#' @return an **`xds`** model object
#' @keywords internal
#' @export
check_MY.SEI = function(xds_obj, s){
  return(xds_obj)
}

#' @title Get the skill set for `SEI` (**MY**)
#'
#' @description Return the capabilities and compatibilities of the `SEI` module as a list.
#'
#' @inheritParams ramp.xds::get_skillset_MY
#'
#' @return the **MY** module skill set, as a list
#'
#' @keywords internal
#' @export
get_skillset_MY.SEI = function(xds_obj, s){
  skills = list(
    f_is_port = TRUE,
    q_is_port = TRUE,
    g_is_port = TRUE,
    sigma_is_port = TRUE,
    sigma_b_is_port = FALSE,
    sigma_q_is_port = FALSE,
    mu_is_port = TRUE,
    nu_is_port = TRUE,
    eip_is_port = TRUE,
    K_matrix_is_port = TRUE,
    delay = FALSE
  )
  return(skills)
}

#' @title Compute steady states for `SEI` (**MY**)
#' @description Compute steady-state mosquito densities for the `SEI` module
#' given adult mosquito emergence `Lambda` and net infectiousness `kappa`
#' for the differential-equation system.
#' @note The current implementation references undefined `M_eq` and requires
#' repair before use. Its intended return includes `M` and `Z` but not `Y`.
#' @inheritSection SEI State Variables
#' @inheritSection SEI Parameters
#' @inheritSection SEI Dynamics
#' @inheritParams ramp.xds::steady_state_MY
#' @seealso [SEI]
#' @return a named [numeric] vector containing `M` followed by `Z`, once the
#' undefined `M_eq` reference is repaired
#' @keywords internal
#' @export
steady_state_MY.SEI = function(Lambda, kappa, xds_obj, s=1){
  with(xds_obj$MY_obj[[s]],{
    Omega_inv <- solve(Omega)
    Z_eq <- as.vector(solve(diag(f*q*kappa) + Omega) %*% diag(f*q*kappa) %*% M_eq)
    return(c(M=M_eq, Z=Z_eq))
})}
