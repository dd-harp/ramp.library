# specialized methods for the adult mosquito BQm model

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
#'   \item{`nu` or \eqn{\nu}}{egg batch rate --- use [ramp.xds::setup_nu] or [ramp.xds::setup_egg_batch_rates]}
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
#' Gravid mosquitoes lay eggs at rate \eqn{\nu}; the egg batch rate and
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

#' @title Compute derivatives for `BQm` (**MY**)
#' @description Implements [dMYdt] for the BQm ODE model.
#' @inheritSection BQm State Variables
#' @inheritSection BQm Dynamical Terms
#' @inheritSection BQm Parameters
#' @inheritSection BQm Dynamics
#' @inheritParams ramp.xds::dMYdt
#' @seealso [BQm], [F_fqZ.BQm], [F_eggs.BQm]
#' @return a [numeric] vector
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

#' @title The skill set
#'
#' @inheritParams ramp.xds::get_skillset_MY
#'
#' @return the **MY** module skill set, as a list
#'
#' @keywords internal
#' @export
get_skillset_MY.BQm = function(xds_obj, s){
  skills = list(
    delay = FALSE,
    not_implemented = TRUE
  )
  return(skills)
}

#' Run a check before solving
#'
#' @inheritParams ramp.xds::check_MY
#'
#' @return an **`xds`** model object
#' @keywords internal
#' @export
check_MY.BQm = function(xds_obj, s){
  return(xds_obj)
}

#' @title Mosquito bionomics for `BQm` (**MY**)
#' @description Implements [MBionomics] for the BQm model
#' @inheritSection BQm Parameters
#' @inheritSection BQm Survival and Dispersal
#' @inheritParams ramp.xds::MBionomics
#' @seealso [BQm], [MEffectSizes.BQm]
#' @return an **`xds`** model object
#' @keywords internal
#' @export
MBionomics.BQm <- function(t, y, xds_obj, s) {with(xds_obj$MY_obj[[s]],{
  xds_obj$MY_obj[[s]]$f_t        <- F_f(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$q_t        <- F_q(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$g_t        <- F_g(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$sigma_b_t  <- F_sigma_b(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$sigma_q_t  <- F_sigma_q(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$mu         <- F_mu(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$nu         <- F_nu(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$eip        <- F_eip(t, xds_obj, s)

  xds_obj$MY_obj[[s]]$es_f     <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_q     <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_g     <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_sigma_b <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_sigma_q <- rep(1, xds_obj$nPatches)
  return(xds_obj)
})}


#' @title Apply effect sizes for `BQm` (**MY**)
#' @description Implements [MEffectSizes] for the BQm model
#' @inheritSection BQm Parameters
#' @inheritSection BQm Survival and Dispersal
#' @inheritParams ramp.xds::MEffectSizes
#' @seealso [BQm], [MBionomics.BQm]
#' @return an **`xds`** model object
#' @keywords internal
#' @export
MEffectSizes.BQm <- function(t, y, xds_obj, s) {with(xds_obj$MY_obj[[s]],{
  xds_obj$MY_obj[[s]]$f <- es_f*f_t
  xds_obj$MY_obj[[s]]$q <- es_q*q_t
  xds_obj$MY_obj[[s]]$g <- es_g*g_t
  xds_obj$MY_obj[[s]]$sigma_b <- es_sigma_b*sigma_b_t
  xds_obj$MY_obj[[s]]$sigma_q <- es_sigma_q*sigma_q_t

  xds_obj <- change_Omega_b(xds_obj, s)
  xds_obj <- change_Omega_q(xds_obj, s)

  return(xds_obj)
})}

#' @title Infective mosquito biting rate for `BQm`
#' @description Implements [F_fqZ] for the BQm model.
#' @inheritSection BQm Infective Biting
#' @inheritParams ramp.xds::F_fqZ
#' @seealso [BQm], [F_fqM.BQm]
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
F_fqZ.BQm <- function(t, y, xds_obj, s) {
  f = get_f(xds_obj, s)
  q = get_f(xds_obj, s)
  fqZ = f*q*y[xds_obj$MY_obj[[s]]$ix$Bz_ix]
  return(fqZ)
}

#' @title Mosquito biting rate for `BQm`
#' @description Returns `fqM`, the daily mosquito biting density on humans
#' for each patch in the BQm model.
#' @inheritSection BQm Infective Biting
#' @inheritParams ramp.xds::F_fqM
#' @seealso [BQm], [ramp.xds::F_hbr]
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
F_fqM.BQm <- function(t, y, xds_obj, s) {
  f = get_f(xds_obj, s)
  q = get_f(xds_obj, s)
  M = with(xds_obj$MY_obj[[s]]$ix, y[Bu_ix] + y[By_ix] + y[Bz_ix])
  fqM = f*q*M
  return(fqM)
}

#' @title Egg production rate for `BQm`
#' @description Implements [F_eggs] for the BQm model.
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


#' @title Set up `BQm` (**MY**)
#' @description
#'
#' This sets up a `BQm` model object for the **MY** module for the \eqn{s^{th}} vector species
#' + the object is created by [make_MY_obj_BQm]
#' + the initial conditions are set by [setup_MY_inits]
#' + `options` is passed to both to override the defaults
#'
#' @inheritSection BQm State Variables
#' @inheritSection BQm Parameters
#' @seealso [BQm], [make_MY_obj_BQm], [setup_MY_inits.BQm]
#'
#' @note The module object is stored at `xds_obj$MY_obj[[s]]`
#'
#' @inheritParams ramp.xds::setup_MY_obj
#' @return a [list] vector
#' @keywords internal
#' @export
setup_MY_obj.BQm = function(MYname, xds_obj, options=list(), s=1){
  xds_obj$MYname <- "BQm"
  xds_obj$MY_obj[[s]] = make_MY_obj_BQm(xds_obj$nPatches, options)
  xds_obj <- setup_F_circadian("setup", xds_obj, s=s)
  xds_obj$MY_obj[[s]]$skill_set <- get_skillset_MY(xds_obj, s)
  xds_obj <- setup_K_matrix("zero", xds_obj, options=list(which_K = "Kb"), s=s)
  xds_obj <- setup_K_matrix("zero", xds_obj, options=list(which_K = "Kq"), s=s)
  xds_obj <- setup_MY_inits(xds_obj, options=options, s=s)
  return(xds_obj)
}

#' @title Make parameters for `BQm` (**MY**)
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


#' @title Return the parameters as a list
#'
#' @description This method dispatches on the type of `xds_obj$MY_obj[[s]]`.
#' @inheritSection BQm Parameters
#'
#' @inheritParams ramp.xds::get_MY_pars
#'
#' @return a [list]
#' @seealso [BQm], [change_MY_pars.BQm]
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
#' @description For the `BQm` module (**MY**) for the \eqn{s^{th}} species,
#' change parameter values named in `options` leaving other values unchanged.
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
    xds_obj <- MEffectSizes(0, y, xds_obj, s)

    return(xds_obj)
  }))}

#' @title Setup initial values for `BQm` (**MY**)
#' @description Set initial values for all the variables in the `BQm` module
#' (**MY**) for the \eqn{s^{th}} species.
#' + [make_MY_inits_BQm] supplies default values
#' + The values in `options` override the defaults
#' + To inspect variable names, use [get_MY_inits] or
#' see the inline module documentation.
#' @inheritSection BQm State Variables
#' @inheritParams ramp.xds::setup_MY_inits
#' @return a [list]
#' @keywords internal
#' @export
setup_MY_inits.BQm = function(xds_obj, options=list(), s=1){
  xds_obj$MY_obj[[s]]$inits = make_MY_inits_BQm(xds_obj$nPatches, options)
  return(xds_obj)
}

#' @title Make inits for BQm adult mosquito model
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
#' @return a [list]
#' @keywords internal
#' @export
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

#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [setup_MY_ix] for the BQm model.
#' @inheritSection BQm State Variables
#' @inheritParams ramp.xds::setup_MY_ix
#' @return none
#' @importFrom utils tail
#' @keywords internal
#' @export
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

#' @title Return the variables as a list
#' @description This method dispatches on the type of `xds_obj$MY_obj[[s]]`
#' @inheritSection BQm State Variables
#' @inheritParams ramp.xds::get_MY_vars
#' @return a [list]
#' @keywords internal
#' @export
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

#' @title Parse outputs for `BQm` (**MY**)
#' @description Implements [parse_MY_orbits] for the BQm model
#' @inheritSection BQm State Variables
#' @inheritSection BQm Infective Biting
#' @inheritParams ramp.xds::parse_MY_orbits
#' @seealso [BQm], [get_MY_vars.BQm]
#' @return none
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


#' @title Set new MY parameter values
#' @description For the `BQm` module (**MY**) for the \eqn{s^{th}} species,
#' change initial values named in `options` leaving other values unchanged.
#' @inheritParams ramp.xds::change_MY_inits
#' @return an **`xds`** model object
#' @keywords internal
#' @export
change_MY_inits.BQm <- function(xds_obj, options=list(), s=1) {
  with(xds_obj$MY_obj[[s]],
    with(options,{
      xds_obj$MY_obj[[s]]$inits$Bu = Bu
      xds_obj$MY_obj[[s]]$inits$Qu = Qu
      xds_obj$MY_obj[[s]]$inits$By = By
      xds_obj$MY_obj[[s]]$inits$Qy = Qy
      xds_obj$MY_obj[[s]]$inits$Bz = Bz
      xds_obj$MY_obj[[s]]$inits$Qz = Qz
      return(xds_obj)
  }))}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** model object
#' @param s the vector species index
#' @return a [numeric] vector
#' @keywords internal
#' @export
get_f.BQm = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], f_t*es_f)
}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** model object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @keywords internal
#' @export
get_q.BQm = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], q_t*es_q)
}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** model object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @keywords internal
#' @export
get_g.BQm = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], g_t*es_g)
}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** model object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @keywords internal
#' @export
get_sigma.BQm = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], sigma_b_t*es_sigma_b)
}
