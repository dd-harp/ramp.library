# specialized methods for the adult mosquito BQm model

#' @title The `BQm` module for the MY component
#' @description
#' Implements the **MY** component using a BQm (Blood-feeding / egg-laying Queue,
#' mosquito-only) model. This is the uninfected-only version of the BQ model:
#' mosquitoes alternate between a blood-feeding state (**B**) and a gravid /
#' egg-laying state (**Q**), without tracking infection states.
#'
#' @section State Variables:
#' \describe{
#'   \item{`B`}{density of blood-feeding (not gravid) mosquitoes}
#'   \item{`Q`}{density of gravid (egg-laying) mosquitoes}
#' }
#'
#' @section Parameters:
#' \describe{
#'   \item{`f`}{blood feeding rate}
#'   \item{`q`}{human blood fraction}
#'   \item{`g`}{mosquito mortality rate}
#'   \item{`sigma_b`}{emigration rate during blood feeding}
#'   \item{`sigma_q`}{emigration rate during egg laying}
#'   \item{`mu`}{fraction lost during emigration}
#'   \item{`nu`}{oviposition rate (per mosquito)}
#'   \item{`eggsPerBatch`}{eggs laid per oviposition bout}
#' }
#'
#' @section Dynamics:
#' \deqn{
#' \begin{array}{rl}
#' dB/dt &= \Lambda + \nu Q - f B - \Omega_b \cdot B \\
#' dQ/dt &= f B - \nu Q - \Omega_q \cdot Q \\
#' \end{array}
#' }
#'
#' @name BQm
#' @rdname BQm
NULL

#' @title Compute derivatives for `BQm` (**MY**)
#' @description Implements [dMYdt] for the BQm ODE model.
#' @inheritParams ramp.xds::dMYdt
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

#' @title The **BQm** module skill set
#'
#' @description The **MY** skill set is a list of
#' a module's capabilities:
#'
#' + `demography` is
#'
#' @inheritParams ramp.xds::skill_set_MY
#'
#' @return *MY* module skill set, as a list
#'
#' @keywords internal
#' @export
skill_set_MY.BQm = function(MYname){
  return(list())
}

#' Run a check before solving
#'
#' @inheritParams ramp.xds::check_MY
#'
#' @return an **`xds`** object
#' @keywords internal
#' @export
check_MY.BQm = function(xds_obj, s){
  return(xds_obj)
}

#' @title Mosquito bionomics for `BQm` (**MY**)
#' @description Implements [MBionomics] for the BQm model
#' @inheritParams ramp.xds::MBionomics
#' @return an **`xds`** object
#' @keywords internal
#' @export
MBionomics.BQm <- function(t, y, xds_obj, s) {with(xds_obj$MY_obj[[s]],{
  xds_obj$MY_obj[[s]]$es_g       <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_sigma_b <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_sigma_q <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_f       <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_q       <- rep(1, nPatches)
  return(xds_obj)
})}


#' @title Apply effect sizes for `BQm` (**MY**)
#' @description Implements [MEffectSizes] for the BQm model
#' @inheritParams ramp.xds::MEffectSizes
#' @return an **`xds`** object
#' @keywords internal
#' @export
MEffectSizes.BQm <- function(t, y, xds_obj, s) {with(xds_obj$MY_obj[[s]],{
  xds_obj$MY_obj[[s]]$f <- es_f*f_t
  xds_obj$MY_obj[[s]]$q <- es_q*q_t
  g <- es_g*g_t
  sigma_b <- es_sigma_b*sigma_b_t
  sigma_q <- es_sigma_q*sigma_q_t
  xds_obj$MY_obj[[s]]$g <- g
  xds_obj$MY_obj[[s]]$sigma_b <- sigma_b
  xds_obj$MY_obj[[s]]$sigma_q <- sigma_q
  xds_obj$MY_obj[[s]]$Omega_b = make_Omega_xde(g, sigma_b, mu, calKb)
  xds_obj$MY_obj[[s]]$Omega_q = make_Omega_xde(g, sigma_q, mu, calKq)
  return(xds_obj)
})}

#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqZ] for the BQm model.
#' @inheritParams ramp.xds::F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
F_fqZ.BQm <- function(t, y, xds_obj, s) {
  f = get_f(xds_obj, s)
  q = get_f(xds_obj, s)
  fqZ = f*q*y[xds_obj$MY_obj[[s]]$ix$Bz_ix]
  return(fqZ)
}

#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqM] for the BQm model.
#' @inheritParams ramp.xds::F_fqM
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

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the BQm model.
#' @inheritParams ramp.xds::F_eggs
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


#' @title Setup MY_obj for the BQm model
#' @description Implements [setup_MY_obj] for the RM model
#' @inheritParams ramp.xds::setup_MY_obj
#' @return a [list] vector
#' @keywords internal
#' @export
setup_MY_obj.BQm = function(MYname, xds_obj, s, options=list()){
  xds_obj$MY_obj[[s]] = make_MY_obj_BQm(xds_obj$nPatches, options)
  return(xds_obj)
}

#' @title Make parameters for RM ODE adult mosquito model
#' @param nPatches is the number of patches, an integer
#' @param options a [list] of values that overwrites the defaults
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
make_MY_obj_BQm = function(nPatches, options=list(), eip=12,
                          g=1/12, sigma_b=1/8, sigma_q=1/8, mu=0, f=0.5, q=0.95,
                          nu=1, eggsPerBatch=60){

  with(options,{
    MY_obj <- list()
    class(MY_obj) <- "BQm"

    MY_obj$nPatches <- nPatches

    MY_obj$g_t          <- checkIt(g, nPatches)
    MY_obj$es_g         <- rep(1, nPatches)
    MY_obj$sigma_b_t    <- checkIt(sigma_b, nPatches)
    MY_obj$es_sigma_b   <- rep(1, nPatches)
    MY_obj$sigma_q_t    <- checkIt(sigma_q, nPatches)
    MY_obj$es_sigma_q   <- rep(1, nPatches)
    MY_obj$mu           <- checkIt(mu, nPatches)
    MY_obj$f_t          <- checkIt(f, nPatches)
    MY_obj$es_f         <- rep(1, nPatches)
    MY_obj$q_t          <- checkIt(q, nPatches)
    MY_obj$es_q         <- rep(1, nPatches)
    MY_obj$nu           <- checkIt(nu, nPatches)
    MY_obj$eggsPerBatch <- eggsPerBatch

    calK <- diag(nPatches)

    MY_obj$calKb <- calK
    MY_obj$calKq <- calK

    Omega_par <- list()
    class(Omega_par) <- "static"
    MY_obj$Omega_par <- Omega_par
    MY_obj$Omega_b <- with(MY_obj, make_Omega_xde(g, sigma_b, mu, calK))
    MY_obj$Omega_q <- with(MY_obj, make_Omega_xde(g, sigma_q, mu, calK))
    base <- 'BQm'
    class(base) <- 'BQm'
    MY_obj$baseline <- base

    return(MY_obj)
})}


#' @title Return the parameters as a list
#'
#' @description This method dispatches on the type of `xds_obj$MY_obj[[s]]`.
#'
#' @inheritParams ramp.xds::get_MY_pars
#'
#' @return a [list]
#' @keywords internal
#' @export
get_MY_pars.BQm <- function(xds_obj, s=1) {
  with(xds_obj$MY_obj[[s]], list(
    f=f_t, q=q_t, g=g_t, sigma=sigma_t, eip=eip, mu=mu_t,
    nu=nu_t, eggsPerBatch=eggsPerBatch, calK=calK
  ))
}


#' @title Return the parameters as a list
#' @description This method dispatches on the type of `xds_obj$MY_obj[[s]]`.
#' @inheritParams ramp.xds::change_MY_pars
#' @return an **`xds`** object
#' @keywords internal
#' @export
change_MY_pars.BQm <- function(xds_obj, s=1, options=list()) {
  nHabitats <- xds_obj$nHabitats
  with(xds_obj$MY_obj[[s]], with(options,{
    xds_obj$MY_obj[[s]]$f_t = f
    xds_obj$MY_obj[[s]]$q_t = q
    xds_obj$MY_obj[[s]]$g_t = g
    xds_obj$MY_obj[[s]]$sigma_b_t = sigma_b
    xds_obj$MY_obj[[s]]$sigma_q_t = sigma_q
    xds_obj$MY_obj[[s]]$eip = eip
    xds_obj$MY_obj[[s]]$mu = mu
    xds_obj$MY_obj[[s]]$nu_t = nu
    xds_obj$MY_obj[[s]]$eggsPerBatch = eggsPerBatch
    return(xds_obj)
  }))}

#' @title Setup initial values for the BQm model
#' @description Implements [setup_MY_inits] for the RM model
#' @inheritParams ramp.xds::setup_MY_inits
#' @return a [list]
#' @keywords internal
#' @export
setup_MY_inits.BQm = function(xds_obj, s, options=list()){
  xds_obj$MY_obj[[s]]$inits = make_MY_inits_BQm(xds_obj$nPatches, options)
  return(xds_obj)
}

#' @title Make inits for BQm adult mosquito model
#' @param nPatches the number of patches in the model
#' @param options a [list] of values that overwrites the defaults
#' @param Bu0 total uninfected, not gravid mosquito density at each patch
#' @param Qu0 total uninfected, gravid mosquito density at each patch
#' @param By0 infected, not gravid mosquito density at each patch
#' @param Qy0 infected, gravid mosquito density at each patch
#' @param Bz0 infectious, not gravid mosquito density at each patch
#' @param Qz0 infectious, gravid mosquito density at each patch
#' @return a [list]
#' @keywords internal
#' @export
make_MY_inits_BQm = function(nPatches, options = list(),
                             Bu0=5, Qu0=1, By0=1, Qy0=1, Bz0=1, Qz0=1){
  with(options,{
    Bu = checkIt(Bu0, nPatches)
    Qu = checkIt(Qu0, nPatches)
    By = checkIt(By0, nPatches)
    Qy = checkIt(Qy0, nPatches)
    Bz = checkIt(Bz0, nPatches)
    Qz = checkIt(Qz0, nPatches)

    return(list(Bu=Bu, Qu=Qu, By=By, Qy=Qy, Bz=Bz, Qz=Qz))
  })}

#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [setup_MY_ix] for the BQm model.
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

#' @title parse the output of deSolve and return variables for the BQm model
#' @description Implements [parse_MY_orbits] for the BQm model
#' @inheritParams ramp.xds::parse_MY_orbits
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
#' @description This method dispatches on the type of `xds_obj$MY_obj[[s]]`.
#' @inheritParams ramp.xds::change_MY_inits
#' @return an `xds` object
#' @keywords internal
#' @export
change_MY_inits.BQm <- function(xds_obj, s=1, options=list()) {
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
#' @param xds_obj an **`xds`** object
#' @param s the vector species index
#' @return a [numeric] vector
#' @keywords internal
#' @export
get_f.BQm = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], f_t*es_f)
}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @keywords internal
#' @export
get_q.BQm = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], q_t*es_q)
}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @keywords internal
#' @export
get_g.BQm = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], g_t*es_g)
}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @keywords internal
#' @export
get_sigma.BQm = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], sigma_b_t*es_sigma_b)
}

