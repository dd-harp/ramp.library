# specialized methods for the adult mosquito SEI model

#' @title The `SEI` module for the MY component
#' @description
#' Implements the **MY** component using a Susceptible-Exposed-Infectious
#' (SEI) compartmental model of adult mosquito infection dynamics.
#'
#' @section State Variables:
#' \describe{
#'   \item{`M`}{density of adult mosquitoes}
#'   \item{`Y`}{density of infected (exposed) adult mosquitoes}
#'   \item{`Z`}{density of infectious adult mosquitoes}
#' }
#'
#' @section Parameters:
#' \describe{
#'   \item{`f`}{blood feeding rate}
#'   \item{`q`}{human blood fraction}
#'   \item{`g`}{mosquito mortality rate}
#'   \item{`sigma`}{patch emigration rate}
#'   \item{`mu`}{emigration-related loss}
#'   \item{`eip`}{extrinsic incubation period (\eqn{\tau})}
#'   \item{`nu`}{oviposition rate (per mosquito)}
#'   \item{`eggsPerBatch`}{eggs laid per oviposition bout}
#' }
#'
#' @section Dynamics:
#' \deqn{
#' \begin{array}{rl}
#' dM/dt &= \Lambda - \Omega \cdot M \\
#' dY/dt &= fq\kappa(M-Y) - \Omega \cdot Y \\
#' dZ/dt &= (Y-Z)/\tau - \Omega \cdot Z \\
#' \end{array}
#' }
#' where \eqn{\Omega} is the demographic matrix encoding mortality and dispersal,
#' and \eqn{\kappa} is the net infectiousness of humans.
#'
#' @name SEI
#' @rdname SEI
NULL

#' @title The skill set
#'
#' @inheritParams ramp.xds::setup_skillset_MY
#'
#' @return the **`xds`** object
#'
#' @keywords internal
#' @export
setup_skillset_MY.SEI = function(xds_obj, s){
  skills = list(
    not_implemented = TRUE
  )
  xds_obj$MY_obj[[s]]$skill_set = skills
  return(xds_obj)
}

#' Run a check before solving
#'
#' @inheritParams ramp.xds::check_MY
#'
#' @return an **`xds`** object
#' @keywords internal
#' @export
check_MY.SEI = function(xds_obj, s){
  return(xds_obj)
}

#' @title Compute derivatives for `SEI` (**MY**)
#' @description Implements [dMYdt] for the SEI ODE model.
#'
#' The dynamics of adult mosquitoes:
#' \deqn{dM/dt = \Lambda - \Omega \cdot M}
#'
#' The density of infected but not infectious mosquitoes:
#' \deqn{dY/dt = fq\kappa(M-Y) - \Omega \cdot Y}
#'
#' The density of infectious mosquitoes:
#' \deqn{dZ/dt = (Y-Z)/\tau - \Omega \cdot Z}
#'
#' @inheritParams ramp.xds::dMYdt
#' @return a [numeric] vector
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


#' @title Derivatives for adult mosquitoes
#' @description Implements [Update_MYt] for the SEI model.
#' @inheritParams ramp.xds::Update_MYt
#' @return a [numeric] vector
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

#' @title Setup MY_obj for the SEI model
#' @description Implements [setup_MY_obj] for the SEI model
#' @inheritParams ramp.xds::setup_MY_obj
#' @return a [list] vector
#' @keywords internal
#' @export
setup_MY_obj.SEI = function(MYname, xds_obj, s, options=list()){
  xds_obj$MYname <- "SEI"
  MY_obj <- make_MY_obj_SEI(xds_obj$nPatches, options)
  class(MY_obj) <- c("SEI", paste("SEI_", xds_obj$xds, sep=""))
  xds_obj$MY_obj[[s]] <- MY_obj
  xds_obj <- setup_MY_inits(xds_obj, s, options)
  xds_obj <- setup_F_circadian("setup", xds_obj, s=s)
  xds_obj <- setup_Omega_obj("xde", "static", xds_obj, s=s)
  xds_obj <- setup_Upsilon_obj(xds_obj, s=s, FALSE)
  xds_obj <- setup_K_matrix("zero", xds_obj, s=s)
  return(xds_obj)
}


#' @title Return the parameters as a list
#' @description This method dispatches on the type of `xds_obj$MY_obj[[s]]`.
#' @inheritParams ramp.xds::change_MY_pars
#' @return an **`xds`** object
#' @keywords internal
#' @export
change_MY_pars.SEI <- function(xds_obj, s=1, options=list()) {
  nHabitats <- xds_obj$nHabitats
  with(xds_obj$MY_obj[[s]], with(options,{
    xds_obj$MY_obj[[s]]$f_t = f
    xds_obj$MY_obj[[s]]$q_t = q
    xds_obj$MY_obj[[s]]$g_t = g
    xds_obj$MY_obj[[s]]$sigma_t = sigma
    xds_obj$MY_obj[[s]]$eip_t = eip
    xds_obj$MY_obj[[s]]$mu_t = mu
    xds_obj$MY_obj[[s]]$nu_t = nu
    xds_obj$MY_obj[[s]]$eggsPerBatch = eggsPerBatch
    return(xds_obj)
  }))}

#' @title Return the parameters as a list
#'
#' @description This method dispatches on the type of `xds_obj$MY_obj[[s]]`.
#'
#' @inheritParams ramp.xds::get_MY_pars
#'
#' @return a [list]
#' @keywords internal
#' @export
get_MY_pars.SEI <- function(xds_obj, s=1) {
  with(xds_obj$MY_obj[[s]], list(
    f=f_t, q=q_t, g=g_t, sigma=sigma_t, eip=eip, mu=mu_t,
    nu=nu_t, eggsPerBatch=eggsPerBatch, K_matrix=K_matrix
  ))
}

#' @title Set new MY parameter values
#' @description This method dispatches on the type of `xds_obj$MY_obj[[s]]`.
#' @inheritParams ramp.xds::change_MY_inits
#' @return an `xds` object
#' @keywords internal
#' @export
change_MY_inits.SEI <- function(xds_obj, s=1, options=list()) {
  with(xds_obj$MY_obj[[s]]$inits,
    with(options,{
      xds_obj$MY_obj[[s]]$inits$M = M
      xds_obj$MY_obj[[s]]$inits$Y = Y
      xds_obj$MY_obj[[s]]$inits$Z = Z
      return(xds_obj)
}))}

#' @title The net blood feeding rate of the infective mosquito population in a patch
#' @description Implements [F_fqZ] for the SEI model.
#' @inheritParams ramp.xds::F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
F_fqZ.SEI <- function(t, y, xds_obj, s) {
  f = get_f(xds_obj, s)
  q = get_q(xds_obj, s)
  Z = y[xds_obj$MY_obj[[s]]$ix$Z_ix]
  return(f*q*Z)
}

#' @title The net blood feeding rate of the infective mosquito population in a patch
#' @description Implements [F_fqM] for the SEI model.
#' @inheritParams ramp.xds::F_fqM
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
F_fqM.SEI <- function(t, y, xds_obj, s) {
  f = get_f(xds_obj, s)
  q = get_q(xds_obj, s)
  M = y[xds_obj$MY_obj[[s]]$ix$M_ix]
  return(f*q*M)
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the SEI model.
#' @inheritParams ramp.xds::F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
F_eggs.SEI <- function(t, y, xds_obj, s) {
  M <- y[xds_obj$MY_obj[[s]]$ix$M_ix]
  with(xds_obj$MY_obj[[s]],{
    return(M*nu*eggsPerBatch)
  })
}

#' @title Return the variables as a list
#' @description This method dispatches on the type of `xds_obj$MY_obj[[s]]`
#' @inheritParams ramp.xds::get_MY_vars
#' @return a [list]
#' @keywords internal
#' @export
get_MY_vars.SEI <- function(y, xds_obj, s){
  with(xds_obj$MY_obj[[s]]$ix,
       return(list(
         M = y[M_ix],
         Y = y[Y_ix],
         Z = y[Z_ix]
       )))
}


#' @title Setup initial values for the SEI model
#' @description Implements [setup_MY_inits] for the SEI model
#' @inheritParams ramp.xds::setup_MY_inits
#' @return a [list]
#' @keywords internal
#' @export
setup_MY_inits.SEI = function(xds_obj, s, options=list()){
  xds_obj$MY_obj[[s]]$inits = make_MY_inits_SEI(xds_obj$nPatches, options)
  return(xds_obj)
}


#' @title Make inits for SEI adult mosquito model
#' @param nPatches the number of patches in the model
#' @param options a [list] of values that overwrites the defaults
#' @param M total mosquito density at each patch
#' @param Y total infected mosquito density at each patch
#' @param Z infectious mosquito density at each patch
#' @return a [list]
#' @keywords internal
#' @export
make_MY_inits_SEI = function(nPatches, options = list(),
                               M=5, Y=1, Z=0){
  with(options,{
    M = checkIt(M, nPatches)
    Y = checkIt(Y, nPatches)
    Z = checkIt(Z, nPatches)
    return(list(M=M, Y=Y, Z=Z))
  })
}



#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [setup_MY_ix] for the SEI model.
#' @inheritParams ramp.xds::setup_MY_ix
#' @return a [list]
#' @importFrom utils tail
#' @keywords internal
#' @export
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


#' @title parse the output of deSolve and return variables for the SEI model
#' @description Implements [parse_MY_orbits] for the SEI model
#' @inheritParams ramp.xds::parse_MY_orbits
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

#' @title Return the parameters as a list
#' @description This method dispatches on the type of `xds_obj$MY_obj[[s]]`.
#' @inheritParams change_MY_pars
#' @return an **`xds`** object
#' @keywords internal
#' @export
change_MY_pars.SEI <- function(xds_obj, s=1, options=list()) {
  nHabitats <- xds_obj$nHabitats
  with(xds_obj$MY_obj[[s]], with(options,{
    xds_obj$MY_obj[[s]]$f_obj$f = checkIt(f, nPatches)
    xds_obj$MY_obj[[s]]$q_obj$q = checkIt(q, nPatches)
    xds_obj$MY_obj[[s]]$g_obj$g = checkIt(g, nPatches)
    xds_obj$MY_obj[[s]]$sigma_obj$sigma = checkIt(sigma, nPatches)
    xds_obj$MY_obj[[s]]$eip_obj$eip = checkIt(eip, nPatches)
    xds_obj$MY_obj[[s]]$mu_obj$mu = checkIt(mu, nPatches)
    xds_obj$MY_obj[[s]]$nu_obj$nu = checkIt(nu, nPatches)
    xds_obj$MY_obj[[s]]$eggsPerBatch = checkIt(eggsPerBatch, nPatches)

    xds_obj$MY_obj[[s]]$f_t = checkIt(f, nPatches)
    xds_obj$MY_obj[[s]]$q_t = checkIt(q, nPatches)
    xds_obj$MY_obj[[s]]$g_t = checkIt(g, nPatches)
    xds_obj$MY_obj[[s]]$sigma_t = checkIt(sigma, nPatches)
    xds_obj$MY_obj[[s]]$eip_t = checkIt(eip, nPatches)
    xds_obj$MY_obj[[s]]$mu_t = checkIt(mu, nPatches)
    xds_obj$MY_obj[[s]]$nu_t = checkIt(nu, nPatches)
    xds_obj$MY_obj[[s]]$eggsPerBatch = eggsPerBatch
    return(xds_obj)
  }))}


#' @title Make parameters for SEI ODE adult mosquito model
#' @param nPatches is the number of patches, an integer
#' @param options a [list] of values that overwrites the defaults
#' @param eip the extrinsic incubation period
#' @param g mosquito mortality rate
#' @param sigma emigration rate
#' @param mu emigration loss
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @return a [list]
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


#' @title Mosquito bionomics for `SEI` (**MY**)
#' @description Implements [MBionomics] for the SEI model
#' @inheritParams ramp.xds::MBionomics
#' @return an **`xds`** object
#' @keywords internal
#' @export
MBionomics.SEI <- function(t, y, xds_obj, s){with(xds_obj$MY_obj[[s]],{
  # Baseline parameters
  xds_obj$MY_obj[[s]]$f_t      <- F_f(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$q_t      <- F_q(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$g_t      <- F_g(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$sigma_t  <- F_sigma(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$mu       <- F_mu(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$nu       <- F_nu(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$eip      <- F_eip(t, xds_obj, s)
  xds_obj                      <- F_K_matrix(t, xds_obj, s)
  xds_obj$MY_obj[[s]]$eggsPerBatch <- eggsPerBatch

  # Reset Effect Sizes
  xds_obj$MY_obj[[s]]$es_f     < rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_q     <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_g     <- rep(1, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$es_sigma <- rep(1, xds_obj$nPatches)
  return(xds_obj)
})}

#' @title Apply effect sizes for `SEI` (**MY**)
#' @description Implements [MEffectSizes] for the SEI model
#' @inheritParams ramp.xds::MEffectSizes
#' @return an **`xds`** object
#' @keywords internal
#' @export
MEffectSizes.SEI <- function(t, y, xds_obj, s) {
  with(xds_obj$MY_obj[[s]],{
    xds_obj$MY_obj[[s]]$f <- es_f*f_t
    xds_obj$MY_obj[[s]]$q <- es_q*q_t
    xds_obj$MY_obj[[s]]$g <- es_g*g_t
    xds_obj$MY_obj[[s]]$sigma <- es_sigma*sigma_t
    xds_obj <- update_Omega(xds_obj, s)
    return(xds_obj)
})}


#' @title Get the feeding rate
#' @param xds_obj an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @keywords internal
#' @export
get_f.SEI = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], f_t*es_f)
}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @keywords internal
#' @export
get_q.SEI = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], q_t*es_q)
}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @keywords internal
#' @export
get_g.SEI = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], g_t*es_g)
}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @keywords internal
#' @export
get_sigma.SEI = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], sigma_t*es_sigma)
}

#' @title Steady States: MY-SEI
#' @description This method dispatches on the type of `MY_obj`.
#' @inheritParams ramp.xds::steady_state_MY
#' @return none
#' @keywords internal
#' @export
steady_state_MY.SEI = function(Lambda, kappa, xds_obj, s=1){
  with(xds_obj$MY_obj[[s]],{
    Omega_inv <- solve(Omega)
    Z_eq <- as.vector(solve(diag(f*q*kappa) + Omega) %*% diag(f*q*kappa) %*% M_eq)
    return(c(M=M_eq, Z=Z_eq))
})}
