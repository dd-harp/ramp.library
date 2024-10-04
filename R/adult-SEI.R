# specialized methods for the adult mosquito SEI model

#' @title \eqn{\cal MYZ} Component Derivatives for the `SEI` Mosquito Model
#' @description Implements [dMYZdt] for the SEI ODE model.
#' @details
#' The dynamics of adult mosquitoes:
#' \deqn{\frac{dM}{dt} = \Lambda - \Omega \cdot M}
#' The density of infected but not infectious mosquitoes:
#' \deqn{\frac{dY}{dt} = fq\kappa(M-Y) - \Omega \cdot Y}
#' The density of infectious mosquitoes:
#' \deqn{\frac{dZ}{dt} = (Y-Z)/\tau - \Omega \cdot Z}
#' @inheritParams ramp.xds::dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.SEI <- function(t, y, pars, s) {
  Lambda = pars$Lambda[[s]]
  kappa = pars$kappa[[s]]

  with(list_MYZvars(y, pars, s),{
    with(pars$MYZpar[[s]],{

      dM <- Lambda - (Omega %*% M)
      dY <- f*q*kappa*(M-Y) - Omega %*% Y
      dZ <- (Y-Z)/eip - (Omega %*% Z)

      return(c(dM, dY, dZ))
    })
  })
}

#' @title Steady States: MYZ-SEI
#' @description This method dispatches on the type of `MYZpar`.
#' @inheritParams ramp.xds::xde_steady_state_MYZ
#' @return none
#' @export
xde_steady_state_MYZ.SEI = function(Lambda, kappa, MYZpar){with(MYZpar,{
  Omega_inv <- solve(Omega)
  Z_eq <- as.vector(solve(diag(f*q*kappa) + Omega) %*% diag(f*q*kappa) %*% M_eq)
  return(c(M=M_eq, Z=Z_eq))
})}

#' @title Derivatives for adult mosquitoes
#' @description Implements [Update_MYZt] for the SEI model.
#' @inheritParams ramp.xds::Update_MYZt
#' @return a [numeric] vector
#' @export
Update_MYZt.SEI <- function(t, y, pars, s) {
  Lambda = pars$Lambda[[s]]
  kappa = pars$kappa[[s]]
  D = pars$MYZday

  with(list_MYZvars(y, pars, s),{
    with(pars$MYZpar[[s]],{
      Mt <- ccc*Lambda*D + Omega %*% M
      Zt <- Omega%*%((1-exp(-f*q*kappa*pars$MYZday))*(M-Z)) + (Omega %*% Z)
      Zt <- Zt + (1-exp(-ccc))*(1-exp(-f*q*kappa*pars$MYZday))*Lambda*ccc*D
      return(list(M=unname(Mt), Z=unname(Zt)))
    })
  })
}

#' @title Setup MYZpar for the SEI model
#' @description Implements [setup_MYZpar] for the SEI model
#' @inheritParams ramp.xds::setup_MYZpar
#' @return a [list] vector
#' @export
setup_MYZpar.SEI = function(MYZname, pars, s, MYZopts=list()){
  MYZpar <- make_MYZpar_SEI(pars$nPatches, MYZopts)
  class(MYZpar) <- 'SEI'
  pars$MYZpar[[s]] <- MYZpar
  return(pars)
}


#' @title Return the parameters as a list
#' @description This method dispatches on the type of `pars$MYZpar[[s]]`.
#' @inheritParams ramp.xds::set_MYZpars
#' @return an **`xds`** object
#' @export
set_MYZpars.SEI <- function(pars, s=1, MYZopts=list()) {
  nHabitats <- pars$nHabitats
  with(pars$MYZpar[[s]], with(MYZopts,{
    pars$MYZpar[[s]]$f_t = f
    pars$MYZpar[[s]]$q_t = q
    pars$MYZpar[[s]]$g_t = g
    pars$MYZpar[[s]]$sigma_t = sigma
    pars$MYZpar[[s]]$eip_t = eip
    pars$MYZpar[[s]]$mu_t = mu
    pars$MYZpar[[s]]$nu_t = nu
    pars$MYZpar[[s]]$eggsPerBatch = eggsPerBatch
    return(pars)
  }))}

#' @title Return the parameters as a list
#' @description This method dispatches on the type of `pars$MYZpar[[s]]`.
#' @param pars an **`xds`** object
#' @param s the vector species index
#' @return a [list]
#' @export
get_MYZpars.SEI <- function(pars, s=1) {
  with(pars$MYZpar[[s]], list(
    f=f_t, q=q_t, g=g_t, sigma=sigma_t, eip=eip, mu=mu_t,
    nu=nu_t, eggsPerBatch=eggsPerBatch, calK=calK
  ))
}
#' @title Set new MYZ parameter values
#' @description This method dispatches on the type of `pars$MYZpar[[s]]`.
#' @inheritParams ramp.xds::set_MYZinits
#' @return an `xds` object
#' @export
set_MYZinits.SEI <- function(pars, s=1, MYZopts=list()) {
  with(pars$MYZpar[[s]], with(MYZopts,{
    pars$MYZinits[[s]]$M = M
    pars$MYZinits[[s]]$Y = Y
    pars$MYZinits[[s]]$Z = Z
    return(pars)
  }))}

#' @title The net blood feeding rate of the infective mosquito population in a patch
#' @description Implements [F_fqZ] for the SEI model.
#' @inheritParams ramp.xds::F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ.SEI <- function(t, y, pars, s) {
  f = get_f(pars, s)
  q = get_q(pars, s)
  Z = y[pars$ix$MYZ[[s]]$Z_ix]
  return(f*q*Z)
}

#' @title The net blood feeding rate of the infective mosquito population in a patch
#' @description Implements [F_fqM] for the SEI model.
#' @inheritParams ramp.xds::F_fqM
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqM.SEI <- function(t, y, pars, s) {
  f = get_f(pars, s)
  q = get_q(pars, s)
  M = y[pars$ix$MYZ[[s]]$M_ix]
  return(f*q*M)
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the SEI model.
#' @inheritParams ramp.xds::F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.SEI <- function(t, y, pars, s) {
  M <- y[pars$ix$MYZ[[s]]$M_ix]
  with(pars$MYZpar[[s]],{
    return(M*nu*eggsPerBatch)
  })
}

#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$MYZpar[[s]]`
#' @inheritParams ramp.xds::list_MYZvars
#' @return a [list]
#' @export
list_MYZvars.SEI <- function(y, pars, s){
  with(pars$ix$MYZ[[s]],
       return(list(
         M = y[M_ix],
         Y = y[Y_ix],
         Z = y[Z_ix]
       )))
}

#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$MYZpar[[s]]`
#' @inheritParams ramp.xds::put_MYZvars
#' @return a [list]
#' @export
put_MYZvars.SEI <- function(MYZvars, y, pars, s){
  with(pars$ix$MYZ[[s]],
       with(MYZvars,{
         y[M_ix] = M
         y[Y_ix] = Y
         y[Z_ix] = Z
         return(y)
       }))
}





#' @title Setup initial values for the SEI model
#' @description Implements [setup_MYZinits] for the SEI model
#' @inheritParams ramp.xds::setup_MYZinits
#' @return a [list]
#' @export
setup_MYZinits.SEI = function(pars, s, MYZopts=list()){
  pars$MYZinits[[s]] = with(pars$MYZpar[[s]], make_MYZinits_SEI(nPatches, MYZopts))
  return(pars)
}


#' @title Make inits for SEI adult mosquito model
#' @param nPatches the number of patches in the model
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param M total mosquito density at each patch
#' @param Y total infected mosquito density at each patch
#' @param Z infectious mosquito density at each patch
#' @return a [list]
#' @export
make_MYZinits_SEI = function(nPatches, MYZopts = list(),
                               M=5, Y=1, Z=0){
  with(MYZopts,{
    M = checkIt(M, nPatches)
    Y = checkIt(Y, nPatches)
    Z = checkIt(Z, nPatches)
    return(list(M=M, Y=Y, Z=Z))
  })
}



#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [setup_MYZix] for the SEI model.
#' @inheritParams ramp.xds::setup_MYZix
#' @return a [list]
#' @importFrom utils tail
#' @export
setup_MYZix.SEI <- function(pars, s) {with(pars,{

  M_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(M_ix, 1)

  Y_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Y_ix, 1)

  Z_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Z_ix, 1)

  pars$max_ix = max_ix
  pars$ix$MYZ[[s]] = list(M_ix=M_ix, Y_ix=Y_ix, Z_ix=Z_ix)
  return(pars)
})}


#' @title Parse the output of deSolve and return variables for the SEI model
#' @description Implements [parse_MYZorbits] for the SEI model
#' @inheritParams ramp.xds::parse_MYZorbits
#' @return a [list]
#' @export
parse_MYZorbits.SEI <- function(outputs, pars, s) {with(pars$ix$MYZ[[s]],{
  M = outputs[,M_ix]
  Y = outputs[,Y_ix]
  Z = outputs[,Z_ix]
  f = get_ft(pars,s)
  q = get_ft(pars,s)
  y = Y/M
  z = Z/M
  return(list(M=M, Z=Z, Y=Y, y=y, z=z, fqZ=f*q*Z, fqM=f*q*M))
})}

#' @title Return initial values as a vector
#' @description Implements [get_MYZinits] for the SEI model.
#' @inheritParams ramp.xds::get_MYZinits
#' @return [numeric]
#' @export
get_MYZinits.SEI <- function(pars, s) {
  pars$MYZinits[[s]]
}


#' @title Return the parameters as a list
#' @description This method dispatches on the type of `pars$MYZpar[[s]]`.
#' @inheritParams set_MYZpars
#' @return an **`xds`** object
#' @export
set_MYZpars.SEI <- function(pars, s=1, MYZopts=list()) {
  nHabitats <- pars$nHabitats
  with(pars$MYZpar[[s]], with(MYZopts,{
    pars$MYZpar[[s]]$f_t = f
    pars$MYZpar[[s]]$q_t = q
    pars$MYZpar[[s]]$g_t = g
    pars$MYZpar[[s]]$sigma_t = sigma
    pars$MYZpar[[s]]$eip_t = eip
    pars$MYZpar[[s]]$mu_t = mu
    pars$MYZpar[[s]]$nu_t = nu
    pars$MYZpar[[s]]$eggsPerBatch = eggsPerBatch
    return(pars)
  }))}
#' @title Make inits for SEI adult mosquito model
#' @inheritParams ramp.xds::update_MYZinits
#' @return a [list]
#' @export
update_MYZinits.SEI <- function(pars, y0, s) {with(pars$ix$MYZ[[s]],{
  M = y0[M_ix]
  Y = y0[Y_ix]
  Z = y0[Z_ix]
  pars$MYZinits[[s]] = setup_MYZinits_SEI(pars$nPatches, list(), M=M, Y=Y, Z=Z)
  return(pars)
})}


#' @title Make parameters for SEI ODE adult mosquito model
#' @param nPatches is the number of patches, an integer
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param eip the extrinsic incubation period
#' @param g mosquito mortality rate
#' @param sigma emigration rate
#' @param mu emigration loss
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @return a [list]
#' @export
make_MYZpar_SEI = function(nPatches, MYZopts=list(), eip =12,
                             g=1/12,  sigma=1/8,  mu=0,
                             f=0.3,  q=0.95,
                             nu=1,  eggsPerBatch=60){
  with(MYZopts,{

    MYZpar <- list()

    MYZpar$nPatches <- nPatches

    eip_par <- list()
    class(eip_par) <- 'static'
    eip_par$eip = eip
    MYZpar$eip_par <- eip_par
    MYZpar$eip=eip

    f=checkIt(f, nPatches)
    MYZpar$f_par <- list()
    class(MYZpar$f_par) <- "static"
    MYZpar$f_par$f = f
    MYZpar$f_t = f
    MYZpar$es_f = 1

    q=checkIt(q, nPatches)
    MYZpar$q_par <- list()
    class(MYZpar$q_par) <- "static"
    MYZpar$q_par$q = q
    MYZpar$q_t = q
    MYZpar$es_q = 1

    g=checkIt(g, nPatches)
    MYZpar$g_par <- list()
    class(MYZpar$g_par) <- "static"
    MYZpar$g_par$g = g
    MYZpar$g_t = g
    MYZpar$es_g = 1

    mu=checkIt(mu, nPatches)
    MYZpar$mu_par <- list()
    class(MYZpar$mu_par) <- "static"
    MYZpar$mu_par$mu = mu
    MYZpar$mu = mu

    sigma=checkIt(sigma, nPatches)
    MYZpar$sigma_par <- list()
    class(MYZpar$sigma_par) <- "static"
    MYZpar$sigma_par$sigma = sigma
    MYZpar$sigma_t = sigma
    MYZpar$es_sigma = 1

    nu=checkIt(nu, nPatches)
    MYZpar$nu_par <- list()
    class(MYZpar$nu_par) <- "static"
    MYZpar$nu_par$nu = nu
    MYZpar$nu=nu

    calK = diag(nPatches)
    MYZpar$calK_par <- list()
    class(MYZpar$calK_par) <- "static"
    MYZpar$calK_par$calK = calK
    MYZpar$calK=calK

    Omega <- diag(g, nPatches)
    MYZpar$Omega <- Omega
    MYZpar$Upsilon <- expm::expm(-Omega*eip)
    MYZpar$nPatches <- nPatches

    Omega <- diag(g, nPatches)
    MYZpar$Omega <- Omega
    MYZpar$Upsilon <- expm::expm(-Omega*eip)
    MYZpar$nPatches <- nPatches

    MYZpar$eggsPerBatch <- eggsPerBatch

    MYZpar$baseline <- MYZpar
    class(MYZpar$baseline) <- 'SEI'

    return(MYZpar)
})}


#' @title Set mosquito bionomics to baseline
#' @description Implements [MBaseline] for models with no forcing on the baseline
#' @inheritParams ramp.xds::MBaseline
#' @return the model as a [list]
#' @export
MBaseline.SEI <- function(t, y, pars, s){with(pars$MYZpar[[s]],{
  # Baseline parameters
  pars$MYZpar[[s]]$f_t      <- F_f(t, vars, f_par)
  pars$MYZpar[[s]]$q_t      <- F_q(t, vars, q_par)
  pars$MYZpar[[s]]$g_t      <- F_g(t, vars, g_par)
  pars$MYZpar[[s]]$sigma_t  <- F_sigma(t, vars, sigma_par)
  pars$MYZpar[[s]]$mu       <- F_mu(t, vars, mu_par)
  pars$MYZpar[[s]]$nu       <- F_nu(t, vars, nu_par)
  pars$MYZpar[[s]]$eip      <- F_eip(t, vars, eip_par)
  pars$MYZpar[[s]]$calK     <- F_calK(t, vars, calK_par)
  pars$MYZpar[[s]]$eggsPerBatch <- eggsPerBatch
  # Reset Effect Sizes
  pars$MYZpar[[s]]$es_f     < rep(1, pars$nPatches)
  pars$MYZpar[[s]]$es_q     <- rep(1, pars$nPatches)
  pars$MYZpar[[s]]$es_g     <- rep(1, pars$nPatches)
  pars$MYZpar[[s]]$es_sigma <- rep(1, pars$nPatches)
  return(pars)
})}

#' @title Set mosquito bionomics to baseline
#' @description Implements [MBionomics] for models with no forcing on the baseline
#' @inheritParams ramp.xds::MBionomics
#' @return the model as a [list]
#' @export
MBionomics.SEI <- function(t, y, pars, s) {
  with(pars$MYZpar[[s]],{
    pars$MYZpar[[s]]$f <- es_f*f_t
    pars$MYZpar[[s]]$q <- es_q*q_t
    pars$MYZpar[[s]]$g <- es_g*g_t
    pars$MYZpar[[s]]$sigma <- es_sigma*sigma_t
    pars <- make_Omega(pars, s)
    return(pars)
})}


#' @title Get the feeding rate
#' @param pars an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_f.SEI = function(pars, s=1){
  with(pars$MYZpar[[s]], f_t*es_f)
}

#' @title Get the feeding rate
#' @param pars an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_q.SEI = function(pars, s=1){
  with(pars$MYZpar[[s]], q_t*es_q)
}

#' @title Get the feeding rate
#' @param pars an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_g.SEI = function(pars, s=1){
  with(pars$MYZpar[[s]], g_t*es_g)
}

#' @title Get the feeding rate
#' @param pars an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_sigma.SEI = function(pars, s=1){
  with(pars$MYZpar[[s]], sigma_t*es_sigma)
}
