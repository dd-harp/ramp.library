# specialized methods for the adult mosquito RMG model

#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MBaseline] for the RMG model
#' @inheritParams ramp.xds::MBaseline
#' @return a named [list]
#' @export
MBaseline.RMG <- function(t, y, pars, s) {with(pars$MYZpar[[s]],{
  pars$MYZpar[[s]]$es_g       <- rep(1, nPatches)
  pars$MYZpar[[s]]$es_sigma_b <- rep(1, nPatches)
  pars$MYZpar[[s]]$es_sigma_q <- rep(1, nPatches)
  pars$MYZpar[[s]]$es_f       <- rep(1, nPatches)
  pars$MYZpar[[s]]$es_q       <- rep(1, nPatches)
  return(pars)
})}


#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MBionomics] for the RMG model
#' @inheritParams ramp.xds::MBionomics
#' @return a named [list]
#' @export
MBionomics.RMG <- function(t, y, pars, s) {with(pars$MYZpar[[s]],{
  pars$MYZpar[[s]]$f <- es_f*f_t
  pars$MYZpar[[s]]$q <- es_q*q_t
  g <- es_g*g_t
  sigma_b <- es_sigma_b*sigma_b_t
  sigma_q <- es_sigma_q*sigma_q_t
  pars$MYZpar[[s]]$g <- g
  pars$MYZpar[[s]]$sigma_b <- sigma_b
  pars$MYZpar[[s]]$sigma_q <- sigma_q
  pars$MYZpar[[s]]$Omega_b = compute_Omega_xde(g, sigma_b, mu, calKb)
  pars$MYZpar[[s]]$Omega_q = compute_Omega_xde(g, sigma_q, mu, calKq)
  return(pars)
})}

#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqZ] for the RMG model.
#' @inheritParams ramp.xds::F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ.RMG <- function(t, y, pars, s) {
  f = get_f(pars, s)
  q = get_f(pars, s)
  fqZ = f*q*y[pars$ix$MYZ[[s]]$Z_b_ix]
  return(fqZ)
}

#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqM] for the RMG model.
#' @inheritParams ramp.xds::F_fqM
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqM.RMG <- function(t, y, pars, s) {
  f = get_f(pars, s)
  q = get_f(pars, s)
  M = with(pars$ix$MYZ[[s]], y[U_b_ix] + y[Y_b_ix] + y[Z_b_ix])
  fqM = f*q*M
  return(fqM)
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the RMG model.
#' @inheritParams ramp.xds::F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.RMG <- function(t, y, pars, s) {

  G <- with(pars$ix$MYZ[[s]], y[U_g_ix] + y[Y_g_ix] + y[Z_g_ix])
  eggs = with(pars$MYZpar[[s]], {
    return(G*nu*eggsPerBatch)
  })
  return(eggs)
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the RMG ODE model.
#' @inheritParams ramp.xds::dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.RMG <- function(t, y, pars, s){

  Lambda = as.vector(pars$Lambda[[s]])
  kappa = as.vector(pars$kappa[[s]])

  with(pars$ix$MYZ[[s]],{
    U_b <- y[U_b_ix]
    U_g <- y[U_g_ix]
    Y_b <- y[Y_b_ix]
    Y_g <- y[Y_g_ix]
    Z_b <- y[Z_b_ix]
    Z_g <- y[Z_g_ix]

    with(pars$MYZpar[[s]],{

      dU_bdt <- Lambda + nu*U_g - f*U_b - (Omega_b %*% U_b)
      dU_gdt <- f*(1-q*kappa)*U_b - nu*U_g - (Omega_q %*% U_g)
      dY_bdt <- nu*Y_g - f*Y_b - phi*Y_b - (Omega_b %*% Y_b)
      dY_gdt <- f*q*kappa*U_b + f*Y_b - nu*Y_g - phi*Y_g - (Omega_q %*% Y_g)
      dZ_bdt <- phi*Y_b + nu*Z_g - f*Z_b - (Omega_b %*% Z_b)
      dZ_gdt <- phi*Y_g + f*Z_b - nu*Z_g - (Omega_q %*% Z_g)

      return(c(dU_bdt, dU_gdt, dY_bdt, dY_gdt, dZ_bdt, dZ_gdt))
    })
  })
}

#' @title Setup MYZpar for the RMG model
#' @description Implements [setup_MYZpar] for the RM model
#' @inheritParams ramp.xds::setup_MYZpar
#' @return a [list] vector
#' @export
setup_MYZpar.RMG = function(MYZname, pars, s, MYZopts=list()){
  pars$MYZpar[[s]] = make_MYZpar_RMG(pars$nPatches, MYZopts)
  return(pars)
}

#' @title Make parameters for RM ODE adult mosquito model
#' @param nPatches is the number of patches, an integer
#' @param MYZopts a [list] of values that overwrites the defaults
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
#' @export
make_MYZpar_RMG = function(nPatches, MYZopts=list(), eip=12,
                          g=1/12, sigma_b=1/8, sigma_q=1/8, mu=0, f=0.5, q=0.95,
                          nu=1, eggsPerBatch=60){

  with(MYZopts,{
    MYZpar <- list()
    class(MYZpar) <- "RMG"

    MYZpar$nPatches <- nPatches

    eip_par <- list()
    class(eip_par) <- 'static'
    MYZpar$eip_par <- eip_par
    MYZpar$eip     <- eip

    MYZpar$g_t          <- checkIt(g, nPatches)
    MYZpar$es_g         <- rep(1, nPatches)
    MYZpar$sigma_b_t    <- checkIt(sigma_b, nPatches)
    MYZpar$es_sigma_b   <- rep(1, nPatches)
    MYZpar$sigma_q_t    <- checkIt(sigma_q, nPatches)
    MYZpar$es_sigma_q   <- rep(1, nPatches)
    MYZpar$mu           <- checkIt(mu, nPatches)
    MYZpar$f_t          <- checkIt(f, nPatches)
    MYZpar$es_f         <- rep(1, nPatches)
    MYZpar$q_t          <- checkIt(q, nPatches)
    MYZpar$es_q         <- rep(1, nPatches)
    MYZpar$nu           <- checkIt(nu, nPatches)
    MYZpar$eggsPerBatch <- eggsPerBatch
    MYZpar$phi <- 1/MYZpar$eip

    calK <- diag(nPatches)

    MYZpar$calKb <- calK
    MYZpar$calKq <- calK

    Omega_par <- list()
    class(Omega_par) <- "static"
    MYZpar$Omega_par <- Omega_par
    MYZpar$Omega_b <- with(MYZpar, compute_Omega_xde(g, sigma_b, mu, calK))
    MYZpar$Omega_q <- with(MYZpar, compute_Omega_xde(g, sigma_q, mu, calK))
    base <- 'RMG'
    class(base) <- 'RMG'
    MYZpar$baseline <- base

    return(MYZpar)
})}

#' @title Return the parameters as a list
#' @description This method dispatches on the type of `pars$MYZpar[[s]]`.
#' @param pars an **`xds`** object
#' @param s the vector species index
#' @return a [list]
#' @export
get_MYZpars.RMG <- function(pars, s=1) {
  with(pars$MYZpar[[s]], list(
    f=f_t, q=q_t, g=g_t, sigma_b=sigma_b_t, sigma_q = sigma_q_t, eip=eip, mu=mu_t,
    nu=nu_t, eggsPerBatch=eggsPerBatch, calK=calK
  ))
}

#' @title Return the parameters as a list
#' @description This method dispatches on the type of `pars$MYZpar[[s]]`.
#' @inheritParams ramp.xds::set_MYZpars
#' @return an **`xds`** object
#' @export
set_MYZpars.RMG <- function(pars, s=1, MYZopts=list()) {
  nHabitats <- pars$nHabitats
  with(pars$MYZpar[[s]], with(MYZopts,{
    pars$MYZpar[[s]]$f_t = f
    pars$MYZpar[[s]]$q_t = q
    pars$MYZpar[[s]]$g_t = g
    pars$MYZpar[[s]]$sigma_b_t = sigma_b
    pars$MYZpar[[s]]$sigma_q_t = sigma_q
    pars$MYZpar[[s]]$eip = eip
    pars$MYZpar[[s]]$mu = mu
    pars$MYZpar[[s]]$nu_t = nu
    pars$MYZpar[[s]]$eggsPerBatch = eggsPerBatch
    return(pars)
  }))}

#' @title Setup initial values for the RMG model
#' @description Implements [setup_MYZinits] for the RM model
#' @inheritParams ramp.xds::setup_MYZinits
#' @return a [list]
#' @export
setup_MYZinits.RMG = function(pars, s, MYZopts=list()){
  pars$MYZinits[[s]] = with(pars$MYZpar[[s]], make_MYZinits_RMG(nPatches, MYZopts))
  return(pars)
}

#' @title Make inits for RMG adult mosquito model
#' @param nPatches the number of patches in the model
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param U_b0 total uninfected, not gravid mosquito density at each patch
#' @param U_g0 total uninfected, gravid mosquito density at each patch
#' @param Y_b0 infected, not gravid mosquito density at each patch
#' @param Y_g0 infected, gravid mosquito density at each patch
#' @param Z_b0 infectious, not gravid mosquito density at each patch
#' @param Z_g0 infectious, gravid mosquito density at each patch
#' @return a [list]
#' @export
make_MYZinits_RMG = function(nPatches, MYZopts = list(),
                             U_b0=5, U_g0=1, Y_b0=1, Y_g0=1, Z_b0=1, Z_g0=1){
  with(MYZopts,{
    U_b = checkIt(U_b0, nPatches)
    U_g = checkIt(U_g0, nPatches)
    Y_b = checkIt(Y_b0, nPatches)
    Y_g = checkIt(Y_g0, nPatches)
    Z_b = checkIt(Z_b0, nPatches)
    Z_g = checkIt(Z_g0, nPatches)

    return(list(U_b=U_b, U_g=U_g, Y_b=Y_b, Y_g=Y_g, Z_b=Z_b, Z_g=Z_g))
  })}

#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [setup_MYZix] for the RMG model.
#' @inheritParams ramp.xds::setup_MYZix
#' @return none
#' @importFrom utils tail
#' @export
setup_MYZix.RMG <- function(pars, s) {with(pars,{

  U_b_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(U_b_ix, 1)

  U_g_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(U_g_ix, 1)

  Y_b_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Y_b_ix, 1)

  Y_g_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Y_g_ix, 1)

  Z_b_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Z_b_ix, 1)

  Z_g_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Z_g_ix, 1)

  pars$max_ix = max_ix
  pars$ix$MYZ[[s]] = list(U_b_ix=U_b_ix, U_g_ix=U_g_ix,
                          Y_b_ix=Y_b_ix, Y_g_ix=Y_g_ix,
                          Z_b_ix=Z_b_ix, Z_g_ix=Z_g_ix)
  return(pars)
})}


#' @title Parse the output of deSolve and return variables for the RMG model
#' @description Implements [parse_MYZorbits] for the RMG model
#' @inheritParams ramp.xds::parse_MYZorbits
#' @return none
#' @export
parse_MYZorbits.RMG <- function(outputs, pars, s) {with(pars$ix$MYZ[[s]],{
    U_b = outputs[,U_b_ix]
    U_g = outputs[,U_g_ix]
    Y_b = outputs[,Y_b_ix]
    Y_g = outputs[,Y_g_ix]
    Z_b = outputs[,Z_b_ix]
    Z_g = outputs[,Z_g_ix]
    M =  U_b+Y_b+Z_b+U_g+Y_g+Z_g
    Y =  Y_b+Y_g
    Z =  Z_b+Z_g
    y = Y/M
    z = Z/M
  return(list(U_b=U_b, U_g=U_g, Y_b=Y_b, Y_g=Y_g, Z_b=Z_b, Z_g=Z_g, M=M, Y=Y, Z=Z, y=y, z=z))
})}

#' @title Make inits for RMG adult mosquito model
#' @inheritParams ramp.xds::update_MYZinits
#' @return none
#' @export
update_MYZinits.RMG <- function(pars, y0, s) {
  with(pars$ix$MYZ[[s]],{
    U_b = y[U_b_ix]
    U_g = y[U_g_ix]
    Y_b = y[Y_b_ix]
    Y_g = y[Y_g_ix]
    Z_b = y[Z_b_ix]
    Z_g = y[Z_g_ix]
    pars = setup_MYZinits_RMG(pars$nPatches, U_b0=U_b, U_g0=U_g, Y_b0=Y_b, Y_g0=Y_g, Z_b0=Z_b, Z_g0=Z_g)
    return(pars)
})}

#' @title Set new MYZ parameter values
#' @description This method dispatches on the type of `pars$MYZpar[[s]]`.
#' @inheritParams ramp.xds::set_MYZinits
#' @return an `xds` object
#' @export
set_MYZinits.RMG <- function(pars, s=1, MYZopts=list()) {
  with(pars$MYZpar[[s]], with(MYZopts,{
    pars$MYZinits[[s]]$U_b = U_b
    pars$MYZinits[[s]]$U_g = U_g
    pars$MYZinits[[s]]$Y_b = Y_b
    pars$MYZinits[[s]]$Y_g = Y_g
    pars$MZZinits[[s]]$Z_b = Z_b
    pars$MYZinits[[s]]$Y_g = Y_g
    return(pars)
  }))}

#' @title Return initial values as a vector
#' @description Implements [get_MYZinits] for the RMG model.
#' @inheritParams ramp.xds::get_MYZinits
#' @return none
#' @export
get_MYZinits.RMG <- function(pars, s) {pars$MYZinits[[s]]}

#' @title Get the feeding rate
#' @param pars an **`xds`** object
#' @param s the vector species index
#' @return a [numeric] vector
#' @export
get_f.RMG = function(pars, s=1){
  with(pars$MYZpar[[s]], f_t*es_f)
}

#' @title Get the feeding rate
#' @param pars an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_q.RMG = function(pars, s=1){
  with(pars$MYZpar[[s]], q_t*es_q)
}

#' @title Get the feeding rate
#' @param pars an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_g.RMG = function(pars, s=1){
  with(pars$MYZpar[[s]], g_t*es_g)
}

#' @title Get the feeding rate
#' @param pars an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_sigma.RMG = function(pars, s=1){
  with(pars$MYZpar[[s]], sigma_b_t*es_sigma_b)
}
