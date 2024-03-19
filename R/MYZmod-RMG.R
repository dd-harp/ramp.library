# specialized methods for the adult mosquito RMG model

#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MBionomics] for the RMG model
#' @inheritParams ramp.xde::MBionomics
#' @return a named [list]
#' @export
MBionomics.RMG <- function(t, y, pars, s) {
  with(pars$MYZpar[[s]],{
    pars$MYZpar[[s]]$f <- f0
    pars$MYZpar[[s]]$q <- q0
    pars$MYZpar[[s]]$g <- g0
    pars$MYZpar[[s]]$phi <- phi0
    pars$MYZpar[[s]]$sigma <- sigma0
    pars$MYZpar[[s]]$nu <- nu0

  return(pars)
})}

#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqZ] for the RMG model.
#' @inheritParams ramp.xde::F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ.RMG <- function(t, y, pars, s) {
  fqZ = with(pars$MYZpar[[s]], f*q)*y[pars$ix$MYZ[[s]]$Z_b_ix]
  return(fqZ)
}

#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqM] for the RMG model.
#' @inheritParams ramp.xde::F_fqM
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqM.RMG <- function(t, y, pars, s) {
  M = with(pars$ix$MYZ[[s]], y[U_b_ix] + y[Y_b_ix] + y[Z_b_ix])
  fqM = with(pars$MYZpar[[s]], f*q)*M
  return(fqM)
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the RMG model.
#' @inheritParams ramp.xde::F_eggs
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
#' @inheritParams ramp.xde::dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.RMG <- function(t, y, pars, s){

  Lambda = pars$Lambda[[s]]
  kappa = pars$kappa[[s]]

  with(pars$ix$MYZ[[s]],{
    U_b <- y[U_b_ix]
    U_g <- y[U_g_ix]
    Y_b <- y[Y_b_ix]
    Y_g <- y[Y_g_ix]
    Z_b <- y[Z_b_ix]
    Z_g <- y[Z_g_ix]

    with(pars$MYZpar[[s]],{

      Omega_b <- make_Omega(g, sigma, calKb, nPatches)
      Omega_q <- make_Omega(g, sigma, calKq, nPatches)

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
#' @inheritParams ramp.xde::setup_MYZpar
#' @return a [list] vector
#' @export
setup_MYZpar.RMG = function(MYZname, pars, s, MYZopts=list(), EIPmod, calK){
  pars$MYZpar[[s]] = make_MYZpar_RMG(pars$nPatches, MYZopts, EIPmod, calK)
  return(pars)
}

#' @title Make parameters for RM ODE adult mosquito model
#' @param nPatches is the number of patches, an integer
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param EIPmod a [list] that defines the EIP model
#' @param calK a mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param g mosquito mortality rate
#' @param sigma emigration rate
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @return a [list]
#' @export
make_MYZpar_RMG = function(nPatches, MYZopts=list(), EIPmod, calK,
                          g=1/12, sigma=1/8, f=0.5, q=0.95,
                          nu=1, eggsPerBatch=60){

  stopifnot(is.matrix(calK))
  stopifnot(dim(calK) == c(nPatches, nPatches))

  with(MYZopts,{
    MYZpar <- list()
    MYZpar$xde <- "ode"
    class(MYZpar$xde) <- "ode"
    class(MYZpar) <- "RMG"

    MYZpar$nPatches <- nPatches

    MYZpar$g      <- checkIt(g, nPatches)
    MYZpar$sigma  <- checkIt(sigma, nPatches)
    MYZpar$f      <- checkIt(f, nPatches)
    MYZpar$q      <- checkIt(q, nPatches)
    MYZpar$nu     <- checkIt(nu, nPatches)
    MYZpar$eggsPerBatch <- eggsPerBatch
    MYZpar$phi <- 1/MYZpar$eip

    # Store as baseline values
    MYZpar$g0      <- MYZpar$g
    MYZpar$sigma0  <- MYZpar$sigma
    MYZpar$f0      <- MYZpar$f
    MYZpar$q0      <- MYZpar$q
    MYZpar$nu0     <- MYZpar$nu

    # The EIP model and the eip
    MYZpar$EIPmod <- EIPmod
    MYZpar$eip <- EIP(0, EIPmod)
    MYZpar$phi0 <- 1/MYZpar$eip

    MYZpar$calKb <- calK
    MYZpar$calKq <- calK

    MYZpar$Omega_b <- make_Omega(g, sigma, calK, nPatches)
    MYZpar$Omega_q <- make_Omega(g, sigma, calK, nPatches)

    return(MYZpar)
})}

#' @title Setup initial values for the RMG model
#' @description Implements [setup_MYZinits] for the RM model
#' @inheritParams ramp.xde::setup_MYZinits
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
#' @description Implements [make_indices_MYZ] for the RMG model.
#' @inheritParams ramp.xde::make_indices_MYZ
#' @return none
#' @importFrom utils tail
#' @export
make_indices_MYZ.RMG <- function(pars, s) {with(pars,{

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

#' @title Make parameters for RMG ODE adult mosquito model
#' @param pars a [list]
#' @param g mosquito mortality rate
#' @param sigma emigration rate
#' @param calK mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @param eip length of extrinsic incubation period
#' @return none
#' @export
make_parameters_MYZ_RMG <- function(pars, g, sigma, f, q, nu, eggsPerBatch, eip, calK) {
  stopifnot(is.numeric(g), is.numeric(sigma), is.numeric(f), is.numeric(q), is.numeric(nu), is.numeric(eggsPerBatch))

   nPatches = pars$nPatches

  MYZpar <- list()
  MYZpar$xde = 'ode'
  class(MYZpar$xde) <- 'ode'
  class(MYZpar) <- "RMG"

  MYZpar$g      <- checkIt(g, nPatches)
  MYZpar$sigma  <- checkIt(sigma, nPatches)
  MYZpar$f      <- checkIt(f, nPatches)
  MYZpar$q      <- checkIt(q, nPatches)
  MYZpar$nu     <- checkIt(nu, nPatches)
  MYZpar$eggsPerBatch <- eggsPerBatch

  # Store as baseline values
  MYZpar$g0      <- MYZpar$g
  MYZpar$sigma0  <- MYZpar$sigma
  MYZpar$f0      <- MYZpar$f
  MYZpar$q0      <- MYZpar$q
  MYZpar$nu0     <- MYZpar$nu

  pars$MYZpar[[1]] <- MYZpar

  return(pars)
}

#' @title Parse the output of deSolve and return variables for the RMG model
#' @description Implements [parse_deout_MYZ] for the RMG model
#' @inheritParams ramp.xde::parse_deout_MYZ
#' @return none
#' @export
parse_deout_MYZ.RMG <- function(deout, pars, s) {
  time = deout[,1]
  with(pars$ix$MYZ[[s]],{
    U_b = deout[,U_b_ix+1]
    U_g = deout[,U_g_ix+1]
    Y_b = deout[,Y_b_ix+1]
    Y_g = deout[,Y_g_ix+1]
    Z_b = deout[,Z_b_ix+1]
    Z_g = deout[,Z_g_ix+1]
    M =  U_b+Y_b+Z_b+U_g+Y_g+Z_g
    Y =  Y_b+Y_g
    Z =  Z_b+Z_g
    y = Y/M
    z = Z/M
  return(list(time=time, U_b=U_b, U_g=U_g, Y_b=Y_b, Y_g=Y_g, Z_b=Z_b, Z_g=Z_g, M=M, Y=Y, Z=Z, y=y, z=z))
})}

#' @title Make inits for RMG adult mosquito model
#' @inheritParams ramp.xde::update_inits_MYZ
#' @return none
#' @export
update_inits_MYZ.RMG <- function(pars, y0, s) {
  with(pars$ix$MYZ[[s]],{
    U_b = y[U_b_ix]
    U_g = y[U_g_ix]
    Y_b = y[Y_b_ix]
    Y_g = y[Y_g_ix]
    Z_b = y[Z_b_ix]
    Z_g = y[Z_g_ix]
    pars = make_MYZinits_RMG(pars$nPatches, U_b0=U_b, U_g0=U_g, Y_b0=Y_b, Y_g0=Y_g, Z_b0=Z_b, Z_g0=Z_g)
    return(pars)
})}

#' @title Make inits for RMG adult mosquito model
#' @param pars a [list]
#' @param U_b0 total mosquito density at each patch
#' @param U_g0 total gravid uninfected mosquito density at each patch
#' @param Y_b0 infected mosquito density at each patch
#' @param Y_g0 total gravid infected mosquito density at each patch
#' @param Z_b0 infectious mosquito density at each patch
#' @param Z_g0 total gravid infectious mosquito density at each patch
#' @return none
#' @export
make_inits_MYZ_RMG <- function(pars, U_b0, U_g0, Y_b0, Y_g0, Z_b0, Z_g0) {
  pars$MYZinits[[1]] = list(U_b=U_b0, U_g=U_g0, Y_b=Y_b0, Y_g=Y_g0, Z_b=Z_b0, Z_g=Z_g0)
  return(pars)
}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the RMG model.
#' @inheritParams ramp.xde::get_inits_MYZ
#' @return none
#' @export
get_inits_MYZ.RMG <- function(pars, s) {with(pars$MYZinits[[s]],{
  c(U_b, U_g, Y_b, Y_g, Z_b, Z_g)
})}

