# specialized methods for the adult mosquito macdonald model

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for Macdonald's DDE model.
#' @inheritParams ramp.xds::dMYZdt
#' @return a [numeric] vector
#' @importFrom deSolve lagvalue
#' @importFrom deSolve lagderiv
#' @export
dMYZdt.macdonald <- function(t, y, pars, s){
  Lambda = pars$Lambda[[s]]
  kappa = pars$kappa[[s]]

  with(list_MYZvars(y, pars, s),{
    with(pars$MYZpar[[s]],{

      if (t <= eip) {
        M_eip <- pars$MYZinits[[s]]$M
        Y_eip <- pars$MYZinits[[s]]$Y
        kappa_eip <- kappa
      } else {
        M_eip <- lagvalue(t=t-eip, nr = pars$ix$MYZ[[s]]$M_ix)
        Y_eip <- lagvalue(t=t-eip, nr = pars$ix$MYZ[[s]]$Y_ix)
        kappa_eip <- lagderiv(t=t-eip, nr = pars$ix$MYZ[[s]]$kappa_ix)
      }

      Omega_eip <- compute_Omega_xde(g_eip, sigma_eip, mu, calK)
      dMdt <- Lambda - (Omega %*% M)
      dPdt <- f*(M - P) - (Omega %*% P)
      dYdt <- f*q*kappa*(M - Y) - (Omega %*% Y)
      dZdt <- Upsilon %*% (f*q*kappa_eip * (M_eip - Y_eip)) - (Omega %*% Z)

      return(c(dMdt, dPdt, dYdt, dZdt, kappa))
    })
  })
}



#' @title Make inits for macdonald adult mosquito model
#' @param nPatches the number of patches in the model
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param M total mosquito density at each patch
#' @param P total parous mosquito density at each patch
#' @param Y infected mosquito density at each patch
#' @param Z infectious mosquito density at each patch
#' @return a [list]
#' @export
create_MYZinits_macdonald = function(nPatches, MYZopts = list(),
                              M=5, P=1, Y=1, Z=1){
  with(MYZopts,{
    M = unname(checkIt(M, nPatches))
    P = unname(checkIt(P, nPatches))
    Y = unname(checkIt(Y, nPatches))
    Z = unname(checkIt(Z, nPatches))
    dd = rep(0, nPatches)
    return(list(M=M, P=P, Y=Y, Z=Z, dd=dd))
  })
}

#' @title The net blood feeding rate of the infective mosquito population in a patch
#' @description Implements [F_fqZ] for the macdonald model.
#' @inheritParams ramp.xds::F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ.macdonald <- function(t, y, pars, s) {
  f = get_f(pars, s)
  q = get_q(pars, s)
  f*q*y[pars$ix$MYZ[[s]]$Z_ix]
}

#' @title The net blood feeding rate of the infective mosquito population in a patch
#' @description Implements [F_fqM] for the macdonald model.
#' @inheritParams ramp.xds::F_fqM
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqM.macdonald <- function(t, y, pars, s) {
  f = get_f(pars, s)
  q = get_q(pars, s)
  f*q*y[pars$ix$MYZ[[s]]$M_ix]
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the macdonald model.
#' @inheritParams ramp.xds::F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.macdonald <- function(t, y, pars, s) {
  M <- y[pars$ix$MYZ[[s]]$M_ix]
  with(pars$MYZpar[[s]],{
    return(M*nu*eggsPerBatch)
  })
}


#' @title Setup MYZpar for the macdonald model
#' @description Implements [make_MYZpar] for the macdonald model
#' @inheritParams ramp.xds::make_MYZpar
#' @return a [list] vector
#' @export
make_MYZpar.macdonald = function(MYZname, pars, s, MYZopts=list()){
  pars$dlay <- 'dde'
  class(pars$dlay) <- 'dde'
  MYZpar <- create_MYZpar_macdonald(pars$nPatches, MYZopts)
  class(MYZpar) <- 'macdonald'
  pars$MYZpar[[s]] = MYZpar
  return(pars)
}

#' @title Make parameters for macdonald ODE adult mosquito model
#' @param nPatches is the number of patches, an integer
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param eip the extrinsic incubation period
#' @param g mosquito mortality rate
#' @param sigma emigration rate
#' @param mu fraction lost during emigration
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @return a [list]
#' @importFrom expm expm
#' @export
create_MYZpar_macdonald = function(nPatches, MYZopts=list(), eip=12,
                            g=1/12, sigma=1/8, mu=0, f=0.3, q=0.95,
                            nu=1, eggsPerBatch=60){

  with(MYZopts,{
    MYZpar <- list()
    MYZpar$nPatches <- nPatches

    MYZpar$eip    <- checkIt(eip, 1)
    MYZpar$f      <- checkIt(f, nPatches)
    MYZpar$q      <- checkIt(q, nPatches)
    MYZpar$g      <- checkIt(g, nPatches)
    MYZpar$sigma  <- checkIt(sigma, nPatches)
    MYZpar$mu     <- checkIt(mu, nPatches)
    MYZpar$nu     <- checkIt(nu, nPatches)

    MYZpar$eggsPerBatch <- eggsPerBatch
    MYZpar$calK <- diag(nPatches)

    Omega <- diag(g, nPatches)
    MYZpar$Omega <- Omega
    MYZpar$Upsilon <- expm::expm(-Omega*eip)

    MYZpar$baseline <- 'macdonald'
    class(MYZpar$baseline) <- 'macdonald'

    return(MYZpar)
})}

#' @title Setup initial values for the macdonald model
#' @description Implements [make_MYZinits] for the macdonald model
#' @inheritParams ramp.xds::make_MYZinits
#' @return a [list]
#' @export
make_MYZinits.macdonald = function(pars, s, MYZopts=list()){with(pars$MYZpar[[s]], {
  pars$MYZinits[[s]] = create_MYZinits_macdonald(nPatches, MYZopts)
  return(pars)
})}

#' @title Return initial values as a list
#' @description Implements [get_MYZinits] for the macdonald model.
#' @inheritParams ramp.xds::get_MYZinits
#' @return [list]
#' @export
get_MYZinits.macdonald <- function(pars, s=1){
  pars$MYZinits[[s]]
}

#' @title Update inits for macdonald adult mosquito model
#' @inheritParams ramp.xds::update_MYZinits
#' @return a [list]
#' @export
update_MYZinits.macdonald <- function(pars, y0, s) {
  with(pars$ix$MYZ[[s]],{
    M = y0[M_ix]
    P = y0[P_ix]
    Y = y0[Y_ix]
    Z = y0[Z_ix]
    pars = make_MYZinits(pars, s, list(M=M, P=P, Y=Y, Z=Z))
    return(pars)
  })}

#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$MYZpar[[s]]`
#' @inheritParams ramp.xds::list_MYZvars
#' @return a [list]
#' @export
list_MYZvars.macdonald <- function(y, pars, s){
  with(pars$ix$MYZ[[s]],
       return(list(
         M = y[M_ix],
         P = y[P_ix],
         Y = y[Y_ix],
         Z = y[Z_ix]
       )))
}


#' @title Parse the output of deSolve and return variables for the macdonald model
#' @description Implements [parse_MYZorbits] for the macdonald model
#' @inheritParams ramp.xds::parse_MYZorbits
#' @return a [list]
#' @export
parse_MYZorbits.macdonald <- function(outputs, pars, s) {with(pars$ix$MYZ[[s]],{
  M = outputs[,M_ix]
  P = outputs[,P_ix]
  Y = outputs[,Y_ix]
  Z = outputs[,Z_ix]
  f = get_ft(pars, s)
  q = get_qt(pars, s)
  y = Y/M
  z = Z/M
  parous = P/M
  return(list(M=M, P=P, Y=Y, Z=Z, y=y, z=z, parous=parous, fqZ=f*q*Z, fqM=f*q*M))
})}

#' @title Compute the steady states as a function of the daily EIR
#' @description This method dispatches on the type of `MYZpar`.
#' @inheritParams ramp.xds::xde_steady_state_MYZ
#' @return none
#' @importFrom MASS ginv
#' @importFrom expm expm
#' @export
xde_steady_state_MYZ.macdonald = function(Lambda, kappa, MYZpar){with(MYZpar,{
  kappa = as.vector(kappa); Lambda = as.vector(Lambda)
  Omega_inv <- ginv(Omega)
  Upsilon = expm(-Omega*eip)
  M_eq <- as.vector(Omega_inv %*% Lambda)
  P_eq <- as.vector(ginv(diag(f, nPatches) + Omega) %*% diag(f, nPatches) %*% M_eq)
  Y_eq <- as.vector(ginv(diag(f*q*kappa) + Omega) %*% diag(f*q*kappa) %*% M_eq)
  Z_eq <- as.vector(Omega_inv %*% Upsilon %*% diag(f*q*kappa) %*% (M_eq - Y_eq))
  return(list(M=M_eq, P=P_eq, Y=Y_eq, Z=Z_eq))
})}


#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for the macdonald model.
#' @inheritParams ramp.xds::make_indices_MYZ
#' @return a [list]
#' @importFrom utils tail
#' @export
make_indices_MYZ.macdonald <- function(pars, s) {with(pars,{

  M_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(M_ix, 1)

  P_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(P_ix, 1)

  Y_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Y_ix, 1)

  Z_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Z_ix, 1)

  fqkappa_ix <- seq(from = max_ix+1, length.out = nPatches)
  max_ix <- tail(fqkappa_ix, 1)

  g_ix <- seq(from = max_ix+1, length.out = nPatches)
  max_ix <- tail(g_ix, 1)

  sigma_ix <- seq(from = max_ix+1, length.out = nPatches)
  max_ix <- tail(sigma_ix, 1)

  pars$max_ix = max_ix
  pars$ix$MYZ[[s]] = list(M_ix=M_ix, P_ix=P_ix, Y_ix=Y_ix, Z_ix=Z_ix,
                          fqkappa_ix=fqkappa_ix,
                          g_ix=g_ix, sigma_ix=sigma_ix)
  return(pars)
})}


#' @title Macdonald-style adult mosquito bionomics
#' @description Reset the effect sizes for static models
#' @description
#' When modules are added to compute effect sizes
#' from baseline parameters, those functions store
#' an effect size. The total effect size is the
#' product of the effect sizes for each intervention.
#' Since coverage could be changing dynamically, these
#' must be reset each time the derivatives are computed.
#' @inheritParams ramp.xds::MBionomics
#' @return an **`xds`** object
#' @export
MBionomics.macdonald <- function(t, y, pars, s) {
  return(pars)
}


#' @title Get the feeding rate
#' @param pars an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_f.macdonald = function(pars, s=1){
  with(pars$MYZpar[[s]], f)
}

#' @title Get the feeding rate
#' @param pars an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_q.macdonald = function(pars, s=1){
  with(pars$MYZpar[[s]], q)
}

#' @title Get the feeding rate
#' @param pars an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_g.macdonald = function(pars, s=1){
  with(pars$MYZpar[[s]], g)
}

#' @title Get the feeding rate
#' @param pars an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_sigma.macdonald = function(pars, s=1){
  with(pars$MYZpar[[s]], sigma)
}
