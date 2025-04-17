# specialized methods for the adult mosquito SBRQ model


#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the SBRQ ODE model.
#' @inheritParams ramp.xds::dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.SBRQ <- function(t, y, pars, s){

  Lambda = as.vector(pars$Lambda[[s]])
  kappa = as.vector(pars$kappa[[s]])

  with(pars$ix$MYZ[[s]],{
    Su <- y[Su_ix]
    Bu <- y[Bu_ix]
    Ru <- y[Ru_ix]
    Qu <- y[Qu_ix]
    Sy <- y[Sy_ix]
    By <- y[By_ix]
    Ry <- y[Ry_ix]
    Qy <- y[Qy_ix]
    Sz <- y[Sz_ix]
    Bz <- y[Bz_ix]
    Rz <- y[Rz_ix]
    Qz <- y[Qz_ix]

    with(pars$MYZpar[[s]],{

      dSudt <- Lambda + zeta*Bu + nu*(1-p)*Qu-rho*Su-Omega_s %*% Su
      dBudt <- rho*Su + (nu*p + theta)*Qu + xi*R - (f+zeta)*Bu - Omega_b %*% Bu
      dRudt <- f*(1-q*kappa)*Bu - (xi + eta + g)*Ru
      dQudt <- eta*Ru - (nu + theta)*Qu - Omega_q %*% Qu

      dSydt <- zeta*By + nu*(1-p)*Qy-rho*Sy-Omega_s %*% Sy - tau*Sy
      dBydt <- rho*Sy + (nu*p + theta)*Qy + xi*R - (f+zeta)*By - Omega_b %*% By - tau*By
      dRydt <- f*q*kappa*Bu + f*By - (xi + eta + g)*Ry - tau*Ry
      dQydt <- eta*Ry - (nu + theta)*Qy - Omega_q %*% Qy - tau*Qy

      dSzdt <- zeta*Bz + nu*(1-p)*Qz-rho*Sz-Omega_s %*% Sz + tau*Sy
      dBzdt <- rho*Sz + (nu*p + theta)*Qz + xi*R - (f+zeta)*Bz - Omega_b %*% Bz + tau*By
      dRzdt <- f*Bz - (xi + eta + g)*Rz + tau*Ry
      dQzdt <- eta*Rz - (nu + theta)*Qz - Omega_q %*% Qz + tau*Qy

      return(c(dSudt, dBudt, dRudt, dQudt,
               dSydt, dBydt, dRydt, dQydt,
               dSzdt, dBzdt, dRzdt, dQzdt))
    })
  })
}

#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MBaseline] for the SBRQ model
#' @inheritParams ramp.xds::MBaseline
#' @return a named [list]
#' @export
MBaseline.SBRQ <- function(t, y, pars, s) {with(pars$MYZpar[[s]],{
  pars$MYZpar[[s]]$es_f       <- rep(1, nPatches)
  pars$MYZpar[[s]]$es_q       <- rep(1, nPatches)
  pars$MYZpar[[s]]$es_nu      <- rep(1, nPatches)
  pars$MYZpar[[s]]$es_g_s     <- rep(1, nPatches)
  pars$MYZpar[[s]]$es_g_b     <- rep(1, nPatches)
  pars$MYZpar[[s]]$es_g_r     <- rep(1, nPatches)
  pars$MYZpar[[s]]$es_g_q     <- rep(1, nPatches)
  pars$MYZpar[[s]]$es_sigma_b <- rep(1, nPatches)
  pars$MYZpar[[s]]$es_sigma_q <- rep(1, nPatches)
  pars$MYZpar[[s]]$es_sigma_s <- rep(1, nPatches)
  pars$MYZpar[[s]]$es_mu_b    <- rep(1, nPatches)
  pars$MYZpar[[s]]$es_mu_q    <- rep(1, nPatches)
  pars$MYZpar[[s]]$es_mu_s    <- rep(1, nPatches)
  pars$MYZpar[[s]]$es_theta   <- rep(1, nPatches)
  return(pars)
})}


#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MBionomics] for the SBRQ model
#' @inheritParams ramp.xds::MBionomics
#' @return a named [list]
#' @export
MBionomics.SBRQ <- function(t, y, pars, s) {with(pars$MYZpar[[s]],{
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
#' @description Implements [F_fqZ] for the SBRQ model.
#' @inheritParams ramp.xds::F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ.SBRQ <- function(t, y, pars, s) {
  f = get_f(pars, s)
  q = get_f(pars, s)
  fqZ = f*q*y[pars$ix$MYZ[[s]]$Z_b_ix]
  return(fqZ)
}

#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqM] for the SBRQ model.
#' @inheritParams ramp.xds::F_fqM
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqM.SBRQ <- function(t, y, pars, s) {
  f = get_f(pars, s)
  q = get_f(pars, s)
  M = with(pars$ix$MYZ[[s]], y[Bu_ix] + y[By_ix] + y[Bz_ix])
  fqM = f*q*M
  return(fqM)
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the SBRQ model.
#' @inheritParams ramp.xds::F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.SBRQ <- function(t, y, pars, s) {

  G <- with(pars$ix$MYZ[[s]], y[Qu_ix] + y[Qy_ix] + y[Qz_ix])
  eggs = with(pars$MYZpar[[s]], {
    return(G*nu*eggsPerBatch)
  })
  return(eggs)
}


#' @title Setup MYZpar for the SBRQ model
#' @description Implements [setup_MYZpar] for the RM model
#' @inheritParams ramp.xds::setup_MYZpar
#' @return a [list] vector
#' @export
setup_MYZpar.SBRQ = function(MYZname, pars, s, MYZopts=list()){
  pars$MYZpar[[s]] = make_MYZpar_SBRQ(pars$nPatches, MYZopts)
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
make_MYZpar_SBRQ = function(nPatches, MYZopts=list(), eip=12,
                          g=1/12, sigma_b=1/8, sigma_q=1/8, mu=0, f=0.5, q=0.95,
                          nu=1, eggsPerBatch=60){

  with(MYZopts,{
    MYZpar <- list()
    class(MYZpar) <- "SBRQ"

    MYZpar$nPatches <- nPatches

    eip_par <- list()
    class(eip_par) <- 'static'
    MYZpar$eip_par <- eip_par
    MYZpar$eip     <- eip

    MYZpar$f_t          <- checkIt(f, nPatches)
    MYZpar$es_f         <- rep(1, nPatches)

    MYZpar$nu_t         <- checkIt(nu, nPatches)
    MYZpar$es_nu        <- rep(1, nPatches)

    MYZpar$p            <- checkIt(p, nPatches)
    MYZpar$zeta         <- checkIt(zeta, nPatches)

    MYZpar$theta_t        <- checkIt(theta, nPatches)
    MYZpar$es_theta        <- rep(1, nPatches)

    MYZpar$g_b_t         <- checkIt(g_b, nPatches)
    MYZpar$es_g_b        <- rep(1, nPatches)

    MYZpar$g_q_t          <- checkIt(g_q, nPatches)
    MYZpar$es_g_q         <- rep(1, nPatches)

    MYZpar$g_s_t          <- checkIt(g_s, nPatches)
    MYZpar$es_g_s         <- rep(1, nPatches)

    MYZpar$g_r_t          <- checkIt(g_r, nPatches)
    MYZpar$es_g_r         <- rep(1, nPatches)

    MYZpar$sigma_b_t    <- checkIt(sigma_b, nPatches)
    MYZpar$es_sigma_b   <- rep(1, nPatches)

    MYZpar$sigma_q_t    <- checkIt(sigma_q, nPatches)
    MYZpar$es_sigma_q   <- rep(1, nPatches)

    MYZpar$sigma_s_t    <- checkIt(sigma_q, nPatches)
    MYZpar$es_sigma_s   <- rep(1, nPatches)

    MYZpar$mu_b_t    <- checkIt(mu_b, nPatches)
    MYZpar$es_mu_b   <- rep(1, nPatches)

    MYZpar$mu_q_t    <- checkIt(mu_q, nPatches)
    MYZpar$es_mu_q   <- rep(1, nPatches)

    MYZpar$mu_s_t    <- checkIt(mu_q, nPatches)
    MYZpar$es_mu_s   <- rep(1, nPatches)

    MYZpar$xi_t           <- checkIt(xi, nPatches)
    MYZpar$es_xi         <- rep(1, nPatches)

    MYZpar$eta_t           <- checkIt(eta, nPatches)
    MYZpar$es_eta         <- rep(1, nPatches)

    MYZpar$q_t          <- checkIt(q, nPatches)
    MYZpar$es_q         <- rep(1, nPatches)

    MYZpar$eggsPerBatch <- eggsPerBatch

    MYZpar$phi <- 1/MYZpar$eip

    calK <- diag(nPatches)

    MYZpar$calKb <- calK
    MYZpar$calKq <- calK
    MYZpar$calKs <- calK

    Omega_par <- list()
    class(Omega_par) <- "static"
    MYZpar$Omega_par <- Omega_par
    MYZpar$Omega_b <- with(MYZpar, compute_Omega_xde(g, sigma_b, mu, calK))
    MYZpar$Omega_q <- with(MYZpar, compute_Omega_xde(g, sigma_q, mu, calK))
    MYZpar$Omega_s <- with(MYZpar, compute_Omega_xde(g, sigma_s, mu, calK))
    base <- 'SBRQ'
    class(base) <- 'SBRQ'
    MYZpar$baseline <- base

    return(MYZpar)
})}

#' @title Return the parameters as a list
#' @description This method dispatches on the type of `pars$MYZpar[[s]]`.
#' @param pars an **`xds`** object
#' @param s the vector species index
#' @return a [list]
#' @export
get_MYZpars.SBRQ <- function(pars, s=1) {
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
set_MYZpars.SBRQ <- function(pars, s=1, MYZopts=list()) {
  nHabitats <- pars$nHabitats
  with(pars$MYZpar[[s]], with(MYZopts,{
    pars$MYZpar[[s]]$f_t = f
    pars$MYZpar[[s]]$q_t = q
    pars$MYZpar[[s]]$nu_t = nu
    pars$MYZpar[[s]]$g_s_t = g_s
    pars$MYZpar[[s]]$g_b_t = g_b
    pars$MYZpar[[s]]$g_r_t = g_r
    pars$MYZpar[[s]]$g_q_t = g_q
    pars$MYZpar[[s]]$mu_s_t = mu_s
    pars$MYZpar[[s]]$mu_b_t = mu_b
    pars$MYZpar[[s]]$mu_q_t = mu_q
    pars$MYZpar[[s]]$sigma_b_t = sigma_b
    pars$MYZpar[[s]]$sigma_q_t = sigma_q
    pars$MYZpar[[s]]$sigma_s_t = sigma_s
    pars$MYZpar[[s]]$eip = eip
    pars$MYZpar[[s]]$eggsPerBatch = eggsPerBatch
    return(pars)
  }))}

#' @title Setup initial values for the SBRQ model
#' @description Implements [setup_MYZinits] for the RM model
#' @inheritParams ramp.xds::setup_MYZinits
#' @return a [list]
#' @export
setup_MYZinits.SBRQ = function(pars, s, MYZopts=list()){
  pars$MYZinits[[s]] = with(pars$MYZpar[[s]], make_MYZinits_SBRQ(nPatches, MYZopts))
  return(pars)
}

#' @title Make inits for SBRQ adult mosquito model
#' @param nPatches the number of patches in the model
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param Su uninfected, sugar feeding mosquitoes
#' @param Bu uninfected, blood feeding mosquitoes
#' @param Ru uninfected, resting mosquitoes
#' @param Qu uninfected, gravid mosqutioes
#' @param Sy infected, sugar feeding mosquitoes
#' @param By infected, blood feeding mosquitoes
#' @param Ry infected, resting mosquitoes
#' @param Qy infected, gravid mosqutioes
#' @param Sz infectious, sugar feeding mosquitoes
#' @param Bz infectious, blood feeding mosquitoes
#' @param Rz infectious, resting mosquitoes
#' @param Qz infectious, gravid mosqutioes
#' @return a [list]
#' @export
make_MYZinits_SBRQ = function(nPatches, MYZopts = list(),
                              Su=100, Bu=100, Ru=1, Qu=10,
                              Sy=0, By=0, Ry=0, Qy=0,
                              Sz=0, Bz=0, Rz=0, Qz=0){

  with(MYZopts,{
    Su = checkit(Su, nPatches)
    Bu = checkit(Bu, nPatches)
    Ru = checkit(Ru, nPatches)
    Qu = checkit(Qu, nPatches)
    Sy = checkit(Sy, nPatches)
    By = checkit(By, nPatches)
    Ry = checkit(Ry, nPatches)
    Qy = checkit(Qy, nPatches)
    Sz = checkit(Sz, nPatches)
    Bz = checkit(Bz, nPatches)
    Rz = checkit(Rz, nPatches)
    Qz = checkit(Qz, nPatches)

    return(list(Su=Su, Bu=Bu, Ru=Ru, Qu=Qu,
                Sy=Sy, By=By, Ry=Ry, Qy=Qy,
                Sz=Sz, Bz=Bz, Rz=Rz, Qz=Qz))
  })}


#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [setup_MYZix] for the SBRQ model.
#' @inheritParams ramp.xds::setup_MYZix
#' @return none
#' @importFrom utils tail
#' @export
setup_MYZix.SBRQ <- function(pars, s) {with(pars,{


  Su_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Su_ix, 1)

  Bu_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Bu_ix, 1)

  Ru_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Ru_ix, 1)

  Qu_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Qu_ix, 1)

  Sy_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Sy_ix, 1)

  By_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(By_ix, 1)

  Ry_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Ry_ix, 1)

  Qy_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Qy_ix, 1)

  Sz_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Sz_ix, 1)

  Bz_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Bz_ix, 1)

  Rz_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Rz_ix, 1)

  Qz_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Qz_ix, 1)

  pars$max_ix = max_ix
  pars$ix$MYZ[[s]] = list(
    Su_ix = Su_ix, Bu_ix = Bu_ix, Ru_ix=Ru_ix, Qu_ix = Qu_ix,
    Sy_ix = Sy_ix, By_ix = By_ix, Ry_ix=Ry_ix, Qy_ix = Qy_ix,
    Sz_ix = Sz_ix, Bz_ix = Bz_ix, Rz_ix=Rz_ix, Qz_ix = Qz_ix)


  return(pars)
})}


#' @title Parse the output of deSolve and return variables for the SBRQ model
#' @description Implements [parse_MYZorbits] for the SBRQ model
#' @inheritParams ramp.xds::parse_MYZorbits
#' @return none
#' @export
parse_MYZorbits.SBRQ <- function(outputs, pars, s) {with(pars$ix$MYZ[[s]],{
  Su = outputs[,Su_ix]
  Bu = outputs[,Bu_ix]
  Ru = outputs[,Ru_ix]
  Qu = outputs[,Qu_ix]
  Sy = outputs[,Sy_ix]
  By = outputs[,By_ix]
  Ry = outputs[,Ry_ix]
  Qy = outputs[,Qy_ix]
  Sz = outputs[,Sz_ix]
  Bz = outputs[,Bz_ix]
  Rz = outputs[,Rz_ix]
  Qz = outputs[,Qz_ix]
  M =  Su+Bu+Ru+Qu+Sy+By+Ry+Qy+Sz+Bz+Rz+Qz
  Y =  Sy+By+Ry+Qy+Sz+Bz+Rz+Qz
  Z =  Sz+Bz+Rz+Qz
  y = Y/M
  z = Z/M
  return(list(
    Su=Su, Bu=Bu, Ru=Ru, Qu=Qu,
    Sy=Sy, By=By, Ry=Ry, Qy=Qy,
    Sz=Sz, Bz=Bz, Rz=Rz, Qz=Qz, M=M, Y=Y, Z=Z, y=y, z=z))
})}

#' @title Make inits for SBRQ adult mosquito model
#' @inheritParams ramp.xds::update_MYZinits
#' @return none
#' @export
update_MYZinits.SBRQ <- function(pars, y0, s) {
  with(pars$ix$MYZ[[s]],{
    Su = y[Su_ix]
    Bu = y[Bu_ix]
    Ru = y[Ru_ix]
    Qu = y[Qu_ix]
    Sy = y[Sy_ix]
    By = y[By_ix]
    Ry = y[Ry_ix]
    Qy = y[Qy_ix]
    Sz = y[Sz_ix]
    Bz = y[Bz_ix]
    Rz = y[Rz_ix]
    Qz = y[Qz_ix]
    pars = setup_MYZinits_SBRQ(pars$nPatches,
                               Su=Su, Bu=Bu, Ru=Ru, Qu=Qu,
                               Sy=Sy, By=By, Ry=Ry, Qy=Qy,
                               Sz=Sz, Bz=Bz, Rz=Rz, Qz=Qz)
    return(pars)
  })}


#' @title Set new MYZ parameter values
#' @description This method dispatches on the type of `pars$MYZpar[[s]]`.
#' @inheritParams ramp.xds::set_MYZinits
#' @return an `xds` object
#' @export
set_MYZinits.SBRQ <- function(pars, s=1, MYZopts=list()) {
  with(pars$MYZpar[[s]], with(MYZopts,{
    pars$MYZinits[[s]]$Su = Su
    pars$MYZinits[[s]]$Bu = Bu
    pars$MYZinits[[s]]$Ru = Ru
    pars$MYZinits[[s]]$Qu = Qu
    pars$MYZinits[[s]]$Sy = Sy
    pars$MYZinits[[s]]$By = By
    pars$MYZinits[[s]]$Ry = Ry
    pars$MYZinits[[s]]$Qy = Qy
    pars$MYZinits[[s]]$Sz = Sz
    pars$MYZinits[[s]]$Bz = Bz
    pars$MYZinits[[s]]$Rz = Rz
    pars$MYZinits[[s]]$Qz = Qz
    return(pars)
  }))}

#' @title Return initial values as a vector
#' @description Implements [get_MYZinits] for the SBRQ model.
#' @inheritParams ramp.xds::get_MYZinits
#' @return none
#' @export
get_MYZinits.SBRQ <- function(pars, s) {pars$MYZinits[[s]]}

#' @title Get the feeding rate
#' @param pars an **`xds`** object
#' @param s the vector species index
#' @return a [numeric] vector
#' @export
get_f.SBRQ = function(pars, s=1){
  with(pars$MYZpar[[s]], f_t*es_f)
}

#' @title Get the feeding rate
#' @param pars an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_q.SBRQ = function(pars, s=1){
  with(pars$MYZpar[[s]], q_t*es_q)
}

#' @title Get the feeding rate
#' @param pars an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_g.SBRQ = function(pars, s=1){
  with(pars$MYZpar[[s]], g_t*es_g)
}

#' @title Get the feeding rate
#' @param pars an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_sigma.SBRQ = function(pars, s=1){
  with(pars$MYZpar[[s]], sigma_b_t*es_sigma_b)
}
