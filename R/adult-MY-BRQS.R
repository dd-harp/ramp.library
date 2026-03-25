# specialized methods for the adult mosquito BRQS model

#' @title The `BRQS` module for the MY component
#' @description
#' Implements the **MY** component using a BRQS (Blood-feeding / Resting /
#' gravid Queue / Sugar-feeding) model of adult mosquito ecology and infection
#' dynamics. Mosquitoes cycle through four behavioural states, and infection
#' is tracked across all states with rate \eqn{\tau} governing progression from
#' exposed (**y**) to infectious (**z**).
#'
#' @section State Variables:
#' \describe{
#'   \item{`Su`, `Bu`, `Ru`, `Qu`}{uninfected mosquitoes in the sugar-feeding, blood-feeding, resting, and gravid states}
#'   \item{`Sy`, `By`, `Ry`, `Qy`}{infected (exposed) mosquitoes in each behavioural state}
#'   \item{`Sz`, `Bz`, `Rz`, `Qz`}{infectious mosquitoes in each behavioural state}
#' }
#'
#' @section Parameters:
#' \describe{
#'   \item{`f`}{blood feeding rate}
#'   \item{`q`}{human blood fraction}
#'   \item{`g`}{background mortality rate}
#'   \item{`rho`}{rate of leaving the sugar-feeding state to blood-feed}
#'   \item{`zeta`}{rate at which blood-feeding mosquitoes take a sugar meal}
#'   \item{`xi`}{rate of leaving the resting state}
#'   \item{`eta`}{rate of transition from resting to gravid}
#'   \item{`nu`, `theta`}{parameters governing the egg-laying / gravid state}
#'   \item{`delta`}{fraction of gravid mosquitoes returning to blood feeding}
#'   \item{`omega`}{fraction of emerging mosquitoes entering the sugar-feeding state}
#'   \item{`tau`}{rate of progression from exposed to infectious (\eqn{1/\tau} = EIP)}
#'   \item{`sigma_b`, `sigma_q`, `sigma_s`}{emigration rates in each state}
#' }
#'
#' @section Dynamics:
#' \deqn{
#' \begin{array}{rl}
#' dB_u/dt &= (1-\omega)\Lambda + \rho S_u + (\nu+\theta)\delta Q_u + \xi R - (f+\zeta) B_u - \Omega_b \cdot B_u \\
#' dR_u/dt &= f(1-q\kappa) B_u - (\xi + \eta + g) R_u \\
#' dQ_u/dt &= \eta R_u - (\nu + \theta) Q_u - \Omega_q \cdot Q_u \\
#' dS_u/dt &= \omega\Lambda + \zeta B_u + (\nu+\theta)(1-\delta) Q_u - \rho S_u - \Omega_s \cdot S_u \\
#' \end{array}
#' }
#' with analogous equations for the infected (**y**) and infectious (**z**) classes,
#' where the term \eqn{\tau} governs progression between infection stages.
#'
#' @name BRQS
#' @rdname BRQS
NULL

#' @title Compute derivatives for `BRQS` (**MY**)
#' @description Implements [dMYdt] for the BRQS ODE model.
#' @inheritParams ramp.xds::dMYdt
#' @return a [numeric] vector
#' @keywords internal
#' @export
dMYdt.BRQS <- function(t, y, xds_obj, s){

  Lambda = as.vector(xds_obj$Lambda[[s]])
  kappa = as.vector(xds_obj$kappa[[s]])

  with(get_MY_vars(y, xds_obj, s),{
    with(xds_obj$MY_obj[[s]],{

      dBudt <- (1-omega)*Lambda + rho*Su + (nu+theta)*delta*Qu + xi*R - (f+zeta)*Bu - Omega_b %*% Bu
      dRudt <- f*(1-q*kappa)*Bu - (xi + eta + g)*Ru
      dQudt <- eta*Ru - (nu + theta)*Qu - Omega_q %*% Qu
      dSudt <- omega*Lambda + zeta*Bu + (nu+theta)*(1-delta)*Qu - rho*Su - Omega_s %*% Su

      dBydt <- rho*Sy + (nu + theta)*delta*Qy + xi*R - (f+zeta)*By - Omega_b %*% By - tau*By
      dRydt <- f*q*kappa*Bu + f*By - (xi + eta + g)*Ry - tau*Ry
      dQydt <- eta*Ry - (nu + theta)*Qy - Omega_q %*% Qy - tau*Qy
      dSydt <- zeta*By + (nu+theta)*(1-delta)*Qy - rho*Sy - Omega_s%*%Sy - tau*Sy

      dBzdt <- rho*Sz + (nu + theta)*delta*Qz + xi*R - (f+zeta)*Bz - Omega_b %*% Bz + tau*By
      dRzdt <- f*Bz - (xi + eta + g)*Rz + tau*Ry
      dQzdt <- eta*Rz - (nu + theta)*Qz - Omega_q %*% Qz + tau*Qy
      dSzdt <- zeta*Bz + (nu+theta)*(1-delta)*Qz -rho*Sz -Omega_s %*% Sz + tau*Sy

      return(c(dSudt, dBudt, dRudt, dQudt,
               dSydt, dBydt, dRydt, dQydt,
               dSzdt, dBzdt, dRzdt, dQzdt))
    })
  })
}

#' @title Return the variables as a list
#' @description This method dispatches on the type of `xds_obj$MY_obj[[s]]`
#' @inheritParams ramp.xds::get_MY_vars
#' @return a [list]
#' @keywords internal
#' @export
get_MY_vars.BRQS<- function(y, xds_obj, s){
  with(xds_obj$MY_obj[[s]]$ix, return(list(
    Su <- y[Su_ix],
    Bu <- y[Bu_ix],
    Ru <- y[Ru_ix],
    Qu <- y[Qu_ix],
    Sy <- y[Sy_ix],
    By <- y[By_ix],
    Ry <- y[Ry_ix],
    Qy <- y[Qy_ix],
    Sz <- y[Sz_ix],
    Bz <- y[Bz_ix],
    Rz <- y[Rz_ix],
    Qz <- y[Qz_ix],
    B = Bu+By+Bz,
    R = Ru+Ry+Rz,
    Q = Qu+Qy+Qz,
    S = Su+Sy+Sz,
    M = B+R+Q+S
)))}

#' @title Mosquito bionomics for `BRQS` (**MY**)
#' @description Implements [MBionomics] for the BRQS model
#' @inheritParams ramp.xds::MBionomics
#' @return an **`xds`** object
#' @keywords internal
#' @export
MBionomics.BRQS <- function(t, y, xds_obj, s) {with(xds_obj$MY_obj[[s]],{
  xds_obj$MY_obj[[s]]$es_f       <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_q       <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_nu      <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_g_s     <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_g_b     <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_g_r     <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_g_q     <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_sigma_b <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_sigma_q <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_sigma_s <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_mu_b    <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_mu_q    <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_mu_s    <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_theta   <- rep(1, nPatches)
  return(pars)
})}


#' @title Apply effect sizes for `BRQS` (**MY**)
#' @description Implements [MEffectSizes] for the BRQS model
#' @inheritParams ramp.xds::MEffectSizes
#' @return an **`xds`** object
#' @keywords internal
#' @export
MEffectSizes.BRQS <- function(t, y, xds_obj, s) {with(xds_obj$MY_obj[[s]],{
  xds_obj$MY_obj[[s]]$f <- es_f*f_t
  xds_obj$MY_obj[[s]]$q <- es_q*q_t
  g <- es_g*g_t
  sigma_b <- es_sigma_b*sigma_b_t
  sigma_q <- es_sigma_q*sigma_q_t
  xds_obj$MY_obj[[s]]$g <- g
  xds_obj$MY_obj[[s]]$sigma_b <- sigma_b
  xds_obj$MY_obj[[s]]$sigma_q <- sigma_q
  xds_obj$MY_obj[[s]]$Omega_b = compute_Omega_xde(g, sigma_b, mu, calKb)
  xds_obj$MY_obj[[s]]$Omega_q = compute_Omega_xde(g, sigma_q, mu, calKq)
  return(pars)
})}

#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqZ] for the BRQS model.
#' @inheritParams ramp.xds::F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
F_fqZ.BRQS <- function(t, y, xds_obj, s) {
  f = get_f(xds_obj, s)
  q = get_f(xds_obj, s)
  fqZ = f*q*y[xds_obj$ix$MY[[s]]$Z_b_ix]
  return(fqZ)
}

#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqM] for the BRQS model.
#' @inheritParams ramp.xds::F_fqM
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
F_fqM.BRQS <- function(t, y, xds_obj, s) {
  f = get_f(xds_obj, s)
  q = get_f(xds_obj, s)
  M = with(xds_obj$ix$MY[[s]], y[Bu_ix] + y[By_ix] + y[Bz_ix])
  fqM = f*q*M
  return(fqM)
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the BRQS model.
#' @inheritParams ramp.xds::F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @keywords internal
#' @export
F_eggs.BRQS <- function(t, y, xds_obj, s) {

  G <- with(xds_obj$ix$MY[[s]], y[Qu_ix] + y[Qy_ix] + y[Qz_ix])
  eggs = with(xds_obj$MY_obj[[s]], {
    return(G*nu*eggsPerBatch)
  })
  return(eggs)
}


#' @title Setup MY_obj for the BRQS model
#' @description Implements [setup_MY_obj] for the RM model
#' @inheritParams ramp.xds::setup_MY_obj
#' @return a [list] vector
#' @keywords internal
#' @export
setup_MY_obj.BRQS = function(MYname, xds_obj, s, options=list()){
  xds_obj$MY_obj[[s]] = make_MY_obj_BRQS(xds_obj$nPatches, options)
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
make_MY_obj_BRQS = function(nPatches, options=list(), eip=12,
                          g=1/12, sigma_b=1/8, sigma_q=1/8, mu=0, f=0.5, q=0.95,
                          nu=1, eggsPerBatch=60){

  with(options,{
    MY_obj <- list()
    class(MY_obj) <- "BRQS"

    MY_obj$nPatches <- nPatches

    eip_par <- list()
    class(eip_par) <- 'static'
    MY_obj$eip_par <- eip_par
    MY_obj$eip     <- eip

    MY_obj$f_t          <- checkIt(f, nPatches)
    MY_obj$es_f         <- rep(1, nPatches)

    MY_obj$nu_t         <- checkIt(nu, nPatches)
    MY_obj$es_nu        <- rep(1, nPatches)

    MY_obj$p            <- checkIt(p, nPatches)
    MY_obj$zeta         <- checkIt(zeta, nPatches)

    MY_obj$theta_t        <- checkIt(theta, nPatches)
    MY_obj$es_theta        <- rep(1, nPatches)

    MY_obj$g_b_t         <- checkIt(g_b, nPatches)
    MY_obj$es_g_b        <- rep(1, nPatches)

    MY_obj$g_q_t          <- checkIt(g_q, nPatches)
    MY_obj$es_g_q         <- rep(1, nPatches)

    MY_obj$g_s_t          <- checkIt(g_s, nPatches)
    MY_obj$es_g_s         <- rep(1, nPatches)

    MY_obj$g_r_t          <- checkIt(g_r, nPatches)
    MY_obj$es_g_r         <- rep(1, nPatches)

    MY_obj$sigma_b_t    <- checkIt(sigma_b, nPatches)
    MY_obj$es_sigma_b   <- rep(1, nPatches)

    MY_obj$sigma_q_t    <- checkIt(sigma_q, nPatches)
    MY_obj$es_sigma_q   <- rep(1, nPatches)

    MY_obj$sigma_s_t    <- checkIt(sigma_q, nPatches)
    MY_obj$es_sigma_s   <- rep(1, nPatches)

    MY_obj$mu_b_t    <- checkIt(mu_b, nPatches)
    MY_obj$es_mu_b   <- rep(1, nPatches)

    MY_obj$mu_q_t    <- checkIt(mu_q, nPatches)
    MY_obj$es_mu_q   <- rep(1, nPatches)

    MY_obj$mu_s_t    <- checkIt(mu_q, nPatches)
    MY_obj$es_mu_s   <- rep(1, nPatches)

    MY_obj$xi_t           <- checkIt(xi, nPatches)
    MY_obj$es_xi         <- rep(1, nPatches)

    MY_obj$eta_t           <- checkIt(eta, nPatches)
    MY_obj$es_eta         <- rep(1, nPatches)

    MY_obj$q_t          <- checkIt(q, nPatches)
    MY_obj$es_q         <- rep(1, nPatches)

    MY_obj$eggsPerBatch <- eggsPerBatch

    MY_obj$phi <- 1/MY_obj$eip

    calK <- diag(nPatches)

    MY_obj$calKb <- calK
    MY_obj$calKq <- calK
    MY_obj$calKs <- calK

    Omega_par <- list()
    class(Omega_par) <- "static"
    MY_obj$Omega_par <- Omega_par
    MY_obj$Omega_b <- with(MY_obj, compute_Omega_xde(g, sigma_b, mu, calK))
    MY_obj$Omega_q <- with(MY_obj, compute_Omega_xde(g, sigma_q, mu, calK))
    MY_obj$Omega_s <- with(MY_obj, compute_Omega_xde(g, sigma_s, mu, calK))
    base <- 'BRQS'
    class(base) <- 'BRQS'
    MY_obj$baseline <- base

    return(MY_obj)
})}

#' @title Return the parameters as a list
#' @description This method dispatches on the type of `xds_obj$MY_obj[[s]]`.
#' @param xds_obj an **`xds`** object
#' @param s the vector species index
#' @return a [list]
#' @keywords internal
#' @export
get_MY_pars.BRQS <- function(xds_obj, s=1) {
  with(xds_obj$MY_obj[[s]], list(
    f=f_t, q=q_t, g=g_t, sigma_b=sigma_b_t, sigma_q = sigma_q_t, eip=eip, mu=mu_t,
    nu=nu_t, eggsPerBatch=eggsPerBatch, calK=calK
  ))
}

#' @title Return the parameters as a list
#' @description This method dispatches on the type of `xds_obj$MY_obj[[s]]`.
#' @inheritParams ramp.xds::change_MY_pars
#' @return an **`xds`** object
#' @keywords internal
#' @export
change_MY_pars.BRQS <- function(xds_obj, s=1, options=list()) {
  nHabitats <- xds_obj$nHabitats
  with(xds_obj$MY_obj[[s]], with(options,{
    xds_obj$MY_obj[[s]]$f_t = f
    xds_obj$MY_obj[[s]]$q_t = q
    xds_obj$MY_obj[[s]]$nu_t = nu
    xds_obj$MY_obj[[s]]$g_s_t = g_s
    xds_obj$MY_obj[[s]]$g_b_t = g_b
    xds_obj$MY_obj[[s]]$g_r_t = g_r
    xds_obj$MY_obj[[s]]$g_q_t = g_q
    xds_obj$MY_obj[[s]]$mu_s_t = mu_s
    xds_obj$MY_obj[[s]]$mu_b_t = mu_b
    xds_obj$MY_obj[[s]]$mu_q_t = mu_q
    xds_obj$MY_obj[[s]]$sigma_b_t = sigma_b
    xds_obj$MY_obj[[s]]$sigma_q_t = sigma_q
    xds_obj$MY_obj[[s]]$sigma_s_t = sigma_s
    xds_obj$MY_obj[[s]]$eip = eip
    xds_obj$MY_obj[[s]]$eggsPerBatch = eggsPerBatch
    return(pars)
  }))}

#' @title Setup initial values for the BRQS model
#' @description Implements [setup_MY_inits] for the RM model
#' @inheritParams ramp.xds::setup_MY_inits
#' @return a [list]
#' @keywords internal
#' @export
setup_MY_inits.BRQS = function(xds_obj, s, options=list()){
  xds_obj$MYinits[[s]] = with(xds_obj$MY_obj[[s]], make_MY_inits_BRQS(nPatches, options))
  return(xds_obj)
}

#' @title Make inits for BRQS adult mosquito model
#' @param nPatches the number of patches in the model
#' @param options a [list] of values that overwrites the defaults
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
#' @keywords internal
#' @export
make_MY_inits_BRQS = function(nPatches, options = list(),
                              Su=100, Bu=100, Ru=1, Qu=10,
                              Sy=0, By=0, Ry=0, Qy=0,
                              Sz=0, Bz=0, Rz=0, Qz=0){

  with(options,{
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
#' @description Implements [setup_MY_ix] for the BRQS model.
#' @inheritParams ramp.xds::setup_MY_ix
#' @return none
#' @importFrom utils tail
#' @keywords internal
#' @export
setup_MY_ix.BRQS <- function(xds_obj, s) {with(xds_obj,{

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

  xds_obj$max_ix = max_ix
  xds_obj$MY[[s]]$ix = list(
    Su_ix = Su_ix, Bu_ix = Bu_ix, Ru_ix=Ru_ix, Qu_ix = Qu_ix,
    Sy_ix = Sy_ix, By_ix = By_ix, Ry_ix=Ry_ix, Qy_ix = Qy_ix,
    Sz_ix = Sz_ix, Bz_ix = Bz_ix, Rz_ix=Rz_ix, Qz_ix = Qz_ix)


  return(pars)
})}


#' @title Parse the output of deSolve and return variables for the BRQS model
#' @description Implements [parse_MY_orbits] for the BRQS model
#' @inheritParams ramp.xds::parse_MY_orbits
#' @return none
#' @keywords internal
#' @export
parse_MY_orbits.BRQS <- function(outputs, xds_obj, s) {with(xds_obj$ix$MY[[s]],{
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

#' @title Make inits for BRQS adult mosquito model
#' @inheritParams ramp.xds::change_MY_inits
#' @return none
#' @keywords internal
#' @export
change_MY_inits.BRQS <- function(xds_obj, s, options = list()) {
  inits = with(get_MY_inits(xds_obj, s), with(options,
    list(Su=Su, Bu=Bu, Ru=Ru, Qu=Qu,
         Sy=Sy, By=By, Ry=Ry, Qy=Qy,
         Sz=Sz, Bz=Bz, Rz=Rz, Qz=Qz)))

  xds_obj$MY_obj[[s]]$inits=inits
  return(xds_obj)
}


#' @title Get the feeding rate
#' @param xds_obj an **`xds`** object
#' @param s the vector species index
#' @return a [numeric] vector
#' @keywords internal
#' @export
get_f.BRQS = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], f_t*es_f)
}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @keywords internal
#' @export
get_q.BRQS = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], q_t*es_q)
}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @keywords internal
#' @export
get_g.BRQS = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], g_t*es_g)
}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @keywords internal
#' @export
get_sigma.BRQS = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], sigma_b_t*es_sigma_b)
}
