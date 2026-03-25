# specialized methods for the adult mosquito RMG model

#' @title The **RMG** Module Skill Set
#'
#' @description The **MY** skill set is a list of
#' an module's capabilities:
#'
#' + `demography` is
#'
#' @inheritParams ramp.xds::skill_set_MY
#'
#' @return *MY* module skill set, as a list
#'
#' @export
skill_set_MY.RMG = function(MYname){
  return(list())
}

#' Run a check before solving
#'
#' @inheritParams ramp.xds::check_MY
#'
#' @returns an **`xds`** model object
#' @export
check_MY.RMG = function(xds_obj, s){
  return(xds_obj)
}

#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MBaseline] for the RMG model
#' @inheritParams ramp.xds::MBaseline
#' @return a named [list]
#' @export
MBaseline.RMG <- function(t, y, xds_obj, s) {with(xds_obj$MY_obj[[s]],{
  xds_obj$MY_obj[[s]]$es_g       <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_sigma_b <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_sigma_q <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_f       <- rep(1, nPatches)
  xds_obj$MY_obj[[s]]$es_q       <- rep(1, nPatches)
  return(xds_obj)
})}


#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MBionomics] for the RMG model
#' @inheritParams ramp.xds::MBionomics
#' @return a named [list]
#' @export
MBionomics.RMG <- function(t, y, xds_obj, s) {with(xds_obj$MY_obj[[s]],{
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
#' @description Implements [F_fqZ] for the RMG model.
#' @inheritParams ramp.xds::F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ.RMG <- function(t, y, xds_obj, s) {
  f = get_f(xds_obj, s)
  q = get_f(xds_obj, s)
  fqZ = f*q*y[xds_obj$MY_obj[[s]]$ix$Z_b_ix]
  return(fqZ)
}

#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqM] for the RMG model.
#' @inheritParams ramp.xds::F_fqM
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqM.RMG <- function(t, y, xds_obj, s) {
  f = get_f(xds_obj, s)
  q = get_f(xds_obj, s)
  M = with(xds_obj$MY_obj[[s]]$ix, y[U_b_ix] + y[Y_b_ix] + y[Z_b_ix])
  fqM = f*q*M
  return(fqM)
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the RMG model.
#' @inheritParams ramp.xds::F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.RMG <- function(t, y, xds_obj, s) {

  G <- with(get_MY_vars(y, xds_obj, s), U_g + Y_g + Z_g)
  eggs = with(xds_obj$MY_obj[[s]], {
    return(G*nu*eggsPerBatch)
  })
  return(eggs)
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYdt] for the RMG ODE model.
#' @inheritParams ramp.xds::dMYdt
#' @return a [numeric] vector
#' @export
dMYdt.RMG <- function(t, y, xds_obj, s){

  Lambda = xds_obj$terms$Lambda[[s]]
  kappa = xds_obj$terms$kappa[[s]]

  with(get_MY_vars(y, xds_obj, s),{
    with(xds_obj$MY_obj[[s]],{

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

#' @title Setup MY_obj for the RMG model
#' @description Implements [setup_MY_obj] for the RM model
#' @inheritParams ramp.xds::setup_MY_obj
#' @return a [list] vector
#' @export
setup_MY_obj.RMG = function(MYname, xds_obj, s, options=list()){
  xds_obj$MY_obj[[s]] = make_MY_obj_RMG(xds_obj$nPatches, options)
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
#' @export
make_MY_obj_RMG = function(nPatches, options=list(), eip=12,
                          g=1/12, sigma_b=1/8, sigma_q=1/8, mu=0, f=0.5, q=0.95,
                          nu=1, eggsPerBatch=60){

  with(options,{
    MY_obj <- list()
    class(MY_obj) <- "RMG"

    MY_obj$nPatches <- nPatches

    eip_par <- list()
    class(eip_par) <- 'static'
    MY_obj$eip_par <- eip_par
    MY_obj$eip     <- eip

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
    MY_obj$phi <- 1/MY_obj$eip

    calK <- diag(nPatches)

    MY_obj$calKb <- calK
    MY_obj$calKq <- calK

    Omega_par <- list()
    class(Omega_par) <- "static"
    MY_obj$Omega_par <- Omega_par
    MY_obj$Omega_b <- with(MY_obj, make_Omega_xde(g, sigma_b, mu, calK))
    MY_obj$Omega_q <- with(MY_obj, make_Omega_xde(g, sigma_q, mu, calK))
    base <- 'RMG'
    class(base) <- 'RMG'
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
#' @export
get_MY_pars.RMG <- function(xds_obj, s=1) {
  with(xds_obj$MY_obj[[s]], list(
    f=f_t, q=q_t, g=g_t, sigma=sigma_t, eip=eip, mu=mu_t,
    nu=nu_t, eggsPerBatch=eggsPerBatch, calK=calK
  ))
}


#' @title Return the parameters as a list
#' @description This method dispatches on the type of `xds_obj$MY_obj[[s]]`.
#' @inheritParams ramp.xds::change_MY_pars
#' @return an **`xds`** object
#' @export
change_MY_pars.RMG <- function(xds_obj, s=1, options=list()) {
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

#' @title Setup initial values for the RMG model
#' @description Implements [setup_MY_inits] for the RM model
#' @inheritParams ramp.xds::setup_MY_inits
#' @return a [list]
#' @export
setup_MY_inits.RMG = function(xds_obj, s, options=list()){
  xds_obj$MY_obj[[s]]$inits = make_MY_inits_RMG(xds_obj$nPatches, options)
  return(xds_obj)
}

#' @title Make inits for RMG adult mosquito model
#' @param nPatches the number of patches in the model
#' @param options a [list] of values that overwrites the defaults
#' @param U_b0 total uninfected, not gravid mosquito density at each patch
#' @param U_g0 total uninfected, gravid mosquito density at each patch
#' @param Y_b0 infected, not gravid mosquito density at each patch
#' @param Y_g0 infected, gravid mosquito density at each patch
#' @param Z_b0 infectious, not gravid mosquito density at each patch
#' @param Z_g0 infectious, gravid mosquito density at each patch
#' @return a [list]
#' @export
make_MY_inits_RMG = function(nPatches, options = list(),
                             U_b0=5, U_g0=1, Y_b0=1, Y_g0=1, Z_b0=1, Z_g0=1){
  with(options,{
    U_b = checkIt(U_b0, nPatches)
    U_g = checkIt(U_g0, nPatches)
    Y_b = checkIt(Y_b0, nPatches)
    Y_g = checkIt(Y_g0, nPatches)
    Z_b = checkIt(Z_b0, nPatches)
    Z_g = checkIt(Z_g0, nPatches)

    return(list(U_b=U_b, U_g=U_g, Y_b=Y_b, Y_g=Y_g, Z_b=Z_b, Z_g=Z_g))
  })}

#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [setup_MY_ix] for the RMG model.
#' @inheritParams ramp.xds::setup_MY_ix
#' @return none
#' @importFrom utils tail
#' @export
setup_MY_ix.RMG <- function(xds_obj, s) {with(xds_obj,{

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

  xds_obj$max_ix = max_ix
  xds_obj$MY_obj[[s]]$ix = list(U_b_ix=U_b_ix, U_g_ix=U_g_ix,
                          Y_b_ix=Y_b_ix, Y_g_ix=Y_g_ix,
                          Z_b_ix=Z_b_ix, Z_g_ix=Z_g_ix)
  return(xds_obj)
})}

#' @title Return the variables as a list
#' @description This method dispatches on the type of `xds_obj$MY_obj[[s]]`
#' @inheritParams ramp.xds::get_MY_vars
#' @return a [list]
#' @export
get_MY_vars.RMG <- function(y, xds_obj, s){
  with(xds_obj$MY_obj[[s]]$ix,
       return(list(
         U_b = y[U_b_ix],
         U_g = y[U_g_ix],
         Y_b = y[Y_b_ix],
         Y_g = y[Y_g_ix],
         Z_b = y[Z_b_ix],
         Z_g = y[Z_g_ix]
)))}

#' @title parse the output of deSolve and return variables for the RMG model
#' @description Implements [parse_MY_orbits] for the RMG model
#' @inheritParams ramp.xds::parse_MY_orbits
#' @return none
#' @export
parse_MY_orbits.RMG <- function(outputs, xds_obj, s) {
  with(xds_obj$MY_obj[[s]]$ix,{
    U_b = outputs[,U_b_ix]
    U_g = outputs[,U_g_ix]
    Y_b = outputs[,Y_b_ix]
    Y_g = outputs[,Y_g_ix]
    Z_b = outputs[,Z_b_ix]
    Z_g = outputs[,Z_g_ix]
    M =  U_b+Y_b+Z_b+U_g+Y_g+Z_g
    Y =  Y_b+Y_g+Z_b+Z_g
    Z =  Z_b+Z_g
    y = Y/M
    z = Z/M
  return(list(U_b=U_b, U_g=U_g,
              Y_b=Y_b, Y_g=Y_g,
              Z_b=Z_b, Z_g=Z_g,
              M=M, Y=Y, Z=Z, y=y, z=z))
})}


#' @title Set new MY parameter values
#' @description This method dispatches on the type of `xds_obj$MY_obj[[s]]`.
#' @inheritParams ramp.xds::change_MY_inits
#' @return an `xds` object
#' @export
change_MY_inits.RMG <- function(xds_obj, s=1, options=list()) {
  with(xds_obj$MY_obj[[s]],
    with(options,{
      xds_obj$MY_obj[[s]]$inits$U_b = U_b
      xds_obj$MY_obj[[s]]$inits$U_g = U_g
      xds_obj$MY_obj[[s]]$inits$Y_b = Y_b
      xds_obj$MY_obj[[s]]$inits$Y_g = Y_g
      xds_obj$MY_obj[[s]]$inits$Z_b = Z_b
      xds_obj$MY_obj[[s]]$inits$Y_g = Y_g
      return(xds_obj)
  }))}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** object
#' @param s the vector species index
#' @return a [numeric] vector
#' @export
get_f.RMG = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], f_t*es_f)
}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_q.RMG = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], q_t*es_q)
}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_g.RMG = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], g_t*es_g)
}

#' @title Get the feeding rate
#' @param xds_obj an **`xds`** object
#' @param s the vector species index
#' @return y a [numeric] vector assigned the class "dynamic"
#' @export
get_sigma.RMG = function(xds_obj, s=1){
  with(xds_obj$MY_obj[[s]], sigma_b_t*es_sigma_b)
}

