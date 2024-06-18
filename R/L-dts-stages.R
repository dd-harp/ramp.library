# specialized methods for the aquatic mosquito stages_dts competition model

#' @title Reset aquatic parameters to baseline
#' @description Implements [LBionomics] for the stages_dts model
#' @inheritParams ramp.xde::LBionomics
#' @return a named [list]
#' @export
LBionomics.stages_dts <- function(t, y, pars, s){
  with(list_Lvars(y,pars,s),{
    with(pars$Lpar[[s]], {
      L = L1 + L2 + L3 + L4
      pars$Lpar[[s]]$surv1 = exp(-(phi1+theta1*L))
      pars$Lpar[[s]]$surv2 = exp(-(phi2+theta2*L))
      pars$Lpar[[s]]$surv3 = exp(-(phi3+theta3*L))
      pars$Lpar[[s]]$surv4 = exp(-(phi4+theta4*L))
      pars$Lpar[[s]]$mat1 =  exp(-psi1*exp(-xi1*L))
      pars$Lpar[[s]]$mat2 =  exp(-psi2*exp(-xi2*L))
      pars$Lpar[[s]]$mat3 =  exp(-psi3*exp(-xi3*L))
      pars$Lpar[[s]]$mat4 =  exp(-psi4*exp(-xi4*L))
      return(pars)
    })})}


#' @title Number of newly emerging adults from each larval habitat
#' @description Implements [F_alpha] for the stages_dts competition model.
#' @inheritParams ramp.xde::F_alpha
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_alpha.stages_dts <- function(t, y, pars, s){
  y[pars$ix$L[[s]]$P_ix]
}

#' @title Derivatives for aquatic stage mosquitoes
#' @description Implements [dLdt] for the stages_dts competition model.
#' @inheritParams ramp.xde::dLdt
#' @return a [numeric] vector
#' @export
dLdt.stages_dts <- function(t, y, pars, s) {
  eta <- pars$eggs_laid[[s]]
  with(list_Lvars(y,pars,s),{
    with(pars$Lpar[[s]], {
      L = L1 + L2 + L3 + L4
      surv1 = exp(-(phi1+theta1*L)); surv2 = exp(-(phi2+theta2*L))
      surv3 = exp(-(phi3+theta3*L)); surv4 = exp(-(phi4+theta4*L))
      mat1 =  exp(-psi1*exp(-xi1*L)); mat2 =  exp(-psi2*exp(-xi2*L))
      mat3 =  exp(-psi3*exp(-xi3*L)); mat4 =  exp(-psi4*exp(-xi4*L))

      L1t = eta + (1-mat1)*surv1*L1
      L2t = mat1*surv1*L1 + (1-mat2)*surv2*L2
      L3t = mat2*surv2*L2 + (1-mat3)*surv3*L3
      L4t = mat3*surv3*L3 + (1-mat4)*surv4*L4
      Pt  = mat4*surv4*L4
      return(L1t, L2t, L3t, L4t, Pt)
    })
  })
}

#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Lpar[[s]]`
#' @inheritParams ramp.xde::list_Lvars
#' @return a [list]
#' @export
list_Lvars.stages_dts <- function(y, pars, s){
  with(pars$ix$L[[s]],
       return(list(
         L1 = y[L1_ix],
         L2 = y[L2_ix],
         L3 = y[L3_ix],
         L4 = y[L4_ix],
         P = y[P_ix]
       )))
}

#' @title Setup Lpar for the stages_dts model
#' @description Implements [setup_Lpar] for the stages_dts model
#' @inheritParams ramp.xde::setup_Lpar
#' @return a [list] vector
#' @export
setup_Lpar.stages_dts = function(Lname, pars, s, Lopts=list()){
  pars$Lpar[[s]] = make_Lpar_stages_dts(pars$nHabitats, Lopts)
  return(pars)
}

#' @title Setup the stages_dts model
#' @description Implements [setup_Linits] for the stages_dts model
#' @inheritParams ramp.xde::setup_Linits
#' @return a [list]
#' @export
setup_Linits.stages_dts = function(pars, s, Lopts=list()){
  pars$Linits[[s]] = make_Linits_stages_dts(pars$nHabitats, Lopts)
  return(pars)
}

#' @title Make parameters for stages_dts competition aquatic mosquito model
#' @param nHabitats the number of habitats in the model
#' @param Lopts a [list] that overwrites default values
#' @param psi1 maturation rates for L1 in each aquatic habitat
#' @param psi2 maturation rates for L2 in each aquatic habitat
#' @param psi3 maturation rates for L3 in each aquatic habitat
#' @param psi4 maturation rates for L1 in each aquatic habitat
#' @param xi1 delayed maturation for L1 due to mean crowding
#' @param xi2 delayed maturation for L2 due to mean crowding
#' @param xi3 delayed maturation for L3 due to mean crowding
#' @param xi4 delayed maturation for L4 due to mean crowding
#' @param phi1 density-independent mortality rates for L1 for each aquatic habitat
#' @param phi2 density-independent mortality rates for L2 for each aquatic habitat
#' @param phi3 density-independent mortality rates for L3 for each aquatic habitat
#' @param phi4 density-independent mortality rates for L4 for each aquatic habitat
#' @param theta1 density-dependent mortality terms for L1 for each aquatic habitat
#' @param theta2 density-dependent mortality terms for L2 for each aquatic habitat
#' @param theta3 density-dependent mortality terms for L3 for each aquatic habitat
#' @param theta4 density-dependent mortality terms for L4 for each aquatic habitat
#' @return a [list] with Lpar added
#' @export
make_Lpar_stages_dts = function(nHabitats, Lopts=list(),
                                psi1=1/8,  psi2=1/8,  psi3=1/8,  psi4=1/8,
                                xi1=0,  xi2=0,  xi3=0,  xi4=0,
                                phi1=1/8,  phi2=1/8,  phi3=1/8,  phi4=1/8,
                                theta1=1/100, theta2=1/100, theta3=1/100, theta4=1/100
){with(Lopts,{
  Lpar = list()
  class(Lpar) <- "stages_dts"

  Lpar$psi1 = checkIt(psi1, nHabitats)
  Lpar$psi2 = checkIt(psi2, nHabitats)
  Lpar$psi3 = checkIt(psi3, nHabitats)
  Lpar$psi4 = checkIt(psi4, nHabitats)
  Lpar$xi1 = checkIt(xi1, nHabitats)
  Lpar$xi2 = checkIt(xi2, nHabitats)
  Lpar$xi3 = checkIt(xi3, nHabitats)
  Lpar$xi4 = checkIt(xi4, nHabitats)
  Lpar$psi1 = checkIt(psi1, nHabitats)
  Lpar$psi2 = checkIt(psi2, nHabitats)
  Lpar$psi3 = checkIt(psi3, nHabitats)
  Lpar$psi4 = checkIt(psi4, nHabitats)
  Lpar$theta1 = checkIt(theta1, nHabitats)
  Lpar$theta2 = checkIt(theta2, nHabitats)
  Lpar$theta3 = checkIt(theta3, nHabitats)
  Lpar$theta4 = checkIt(theta4, nHabitats)

  Lpar$psi10 = checkIt(psi1, nHabitats)
  Lpar$psi20 = checkIt(psi2, nHabitats)
  Lpar$psi30 = checkIt(psi3, nHabitats)
  Lpar$psi40 = checkIt(psi4, nHabitats)
  Lpar$xi10 = checkIt(xi1, nHabitats)
  Lpar$xi20 = checkIt(xi2, nHabitats)
  Lpar$xi30 = checkIt(xi3, nHabitats)
  Lpar$xi40 = checkIt(xi4, nHabitats)
  Lpar$psi10 = checkIt(psi1, nHabitats)
  Lpar$psi20 = checkIt(psi2, nHabitats)
  Lpar$psi30 = checkIt(psi3, nHabitats)
  Lpar$psi40 = checkIt(psi4, nHabitats)
  Lpar$theta10 = checkIt(theta1, nHabitats)
  Lpar$theta20 = checkIt(theta2, nHabitats)
  Lpar$theta30 = checkIt(theta3, nHabitats)
  Lpar$theta40 = checkIt(theta4, nHabitats)

  return(Lpar)
})}

#' @title Make inits for stages_dts competition aquatic mosquito model
#' @param nHabitats the number of habitats in the model
#' @param Lopts a [list] that overwrites default values
#' @param L10 initial conditions
#' @param L20 initial conditions
#' @param L30 initial conditions
#' @param L40 initial conditions
#' @param P0 initial conditions
#' @return a [list] with Linits added
#' @export
make_Linits_stages_dts = function(nHabitats, Lopts=list(), L10=1, L20=0, L30=0, L40=0, P0=0){with(Lopts,{
  L10 = checkIt(L10, nHabitats)
  L20 = checkIt(L20, nHabitats)
  L30 = checkIt(L30, nHabitats)
  L40 = checkIt(L40, nHabitats)
  P0 = checkIt(P0, nHabitats)
  return(list(L1=L10, L2=L20, L3=L30, L4=L40, P=P0))
})}

#' @title Add indices for aquatic stage mosquitoes to parameter list
#' @description Implements [make_indices_L] for stages_dts competition model.
#' @inheritParams ramp.xde::make_indices_L
#' @return none
#' @importFrom utils tail
#' @export
make_indices_L.stages_dts <- function(pars, s) {with(pars,{

  L1_ix <- seq(from = max_ix+1, length.out=nHabitats)
  max_ix <- tail(L1_ix, 1)

  L2_ix <- seq(from = max_ix+1, length.out=nHabitats)
  max_ix <- tail(L2_ix, 1)

  L3_ix <- seq(from = max_ix+1, length.out=nHabitats)
  max_ix <- tail(L3_ix, 1)

  L4_ix <- seq(from = max_ix+1, length.out=nHabitats)
  max_ix <- tail(L4_ix, 1)

  P_ix <- seq(from = max_ix+1, length.out=nHabitats)
  max_ix <- tail(P_ix, 1)

  pars$max_ix = max_ix
  pars$ix$L[[s]] = list(L1_ix=L1_ix, L2_ix=L2_ix, L3_ix=L3_ix, L4_ix=L4_ix, P_ix=P_ix)
  return(pars)
})}

#' @title Parse the variable names for the stages_dts model
#' @description Implements [parse_outputs_L] for stages_dts competition model.
#' @inheritParams ramp.xde::parse_outputs_L
#' @return [list]
#' @export
parse_outputs_L.stages_dts <- function(outputs, pars, s) {
  time = outputs[,1]
  L1 = outputs[,pars$ix$L[[s]]$L1_ix+1]
  L2 = outputs[,pars$ix$L[[s]]$L2_ix+1]
  L3 = outputs[,pars$ix$L[[s]]$L3_ix+1]
  L4 = outputs[,pars$ix$L[[s]]$L4_ix+1]
  P = outputs[,pars$ix$L[[s]]$P_ix+1]
  return(list(time=time, L1=L1, L2=L2, L3=L3, L4=L4, P=P))
}


#' @title Make parameters for stages_dts competition aquatic mosquito model
#' @param pars a [list]
#' @param psi1 maturation rates for L1 in each aquatic habitat
#' @param psi2 maturation rates for L2 in each aquatic habitat
#' @param psi3 maturation rates for L3 in each aquatic habitat
#' @param psi4 maturation rates for L1 in each aquatic habitat
#' @param xi1 delayed maturation for L1 due to mean crowding
#' @param xi2 delayed maturation for L2 due to mean crowding
#' @param xi3 delayed maturation for L3 due to mean crowding
#' @param xi4 delayed maturation for L4 due to mean crowding
#' @param phi1 density-independent mortality rates for L1 for each aquatic habitat
#' @param phi2 density-independent mortality rates for L2 for each aquatic habitat
#' @param phi3 density-independent mortality rates for L3 for each aquatic habitat
#' @param phi4 density-independent mortality rates for L4 for each aquatic habitat
#' @param theta1 density-dependent mortality terms for L1 for each aquatic habitat
#' @param theta2 density-dependent mortality terms for L2 for each aquatic habitat
#' @param theta3 density-dependent mortality terms for L3 for each aquatic habitat
#' @param theta4 density-dependent mortality terms for L4 for each aquatic habitat
#' @return a [list] with Lpar added
#' @export
make_parameters_L_stages_dts <- function(pars, psi1, psi2, psi3, psi4,
                                         xi1,  xi2,  xi3,  xi4,
                                         phi1,  phi2,  phi3,  phi4,
                                         theta1,  theta2,  theta3,  theta4) {
  stopifnot(is.numeric(psi1), is.numeric(phi1), is.numeric(xi1), is.numeric(theta1))
  Lpar <- list()
  class(Lpar) <- 'stages_dts'
  Lpar$psi1 <- psi1
  Lpar$psi2 <- psi2
  Lpar$psi3 <- psi3
  Lpar$psi4 <- psi4
  Lpar$xi1 <- xi1
  Lpar$xi2 <- xi2
  Lpar$xi3 <- xi3
  Lpar$xi4 <- xi4
  Lpar$phi1 <- phi1
  Lpar$phi2 <- phi2
  Lpar$phi3 <- phi3
  Lpar$phi4 <- phi4
  Lpar$theta1 <- theta1
  Lpar$theta2 <- theta2
  Lpar$theta3 <- theta3
  Lpar$theta4 <- theta4

  Lpar$psi10 <- psi1
  Lpar$psi20 <- psi2
  Lpar$psi30 <- psi3
  Lpar$psi40 <- psi4
  Lpar$xi10 <- xi1
  Lpar$xi20 <- xi2
  Lpar$xi30 <- xi3
  Lpar$xi40 <- xi4
  Lpar$phi10 <- phi1
  Lpar$phi20 <- phi2
  Lpar$phi30 <- phi3
  Lpar$phi40 <- phi4
  Lpar$theta10 <- theta1
  Lpar$theta20 <- theta2
  Lpar$theta30 <- theta3
  Lpar$theta40 <- theta4

  pars$Lpar[[1]] <- Lpar
  return(pars)
}

#' @title Make inits for stages_dts competition aquatic mosquito model
#' @param pars a [list]
#' @param L10 initial conditions
#' @param L20 initial conditions
#' @param L30 initial conditions
#' @param L40 initial conditions
#' @param P0 initial conditions
#' @return a [list] with Linits added
#' @export
make_inits_L_stages_dts <- function(pars, L10, L20, L30, L40, P0){
  stopifnot(is.numeric(L10))
  stopifnot(is.numeric(L20))
  stopifnot(is.numeric(L30))
  stopifnot(is.numeric(L40))
  stopifnot(is.numeric(P0))
  pars$Linits[[1]] <- list(L1=L10, L2=L20, L3=L30, L4=L40, P=P0)
  return(pars)
}

#' @title Update inits for the stages_dts aquatic mosquito competition model
#' @inheritParams ramp.xde::update_inits_L
#' @return none
#' @export
update_inits_L.stages_dts <- function(pars, y0, s) {
  with(list_Lvars(y0, pars, s),{
    pars = make_Linits_stages_dts(pars$nHabitats, L10=L1, L20=L2, L30=L3, L40=L4, P0=P)
    return(pars)
  })}

#' @title Return initial values as a vector
#' @description Implements [get_inits_L] for the GeRM model.
#' @inheritParams ramp.xde::get_inits_L
#' @return none
#' @export
get_inits_L.stages_dts <- function(pars, s){
  pars$Linits[[s]]$L
}

