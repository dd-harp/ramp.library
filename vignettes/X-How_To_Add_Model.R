## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_X] for the NEW model.
#' @inheritParams ramp.xde::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.NEW <- function(t, y, pars) {
  ########################
  # extract: 
  # VAR <- y[pars$ix$X$...] 
  ########################
  with(pars$Xpar, 
       ########################
       # compute: 
       # X <- ... F(VAR) 
       ########################
  )
  return(X)
}

## -----------------------------------------------------------------------------
#' @title Size of effective infectious human population
#' @description Implements [F_H] for the NEW model.
#' @inheritParams ramp.xde::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.NEW <- function(t, y, pars) {
  ########################
  # extract: 
  # VAR <- y[pars$ix$X$...] 
  ########################
  with(pars$Xpar, 
       ########################
       # compute: 
       # X <- ... F(VAR) 
       ########################
    )
  return(H)
}

## -----------------------------------------------------------------------------
#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SIS model.
#' @inheritParams ramp.xde::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.NEW <- function(t, y, pars) {
  with(pars$Xpar, 
       ########################
       # retrieve it: 
       # ...  
       ########################
  )
  
  #######################
  # return it: 
  # return(...)
  ########################
}

## -----------------------------------------------------------------------------
#' @title Derivatives for human population
#' @description Implements [dXdt] for the NEW model, no demography.
#' @inheritParams ramp.xde::dXdt
#' @return a [numeric] vector
#' @export
dXdt.NEW <- function(t, y, pars, FoI) {

  ###############################
  # get variables by name from y 
  # ...  
  ###############################

  with(pars$Xpar, {
    ###############################
    # Compute the derivatives 
    # 
    # dX <- FoI*(H - X) - r*X + dHdt(t, X, pars)
    # dH <- Births(t, H, pars) + dHdt(t, H, pars) 
    ###############################
    
    ###############################
    # Return the derivatives 
    # ...  
    ###############################
    return(c(dX1, dX2, dX3, dH))
  })
}

## ----eval=F-------------------------------------------------------------------
#  dXdt.SISdH <- function(t, y, pars, FoI) {
#    with(pars$Xpar, {
#  
#      H <- F_H(t, y, pars)
#      X <- y[X_ix]
#  
#      dX <- FoI*(H - X) - r*X + dHdt(t, X, pars)
#      dH <- Births(t, H, pars) + dHdt(t, H, pars)
#  
#      return(c(dX, dH))
#    })
#  }

## -----------------------------------------------------------------------------
#' @title Setup Xpar.NEW
#' @description Implements [setup_X] for the NEW model
#' @inheritParams ramp.xde::setup_X
#' @return a [list] vector
#' @export
setup_X.NEW = function(pars, Xname, Xopts=list()){

  pars$Xname = "NEW"
  pars = make_Xpar_NEW(pars, Xopts)
  pars = make_Xinits_NEW(pars, Xopts)

  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Make parameters for NEW human model, with defaults
#' @param pars a [list]
#' @param Xopts a [list] that could overwrite defaults
#' @param p1 the first parameter 
#' @param p2 the second parameter 
#' @param p3 the third parameter 
#' @return a [list]
#' @export
make_Xpar_NEW = function(pars, Xopts=list(),
                         p1=1, p2=2, p3=3){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("NEW")

    Xpar$p1 = checkIt(p1, pars$nStrata)
    Xpar$p2 = checkIt(p2, pars$nStrata)
    Xpar$p3 = checkIt(p3, pars$nStrata)

    pars$Xpar = Xpar
    return(pars)
  })}

## -----------------------------------------------------------------------------
#' @title Make initial values for the NEW human model, with defaults
#' @param pars a [list]
#' @param Xopts a [list] to overwrite defaults
#' @param X10 the initial values of the parameter X1
#' @param X20 the initial values of the parameter X2
#' @param X30 the initial values of the parameter X3
#' @return a [list]
#' @export
make_Xinits_NEW = function(pars, Xopts = list(), X10=1, X20=2, X30=3){with(Xopts,{
  inits = list()
  inits$X10 = checkIt(X10, pars$nStrata)
  inits$X20 = checkIt(X20, pars$nStrata)
  inits$X30 = checkIt(X30, pars$nStrata)
  pars$Xinits = inits
  return(pars)
})}

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the NEW model.
#' @inheritParams ramp.xde::make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.NEW <- function(pars) {
  
  pars$Xpar$X1_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$Xpar$X1_ix, 1)
  
  pars$Xpar$X2_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$Xpar$X2_ix, 1)
  
  pars$Xpar$X3_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$Xpar$X3_ix, 1)
  
  return(pars)
}

## -----------------------------------------------------------------------------
#' @title Update inits for the SIS human model from a vector of states
#' @inheritParams ramp.xde::update_inits_X 
#' @return none
#' @export
update_inits_X.NEW <- function(pars, y0) {
  with(pars$ix$X,{
    X10 = y0[X1_ix]
    X20 = y0[X2_ix]
    X30 = y0[X3_ix]
    pars = make_Xinits_SIS(pars, list(), X10, X20, X30)
    return(pars)
})}

## -----------------------------------------------------------------------------
#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars a [list]
#' @return none
#' @export
get_inits_X.NEW <- function(pars){
  with(pars$Xinits,{
    return(c(X1, X2, X3)) 
  })
}

## -----------------------------------------------------------------------------
#' @title Parse the output of deSolve and return variables for the NEW model
#' @description Implements [parse_deout_X] for the NEW model
#' @inheritParams ramp.xde::parse_deout_X
#' @return none
#' @export
parse_deout_X.NEW <- function(deout, pars) {
  time = deout[,1]
  Hlist <- parse_deout_H(deout, pars)
  with(Hlist,{
    X1= deout[,pars$ix$X$X1_ix+1]
    X2= deout[,pars$ix$X$X1_ix+1]
    X3= deout[,pars$ix$X$X3_ix+1]
    return(list(time=time, X1=X1, X2=X2, X3=X3, H=H))
})}

## -----------------------------------------------------------------------------
#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the NEW model.
#' @inheritParams ramp.xde::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.NEW<- function(varslist, pars) {
  pr = with(varslist$XH, X1/H)
  return(pr)
}

## -----------------------------------------------------------------------------
#' @title Compute the HTC for the NEW model
#' @description Implements [HTC] for the NEW model with demography.
#' @inheritParams ramp.xde::HTC
#' @return a [numeric] vector
#' @export
HTC.NEW <- function(pars) {
  with(pars$Xpar,
    return(c/r)
  )
}

## -----------------------------------------------------------------------------
#' Plot the density of infected individuals for the NEW model
#'
#' @inheritParams ramp.xde::xde_plot_X
#' @export
xde_plot_X.NEW = function(pars, clrs="black", llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH,
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "Time"))

  xde_lines_X(vars$XH, pars, clrs, llty)
}

## -----------------------------------------------------------------------------
#' Add lines for the density of infected individuals for the NEW model
#'
#' @inheritParams ramp.xde::xde_lines_X
#'
#' @export
xde_lines_X.NEW = function(XH, pars, clrs="black", llty=1){
  with(XH,{
    if(pars$nStrata==1) lines(time, X1, col=clrs[1], lty = llty[1])
    if(pars$nStrata>1){
      if (length(clrs)==1) clrs=rep(clrs, pars$nStrata)
      if (length(llty)==1) llty=rep(llty, pars$nStrata)
      for(i in 1:pars$nStrata){
        lines(time, X1[,i], col=clrs[i], lty = llty[i])
      }
    }
  })}

