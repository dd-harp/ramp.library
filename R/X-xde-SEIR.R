#' @title Compute the derivatives for parasite infection dynamics in human population strata
#' @description Implements [dXdt] for the SEIR model
#' @inheritParams ramp.xds::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SEIR<- function(t, y, pars, i) {
  # do not change this
  foi <- pars$FoI[[i]]
  Hpar <- pars$Hpar[[i]]

  # attach the variables by name
  with(list_Xvars(y, pars, i),{
    # compute H (if it isn't one of the variables)

    # expose the parameters (see create_Xpar_SEIR)
    with(pars$Xpar[[i]], {
      # compute the derivatives
      dS <- Births(t, H, Hpar) - foi*S + dHdt(t, S, Hpar)
      dE <- foi*S - tau*E + dHdt(t, E, Hpar)
      dI <- tau*E - r*I + dHdt(t, I, Hpar)
      dR <- r*I + dHdt(t, R, Hpar)

       # concatenate the derivatives
      derivs = c(dS, dE, dI, dR)

      # return the derivatives
      return(derivs)
    })
  })
}



#' @title Make initial values for the SEIR human model, with defaults
#' @param nStrata the number of strata in the model
#' @param Xopts a [list] to overwrite defaults
#' @param H0 the initial value for H
#' @param S the initial value for S
#' @param E the initial value for E
#' @param I the initial value for I
#' @param R the initial values for R
#' @return a [list]
#' @export
create_Xinits_SEIR = function(nStrata, Xopts = list(), H0= NULL, S=NULL, I=1, E=0,R = 1){with(Xopts,{
  if(is.null(S)) S = H0-(E+I+R)
  stopifnot(is.numeric(S))
  S = checkIt(S, nStrata)
  E = checkIt(E, nStrata)
  I = checkIt(I, nStrata)
  R = checkIt(R, nStrata)
  return(list(S=S, E=E, I=I, R =R))
})}






#' @title Setup Xinits.SEIR
#' @description Implements [make_Xinits] for the SEIR model
#' @inheritParams ramp.xds::make_Xinits
#' @return a [list] vector
#' @export
make_Xinits.SEIR = function(pars, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars, create_Xinits_SEIR(pars$nStrata[i], Xopts, H0=Hpar[[i]]$H))
  return(pars)
}





#' @title Add indices for human population to parameter list
#' @description Implements [make_X_indices] for the SEIR model.
#' @inheritParams ramp.xds::make_X_indices
#' @return none
#' @importFrom utils tail
#' @export
make_X_indices.SEIR <- function(pars, i) {with(pars,{

  S_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(S_ix, 1)

  E_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(E_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(I_ix, 1)

  R_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(R_ix, 1)



  pars$max_ix = max_ix
  pars$ix$X[[i]] = list(S_ix=S_ix,  E_ix=E_ix, I_ix=I_ix, R_ix=R_ix)
  return(pars)
})}





#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xds::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SEIR <- function(y, pars, i) {
  with(pars$ix$X[[i]],{
       S = y[S_ix]
       E = y[E_ix]
       I = y[I_ix]
       R = y[R_ix]
       H =S+E+I+R
       return(list(S=S,E=E,I=I,R=R,H=H))})
}





#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @inheritParams ramp.xds::get_Xinits
#' @return a [numeric] vector
#' @export
get_Xinits.SEIR <- function(pars, i){
  pars$Xinits[[i]]
}






#' @title Update inits for the SEIR human model from a vector of states
#' @inheritParams ramp.xds::update_Xinits
#' @return none
#' @export
update_Xinits.SEIR <- function(pars, y0, i) {
  with(list_Xvars(y0, pars, i),{
    pars = create_Xinits_SEIR(pars, list(), S=S,E= E, I=I, R=R)
    return(pars)
  })}



#' @title Make parameters for SEIR human model, with defaults
#' @param nStrata is the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param tau  incubation rate
#' @param b the proportion of infective bites that cause an infection
#' @param r the the duration of an infection
#' @param c the proportion of bites on infected humans that infect a mosquito
#' @return a [list]
#' @export
create_Xpar_SEIR = function(nStrata, Xopts=list(),
                          b=0.55, r=1/180, c=0.15,tau= 0.5){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SEIR")

    Xpar$b = checkIt(b, nStrata)
    Xpar$tau = checkIt(tau, nStrata)
    Xpar$c = checkIt(c, nStrata)
    Xpar$r = checkIt(r, nStrata)

    return(Xpar)
  })}





#' @title Setup Xpar.SEIR
#' @description Implements [make_Xpar] for the SEIR model
#' @inheritParams ramp.xds::make_Xpar
#' @return a [list] vector
#' @export
make_Xpar.SEIR = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = create_Xpar_SEIR(pars$nStrata[i], Xopts)
  return(pars)
}





#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIS model.
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SEIR <- function(y, pars, i) {
  I = y[pars$ix$X[[i]]$I_ix]
  Y = with(pars$Xpar[[i]], c*I)
  return(Y)
}


#' @title Size of effective infectious human population
#' @description Implements [F_H] for the SEIR model.
#' @inheritParams ramp.xds::F_H
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SEIR <- function(y, pars, i){
  with(list_Xvars(y, pars, i), return(H))
}





#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SEIR model.
#' @inheritParams ramp.xds::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SEIR <- function(y, pars, i) {
  with(pars$Xpar[[i]],return(b))
}


#' @title Parse the output of deSolve and return variables for the SEIR model
#' @description Implements [parse_outputs_X] for the SEIR model
#' @inheritParams ramp.xds::parse_outputs_X
#' @return none
#' @export
parse_outputs_X.SEIR <- function(outputs, pars, i) {
  time = outputs[,1]
  with(pars$ix$X[[i]],{
    S = outputs[,S_ix+1]
    E = outputs[,E_ix+1]
    I = outputs[,I_ix+1]
    R = outputs[,R_ix+1]

    H = S+I+R
    return(list(time=time, S=S, E=E,I=I, R=R, H=H))
  })}





#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIR model.
#' @inheritParams ramp.xds::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIR <- function(vars, Xpar) {
  pr = with(vars, I/H)
  return(pr)
}





#' @title Compute the HTC for the SEIR model
#' @description Implements [HTC] for the SEIR model with demography.
#' @inheritParams ramp.xds::HTC
#' @return a [numeric] vector
#' @export
HTC.SEIR <- function(pars, i) {
  with(pars$Xpar[[i]],
       HTC <- c/r,
       return(HTC)
  )
}





#' Add lines for the density of infected individuals for the SEIR model
#'
#' @param XH a list with the outputs of parse_outputs_X_SEIR
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
xde_lines_X_SEIR = function(XH, nStrata, clrs=c("black","darkblue","darkred","darkgreen"), llty=1){
  with(XH,{
    if(nStrata==1) {
      lines(time, S, col=clrs[1], lty = llty[1])
      lines(time, E, col=clrs[2], lty = llty[1])
      lines(time, I, col=clrs[3], lty = llty[1])
      lines(time, R, col=clrs[4], lty = llty[1])
    }
    if(nStrata>1){
      if (length(clrs)==3) clrs=matrix(clrs, 3, nStrata)
      if (length(llty)==1) llty=rep(llty, nStrata)

      for(i in 1:nStrata){
        lines(time, S[,i], col=clrs[1,i], lty = llty[i])
        lines(time, E[,i], col=clrs[2,i], lty = llty[i])
        lines(time, I[,i], col=clrs[3,i], lty = llty[i])
        lines(time, R[,i], col=clrs[4,i], lty = llty[i])
      }
    }
  })}






#' Plot the density of infected individuals for the SEIR model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @export
xds_plot_X.SEIR = function(pars, i=1, clrs=c("black","darkblue","darkred","darkgreen"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH[[i]],
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "No of. Infected", xlab = "Time"))


  xde_lines_X_SEIR(vars$XH[[i]], pars$nStrata[i], clrs, llty)
}
