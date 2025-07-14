# specialized methods for the human SIPmav3 model, formulated for differential and difference equations


#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIPmav3-xde model.
#' @inheritParams ramp.xds::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SIPmav3 <- function(t, y, pars, i){

  foi <- pars$FoI[[i]]
  Hpar <- pars$Hpar[[i]]
  Xpars = pars$Xpar[[i]]
  with(list_Xvars(y, pars, i),{
    with(Xpars, {

      x = I/H

      # Severe Malaria
      severe_h = F_incidence_h(y, Xpars, Xpars$severe_h)
      severe_x = F_incidence_x(y, Xpars, Xpars$severe_x)
      severe = severe_h*foi + severe_x*x

      # Moderate Malaria
      moderate_h = F_incidence_h(y, Xpars, Xpars$moderate_h)
      moderate_x = F_incidence_x(y, Xpars, Xpars$moderate_x)
      moderate = moderate_h*foi + moderate_x*x

      # Mild Malaria
      mild_h = F_incidence_h(y, Xpars, Xpars$mild_h)
      mild_x = F_incidence_x(y, Xpars, Xpars$mild_x)
      mild = mild_h*foi + mild_x*x

      # Treat Malaria
      pi3 = F_treat(y, Xpars, Xpars$treat3)
      rho = sum(pi3*c(severe_h, moderate_h, mild_h))
      sigma = sum(pi3*c(severe_x, moderate_x, mild_x))


      # Treat Malaria Infections
      rdt_pos = x*F_detect_y(y, Xpars, Xpars$rdt_mod)
      seek_care = F_force_malaria(t, y, Xpars, Xpars$seek_care_mod)
      ofever = F_force_malaria(t, y, Xpars, Xpars$other_fever_mod)
      sigma = sigma + F_force_malaria(t, y, Xpars, Xpars$msat_mod)*rdt_pos
      sigma = sigma + ofever*seek_care*rdt_pos

      # Background Drug Taking
      xi     = F_force_malaria(t, y, Xpars, Xpars$xi_mod)
      xi     = xi + F_force_malaria(t, y, Xpars, Xpars$mda_mod)
      xi     = xi + ofever*(1-seek_care)

      # Total Treatment Rate
      treat = rho*foi*(H-P) + sigma*I + xi*H

      r     = F_clear(y, Xpars)

      dH <- dHdt(t, H, Hpar) + Births(t, H, Hpar)
      dI <- foi*(1-rho)*S - (xi+sigma+r)*I + dHdt(t, I, Hpar)
      dP <- (foi*rho  + xi)*(H-P) + sigma*I - eta*P + dHdt(t, P, Hpar)

      dmoi <- foi*(1-rho)*(1-P/H) - (eta + foi*rho + sigma + xi)*moi + dAdt(t, moi, Hpar)

      daoi <- 1 - aoi*foi*(1-rho)*(1-P/H)/ifelse(moi==0, 1, moi) + dAdt(t, aoi, Hpar)

      dvh <- foi          + dAdt(t, vh, Hpar)

      return(c(dH, dI, dP, dmoi, daoi, dvh, severe, moderate, mild, treat))
    })
  })
}
#' @title Setup the Xpar for the SIPmav3_xde model
#' @description implements [setup_Xpar] for the SIPmav3 model
#' @inheritParams ramp.xds::setup_Xpar
#' @return a [list] vector
#' @export
setup_Xpar.SIPmav3 = function(Xname, pars, i, Xopts=list()){
  pars$Xpar[[i]] = make_Xpar_SIPmav3(pars$nStrata[i], Xopts)
  return(pars)
}

#' @title Make parameters for SIPmav3_xde human model, with defaults
#' @param nStrata the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param mod_b a string, the S3 class name for F_infect
#' @param b_opts a list: options for setup_F_infect
#' @param mod_c a string, the S3 class name for F_transmit
#' @param c_opts a list: options for setup_F_transmit
#' @param mod_r a string, the S3 class name for F_clear
#' @param r_opts a list: options for setup_F_clear
#' @param mod_lm a string, the S3 class name for F_detect for light microscopy
#' @param lm_opts a list: options for setup_F_detect
#' @param mod_pcr a string, the S3 class name for F_detect for PCR
#' @param pcr_opts a list: options for setup_F_detect
#' @param mod_rdt a string, the S3 class name for F_detect for RDT
#' @param rdt_opts a list: options for setup_F_detect
#' @param mod_severe_h a string, the S3 class name for F_incidence_h for severe malaria
#' @param severe_h_opts a list: options for setup_F_incidence_h
#' @param mod_severe_x a string, the S3 class name for F_incidence_x for severe malaria
#' @param severe_x_opts a list: options for setup_F_incidence_x
#' @param mod_moderate_h a string, the S3 class name for F_incidence_h for moderate malaria
#' @param moderate_h_opts a list: options for setup_F_incidence_h
#' @param mod_moderate_x a string, the S3 class name for F_incidence_x for moderate malaria
#' @param moderate_x_opts a list: options for setup_F_incidence_x
#' @param mod_mild_h a string, the S3 class name for F_incidence_h for mild malaria
#' @param mild_h_opts a list: options for setup_F_incidence_h
#' @param mod_mild_x a string, the S3 class name for F_incidence_x for mild malaria
#' @param mild_x_opts a list: options for setup_F_incidence_x
#' @param mod_treat a string, the S3 class name for treating malaria
#' @param treat_opts a list: options for treating malaria
#' @param mod_xi a string, the S3 class name for background drug taking
#' @param xi_opts a list: options for setup_F_incidence_x
#' @param mod_seek_care a string, the S3 class name for care seeking
#' @param care_seeking_opts a list: options to set up care seeking
#' @param mod_other_fever a string, the S3 class name for other fever
#' @param other_fever_opts a list: options to set up other fever
#' @param mod_msat a string, the S3 class name for mass screen and treat
#' @param msat_opts a list: options to set up mass screen and treat
#' @param mod_mda a string, the S3 class name for mass drug administration
#' @param mda_opts a list: options to set up mass drug administration
#' @param eta rate of loss of chemo-protection
#' @return a [list]
#' @export
make_Xpar_SIPmav3= function(nStrata, Xopts=list(),

                             # Infection
                             mod_b          = "b",     b_opts = list(b=0.55),
                             mod_c          = "c",     c_opts = list(c=0.15),
                             mod_r          = "r",     r_opts = list(r=1/180),

                             # Detection
                             mod_lm          = "d",     lm_opts  = list(d=0.7),
                             mod_pcr         = "d",     pcr_opts = list(d=0.9),
                             mod_rdt         = "d",     rdt_opts = list(d=0.7),

                             # Disease Incidence, by Severity
                             mod_severe_h    = "p",     severe_h_opts = list(p=1/500),
                             mod_severe_x    = "rate",  severe_x_opts = list(rate=1/365/5),
                             mod_moderate_h  = "p",     moderate_h_opts = list(p=1/50),
                             mod_moderate_x  = "rate",  moderate_x_opts = list(rate=1/730),
                             mod_mild_h      = "p",     mild_h_opts = list(p=1/10),
                             mod_mild_x      = "rate",  mild_x_opts = list(rate=1/50),

                             # Probability of Treating, by Severity
                             mod_treat       = "p3",    treat_opts = list(p_sev=0.95, p_mod=0.1, p_mild=0.01),

                             # Treatment
                             mod_xi          = "val",    xi_opts = list(scale=0),

                             # Care Seeking
                             mod_seek_care     = "val",    care_seeking_opts = list(val=0),
                             mod_other_fever = "val",    other_fever_opts = list(val=0),

                            # Care Seeking
                             mod_mda     = "val",      mda_opts = list(val=0),
                             mod_msat    = "val",      msat_opts = list(val=0),

                             eta=1/25){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- "SIPmav3"

    Xpar$eta      <- checkIt(eta, nStrata)

    # Dynamics
    Xpar                  <- setup_F_infect(mod_b, b_opts, Xpar, nStrata)
    Xpar                  <- setup_F_transmit(mod_c, c_opts, Xpar, nStrata)
    Xpar                  <- setup_F_clear(mod_r, r_opts, Xpar, nStrata)

    # Detection
    Xpar$lm_mod          <- setup_F_detect(mod_lm, lm_opts, nStrata)
    Xpar$pcr_mod         <- setup_F_detect(mod_pcr, pcr_opts, nStrata)
    Xpar$rdt_mod         <- setup_F_detect(mod_rdt, rdt_opts, nStrata)

    # Disease Incidence
    Xpar$severe_h            <- setup_incidence_h(mod_severe_h, severe_h_opts, nStrata)
    Xpar$severe_x            <- setup_incidence_x(mod_severe_x, severe_x_opts, nStrata)
    Xpar$moderate_h          <- setup_incidence_h(mod_moderate_h, moderate_h_opts, nStrata)
    Xpar$moderate_x          <- setup_incidence_x(mod_moderate_x, moderate_x_opts, nStrata)
    Xpar$mild_h              <- setup_incidence_h(mod_mild_h, mild_h_opts, nStrata)
    Xpar$mild_x              <- setup_incidence_x(mod_mild_x, mild_x_opts, nStrata)

    Xpar$treat3              <- setup_F_treat(mod_treat, treat_opts, nStrata)

    # Treatment
    Xpar$xi_mod              <- setup_force_malaria(mod_xi, nStrata, xi_opts)
    Xpar$seek_care_mod       <- setup_force_malaria(mod_seek_care, nStrata, care_seeking_opts)
    Xpar$other_fever_mod     <- setup_force_malaria(mod_other_fever, nStrata, other_fever_opts)

    # Mass Treatment
    Xpar$mda_mod              <- setup_force_malaria(mod_mda, nStrata, mda_opts)
    Xpar$msat_mod             <- setup_force_malaria(mod_msat, nStrata, msat_opts)

    return(Xpar)
})}


#' @title Compute Infectious Density, \eqn{X}, for `SIPmav3`
#' @description Implements [F_X] for SIPmav3 models
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SIPmav3 <- function(t, y, pars, i) {
  c <- F_transmit_y(y, pars$Xpars[[i]])
  I = y[pars$Xpar[[i]]$ix$I_ix]
  X = c*I
  return(X)
}

#' @title Size of effective infectious human population
#' @description Implements [F_X] for SIPmav3 models
#' @inheritParams ramp.xds::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_H.SIPmav3 <- function(t, y, pars, i){
  with(list_Xvars(y, pars, i),return(H))
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_prevalence] for SIPmav3 models
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @export
F_prevalence.SIPmav3 <- function(vars, Xpar) {
  x = with(vars, I/H)
  return(x)
}

#' @title Compute PR by light microscopy
#' @description Implements [F_prevalence] for SIPmav3 models
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pfpr_by_lm.SIPmav3 <- function(vars, Xpar) {
  x  = with(vars, I/H)
  d = F_detect_y(vars, Xpar, Xpar$mod_d_lm)
  return(x*d)
}

#' @title Compute PR by PCR
#' @description Implements [F_prevalence] for SIPmav3 models
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pfpr_by_pcr.SIPmav3 <- function(vars, Xpar) {
  x = with(vars, I/H)
  d = F_detect_y(vars, Xpar, Xpar$mod_d_pcr)
  return(x*d)
}

#' @title Compute PR by RDT
#' @description Implements [F_prevalence] for SIPmav3 models
#' @inheritParams ramp.xds::F_prevalence
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pfpr_by_rdt.SIPmav3 <- function(vars, Xpar) {
  pr = with(vars, I/H)
  d = F_detect_y(vars, Xpar, Xpar$mod_d_rdt)
  return(pr*d)
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for SIPmav3 models
#' @inheritParams ramp.xds::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SIPmav3 <- function(y, pars, i) {
  F_infect_y(y, pars$Xpar[[i]])
}


#' @title Return the SIPmav3 model variables as a list
#' @description This method dispatches on the type of `pars$Xpar`
#' @inheritParams ramp.xds::list_Xvars
#' @return a [list]
#' @export
list_Xvars.SIPmav3 <- function(y, pars, i) {
  with(pars$Xpar[[i]]$ix,{
       H = y[H_ix]
       I = y[I_ix]
       P = y[P_ix]
       S = H-I-P
       moi = y[moi_ix]
       aoi = y[aoi_ix]
       vh = y[vh_ix]

       return(list(S=S, I=I, P=P, moi=moi, aoi=aoi, vh=vh, H=H))})
}

#' @title Compute the HTC for the SIPmav3_xde model
#' @description Implements [HTC] for the SIPmav3_xde model with demography.
#' @inheritParams ramp.xds::HTC
#' @return a [numeric] vector
#' @export
HTC.SIPmav3 <- function(pars, i) {
  with(pars$Xpar[[i]],
       return((1-rho)*b/(r+xi)*xi/(eta+xi))
  )
}


#' @title Setup Xinits.SIPmav3
#' @description Implements [setup_Xinits] for the SIPmav3 models
#' @inheritParams ramp.xds::setup_Xinits
#' @return a [list] vector
#' @export
setup_Xinits.SIPmav3 = function(pars, H, i, Xopts=list()){
  pars$Xinits[[i]] = with(pars,make_Xinits_SIPmav3(pars$nStrata[i], H, Xopts))
  return(pars)
}

#' @title Make initial values for a SIPmav3 human model, with defaults
#' @param nStrata the number of population strata
#' @param Xopts a [list] that could overwrite defaults
#' @param H the initial human population density
#' @param I the initial values of the variable I
#' @param P the initial values of the variable P
#' @param moi the initial values of MoI
#' @param aoi initial values of AoI
#' @param vh initial values for cumulative exposure
#' @param severe initial values for cumulative severe malaria
#' @param moderate initial values for cumulative moderate malaria
#' @param mild initial values for cumulative moderate malaria
#' @param treat initial values for cumulative treatment
#' @return a [list]
#' @export
make_Xinits_SIPmav3 = function(nStrata, H, Xopts = list(),
                           I=1, P=0, moi=0, aoi=0, vh=0,
                           severe=0, moderate=0, mild=0,
                           treat=0){with(Xopts, {
    I = unname(as.vector(checkIt(I, nStrata)))
    P = unname(as.vector(checkIt(P, nStrata)))
    moi = unname(as.vector(checkIt(moi, nStrata)))
    aoi = unname(as.vector(checkIt(aoi, nStrata)))
    vh = unname(as.vector(checkIt(vh, nStrata)))
    severe = unname(as.vector(checkIt(severe, nStrata)))
    moderate = unname(as.vector(checkIt(moderate, nStrata)))
    mild = unname(as.vector(checkIt(mild, nStrata)))
    treat = unname(as.vector(checkIt(treat, nStrata)))
    return(list(H=H, I=I, P=P, moi=moi, aoi=aoi, vh=vh, severe=severe, moderate=moderate, mild=mild, treat=treat))
})}

#' @title Parse the output of deSolve and return variables for SIPmav3 models
#' @description Implements [parse_Xorbits] for SIPmav3 models
#' @inheritParams ramp.xds::parse_Xorbits
#' @return none
#' @export
parse_Xorbits.SIPmav3 <- function(outputs, pars, i) {with(pars$Xpar[[i]]$ix,{
    H <- outputs[,H_ix]
    I <- outputs[,I_ix]
    P <- outputs[,P_ix]
    moi <- outputs[,moi_ix]
    aoi <- outputs[,aoi_ix]
    vh <- outputs[,vh_ix]
    severe <- outputs[,severe_ix]
    moderate <- outputs[,moderate_ix]
    mild <- outputs[,mild_ix]
    treat <- outputs[,treat_ix]
    S <- H-I-P
    ni <- F_transmit_vars(list(vh=vh,moi=moi,aoi=aoi), pars$Xpar[[i]])*I/H
    prevalence <- I/H
    return(list(time=time, S=S, I=I, P=P, H=H,
                moi=moi, aoi=aoi, vh=vh,
                severe=severe, moderate=moderate, mild=mild,
                treat=treat, prevalence=prevalence))
})}

#' @title Add indices for human population to parameter list
#' @description Implements [setup_Xix] for SIPmav3 models
#' @inheritParams ramp.xds::setup_Xix
#' @return none
#' @importFrom utils tail
#' @export
setup_Xix.SIPmav3 <- function(pars, i) {with(pars,{

  H_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(H_ix, 1)

  I_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(I_ix, 1)

  P_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(P_ix, 1)

  moi_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(moi_ix, 1)

  aoi_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(aoi_ix, 1)

  vh_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(vh_ix, 1)

  severe_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(severe_ix, 1)

  moderate_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(moderate_ix, 1)

  mild_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(mild_ix, 1)

  treat_ix <- seq(from = max_ix+1, length.out=nStrata[i])
  max_ix <- tail(treat_ix, 1)

  pars$max_ix = max_ix
#  pars$ix$X[[i]]
  pars$Xpar[[i]]$ix = list(H_ix=H_ix, I_ix=I_ix, P_ix=P_ix, moi_ix=moi_ix,
                           aoi_ix=aoi_ix, vh_ix=vh_ix,
                           severe_ix=severe_ix, moderate_ix=moderate_ix, mild_ix=mild_ix,
                           treat_ix=treat_ix)
  return(pars)
})}

#' @title Update inits for SIPmav3 models from a vector of states
#' @inheritParams ramp.xds::update_Xinits
#' @return none
#' @export
update_Xinits.SIPmav3 <- function(pars, y, i) {
  with(list_Xvars(y, pars, i),{
    pars$Xinits[[i]] = make_Xinits_SIPmav3(pars$nStrata[i], H, list(), I=I, P=P, moi=moi, aoi=aoi, vh=vh,
                                           severe=severe, moderate=moderate, mild=mild, treat=treat)
    return(pars)
  })}


#' @title Return initial values as a vector from a SIPmav3 model
#' @description This method dispatches on the type of `pars$Xpar[[i]]`
#' @inheritParams ramp.xds::get_Xinits
#' @return none
#' @export
get_Xinits.SIPmav3 <- function(pars, i){pars$Xinits[[i]]}


#' Plot the density of infected individuals for the SIPmav3 model
#'
#' @inheritParams ramp.xds::xds_plot_X
#' @export
xds_plot_X.SIPmav3 = function(pars, i=1, clrs=c("darkblue", "darkred", "darkgreen"), llty=1, add=FALSE){
  XH = pars$outputs$orbits$XH[[i]]
  times = pars$outputs$time

  if(add==FALSE)
    plot(times, 0*times, type = "n", ylim = c(0, max(XH$H)),
         ylab = "# Infected", xlab = "Time")

  xds_lines_X_SIPmav3(times, XH, pars$nStrata[i], clrs, llty)
}


#' Add lines for the density of infected individuals for the SIP model
#'
#' @param times time points for the observations
#' @param XH a list with the outputs of parse_Xorbits.SIP
#' @param nStrata the number of population strata
#' @param clrs a vector of colors
#' @param llty an integer (or integers) to set the `lty` for plotting
#'
#' @export
xds_lines_X_SIPmav3 = function(times, XH, nStrata, clrs=c("darkblue", "darkred", "darkgreen"), llty=1){
  if (length(llty)< nStrata) llty = rep(llty, nStrata)
  with(XH,{
    if(nStrata == 1){
      lines(times, S, col=clrs[1], lty = llty)
      lines(times, I, col=clrs[2], lty = llty)
      lines(times, P, col=clrs[3], lty = llty)
    } else {
      for(i in 1:nStrata)
        lines(times, S[,i], col=clrs[1], lty = llty[i])
        lines(times, I[,i], col=clrs[2], lty = llty[i])
        lines(times, P[,i], col=clrs[3], lty = llty[i])
    }})
}


#' @title Compute the steady states for the SIP model as a function of the daily foi
#' @description Compute the steady state of the SIP model as a function of the daily eir.
#' @inheritParams ramp.xds::xde_steady_state_X
#' @return the steady states as a named vector
#' @export
xde_steady_state_X.SIPmav3 = function(foi, H, Xpar){with(Xpar,{
  Ieq = (foi*H*eta*(1-rho))/((foi+r+xi)*(eta+xi) +foi*(r-eta)*rho)
  Peq  = (H*xi*(foi+r+xi) + (foi*H*r*rho))/((foi+r+xi)*(eta+xi) +foi*(r-eta)*rho)
  Seq = H -Ieq - Peq
  return(list(S=as.vector(Seq), I=as.vector(Ieq), P = as.vector(Peq)))
})}

