
#' @title Exposure and Infection
#' @description This function translates seven days of daily entomological
#' inoculation rate (dEIR) into a multiday attack rate. The daily FoI
#' is the sum of two terms: 1) a function [F_foi] computes the local dFoI;
#' 2) a function [travel_malaria] that computes the FoI resulting from
#' exposure while traveling.
#' @inheritParams ramp.xde::Exposure
#' @return the function modifies **pars** and returns it: the computed FoI are stored as `pars$FoI`
#' @export
Exposure.multiday <- function(t, y, pars){
  dd = t%%pars$D+1
  for(i in 1:pars$nHosts){
    b = F_b(y, pars, i)
    pars$EIRD[[i]] = F_ar(pars$EIR[[i]], b, pars) + travel_malaria(t, pars)
  }
  return(pars)
}
