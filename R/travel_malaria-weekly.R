# Specialized methods to compute attack rates (ar) for travel_malaria models

#' @title A model for the travel Attack Rate
#' @description Implements [travel_malaria] through a model for the travel Attack Rate
#' @inheritParams ramp.xde::travel_malaria
#' @return a [numeric]
#' @export
travel_malaria.multiday <- function(t, pars) {
  with(pars$TRAVEL,{
    return(delta_scale*delta_t(t, pars))
  })}

#' @title Set up parameters and functions to compute the Attack-Rate for a travel_malaria model
#' @param pars a [list]
#' @param travel_Opts a list, overwrites default values
#' @param delta_scale a non-negative numeric value to scale the mean travel related malaria Attack Rate
#' @param delta_t the temporal pattern for travel related malaria Attack-Rate
#' @return none
#' @export
setup_travel_multiday <- function(pars, travel_Opts = list(), delta_scale=0, delta_t=NULL) {with(travel_Opts,{

  TRAVEL <- list()
  class(TRAVEL) <- 'multiday'
  pars$TRAVEL <- TRAVEL

  pars$TRAVEL$delta_scale = delta_scale
  if(is.null(delta_t)) delta_t = function(t, pars){1}
  pars$TRAVEL$delta_t = delta_t

  return(pars)
})}
