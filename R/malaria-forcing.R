
#' @title Mass Treatment
#' @description Set up a function to force some aspect of malaria
#' @param t current simulation time
#' @param y state vector
#' @param Xpars an `Xmod` object
#' @param forcing_mod a function
#' @return a [numeric] vector of length `nStrata`
#' @export
F_force_malaria <- function(t, y, Xpars, forcing_mod) {
  UseMethod("F_force_malaria", forcing_mod)
}

#' @title Make a function to simulate mass treatment
#' @description Mass treatment
#' @param mod_forcing is a string that defines an S3 class
#' @param nStrata  the number of population strata
#' @param opts a [list] that overwrites default values
#' @return none
#' @export
setup_force_malaria = function(mod_forcing, nStrata, opts=list()){
  class(mod_forcing) <- mod_forcing
  UseMethod("setup_force_malaria", mod_forcing)
}

#' @title Mass Treatment
#' @description Implements F_forcing for the trivial model.
#' @inheritParams F_force_malaria
#' @return a [numeric] vector of length ``
#' @export
F_force_malaria.val <- function(t, y, Xpars, forcing_mod) {
  return(forcing_mod$val)
}

#' @title Make a function to simulate mass treatment
#' @description Mass treatment
#' @inheritParams setup_force_malaria
#' @return none
#' @export
setup_force_malaria.val = function(mod_forcing, nStrata, opts=list()){
  with(opts, setup_force_malaria_val(nStrata, opts))
}

#' @title Make a function to simulate mass treatment
#' @description Mass treatment
#' @param nStrata the number of population strata
#' @param opts a [list] that overwrites default values
#' @param value the value to return
#' @return none
#' @export
setup_force_malaria_val = function(nStrata, opts=list(), value=0){
  with(opts,{
    forcing_par = list()
    class(forcing_par) <- "val"
    forcing_par$val = checkIt(value, nStrata)
    return(forcing_par)
})}


#' @title Mass Treatment
#' @description Implements forcing for the trivial model.
#' @inheritParams F_force_malaria
#' @return a [numeric] vector of length ``
#' @export
F_force_malaria.func <- function(t, y, Xpars, forcing_mod) {
  rate = with(forcing_mod, scale*F_season(t)*F_trend(t))
  return(rate)
}

#' @title Make a function to simulate mass treatment
#' @description Mass treatment
#' @inheritParams setup_force_malaria
#' @return a forcing function
#' @export
setup_force_malaria.func = function(mod_forcing, nStrata, opts=list()){
  with(opts, setup_force_malaria_func(nStrata, opts))
}

#' @title Make a function to simulate mass treatment
#' @description Mass treatment
#' @param nStrata the number of population strata
#' @param opts a [list] that overwrites default values
#' @param scale a scaling parameter
#' @param F_season a function describing a seasonal pattern over time
#' @param season_par an object to configure a seasonality function using [make_function]
#' @param F_trend a function describing a temporal trend over time
#' @param trend_par an object to configure a trends function using [make_function]
#' @return none
#' @export
setup_force_malaria_func  = function(nStrata, opts=list(),
                                      scale=1,
                                      F_season=F_flat, season_par = list(),
                                      F_trend=F_flat, trend_par = list()){
  with(opts,{
    forcing_par = list()
    class(forcing_par) <- "func"
    forcing_par$scale = checkIt(scale, nStrata)

    forcing_par$F_season = F_season
    forcing_par$season_par <- season_par
    if(length(season_par)>0)
      forcing_par$F_season <- make_function(season_par)

    forcing_par$F_trend = F_trend
    forcing_par$trend_par <- trend_par
    if(length(trend_par)>0)
      forcing_par$F_trend <- make_function(trend_par)

    return(forcing_par)
  })}

