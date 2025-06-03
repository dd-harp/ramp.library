
#' @title Mass Treatment
#' @description This method dispatches on the type of `pars$mass_treat`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param i the species index
#' @return a [numeric] vector of length `nStrata`
#' @export
F_mass_treat <- function(t, y, pars, i) {
  UseMethod("F_mass_treat", pars$Xpar[[i]]$mass_treat_par)
}



#' @title Mass Treatment
#' @description Implements [F_mass_treat] for the trivial model.
#' @inheritParams F_mass_treat
#' @return a [numeric] vector of length ``
#' @export
F_mass_treat.func <- function(t, y, pars, i) {
  rate = with(pars$Xpar[[i]]$mass_treat_par, scale*F_season(t)*F_trend(t))
  return(rate)
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
make_mass_treat = function(nStrata, opts=list(),
                             scale=1,
                             F_season=F_flat, season_par = list(),
                             F_trend=F_flat, trend_par = list()){
  with(Lopts,{
    mass_treat_par = list()
    class(mass_treat_par) <- "func"
    mass_treat_par$scale = checkIt(scale, nStrata)

    mass_treat_par$F_season = F_season
    mass_treat_par$season_par <- season_par
    if(length(season_par)>0)
      mass_treat_par$F_season <- make_function(season_par)

    mass_treat_par$F_trend = F_trend
    mass_treat_par$trend_par <- trend_par
    if(length(trend_par)>0)
      mass_treat_par$F_trend <- make_function(trend_par)

    return(mass_treat_par)
  })}

#' @title Mass Treatment
#' @description This method dispatches on the type of `pars$mass_test_treat`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param i the species index
#' @return a [numeric] vector of length `nStrata`
#' @export
F_mass_test_treat <- function(t, y, pars, i) {
  UseMethod("F_mass_test_treat", pars$Xpar[[i]]$mass_test_treat_par)
}

#' @title Mass Treatment
#' @description Implements [F_mass_test_treat] for the trivial model.
#' @inheritParams F_mass_test_treat
#' @return a [numeric] vector of length ``
#' @export
F_mass_test_treat.func <- function(t, y, pars, i) {
  rate = with(pars$Xpar[[i]]$mass_test_treat_par, scale*F_season(t)*F_trend(t))
  return(rate)
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
make_mass_test_treat = function(nStrata, opts=list(),
                             scale=1,
                             F_season=F_flat, season_par = list(),
                             F_trend=F_flat, trend_par = list()){
  with(Lopts,{
    mass_test_treat_par = list()
    class(mass_test_treat_par) <- "func"
    mass_test_treat_par$scale = checkIt(scale, nStrata)

    mass_test_treat_par$F_season = F_season
    mass_test_treat_par$season_par <- season_par
    if(length(season_par)>0)
      mass_test_treat_par$F_season <- make_function(season_par)

    mass_test_treat_par$F_trend = F_trend
    mass_test_treat_par$trend_par <- trend_par
    if(length(trend_par)>0)
      mass_test_treat_par$F_trend <- make_function(trend_par)

    return(mass_test_treat_par)
})}
