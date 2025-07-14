
#' @title Mass Treatment
#' @description This method dispatches on the type of `pars$mass_treat`.
#' @param t current simulation time
#' @param y state vector
#' @param Xpars an `Xmod` object
#' @param take_drugs_mod a function
#' @return a [numeric] vector of length `nStrata`
#' @export
F_take_drugs <- function(t, y, Xpars, take_drugs_mod) {
  UseMethod("F_take_drugs", take_drugs_mod)
}

#' @title Mass Treatment
#' @description Implements [F_take_drugs] for the trivial model.
#' @inheritParams F_take_drugs
#' @return a [numeric] vector of length ``
#' @export
F_take_drugs.func <- function(t, y, Xpars, take_drugs_mod) {
  rate = with(take_drugs_mod, scale*F_season(t)*F_trend(t))
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
setup_drug_taking_function = function(nStrata, opts=list(),
                             scale=1,
                             F_season=F_flat, season_par = list(),
                             F_trend=F_flat, trend_par = list()){
  with(opts,{
    take_drugs_par = list()
    class(take_drugs_par) <- "func"
    take_drugs_par$scale = checkIt(scale, nStrata)

    take_drugs_par$F_season = F_season
    take_drugs_par$season_par <- season_par
    if(length(season_par)>0)
      take_drugs_par$F_season <- make_function(season_par)

    take_drugs_par$F_trend = F_trend
    take_drugs_par$trend_par <- trend_par
    if(length(trend_par)>0)
      take_drugs_par$F_trend <- make_function(trend_par)

    return(take_drugs_par)
 })}

