
#' @title Setup a Patch Emigration Bionomic Object
#'
#' @description Set up an object
#' to compute the human fraction, \eqn{sigma}
#'
#' @param sigma_b the mosquito patch emigration rate
#' @param MY_obj an **`MY`** model object
#'
#' @return a **`MY`** model object
#'
#' @keywords internal
#' @export
setup_sigma_b_obj = function(sigma_b, MY_obj){
  MY_obj$sigma_b = sigma_b
  MY_obj$sigma_b_t = sigma_b
  MY_obj$es_sigma_b = 1
  MY_obj$sigma_b_obj <- list()
  class(MY_obj$sigma_b_obj) <- "static"
  MY_obj$sigma_b_obj$sigma_b <- sigma_b
  return(MY_obj)
}

#' @title Setup a Patch Emigration Bionomic Object
#'
#' @description Set up an object
#' to compute the human fraction, \eqn{sigma}
#'
#' @param sigma_q the mosquito patch emigration rate
#' @param MY_obj an **`MY`** model object
#'
#' @return a **`MY`** model object
#'
#' @keywords internal
#' @export
setup_sigma_q_obj = function(sigma_q, MY_obj){
  MY_obj$sigma_q = sigma_q
  MY_obj$sigma_q_t = sigma_q
  MY_obj$es_sigma_q = 1
  MY_obj$sigma_q_obj <- list()
  class(MY_obj$sigma_q_obj) <- "static"
  MY_obj$sigma_q_obj$sigma_q <- sigma_q
  return(MY_obj)
}

#' @title Compute the Mosquito Patch Emigration Rate
#'
#' @description This method dispatches on the type of `sigma_obj`. It should
#' set the values the patch emigration rate, \eqn{\sigma}
#'
#' @inheritParams ramp.xds::F_sigma
#'
#' @return a [numeric] vector osigma length `nPatches`
#'
#' @keywords internal
#' @export
F_sigma_b = function(t, xds_obj, s){
  UseMethod("F_sigma_b", xds_obj$MY_obj[[s]]$sigma_b_obj)
}

#' @title Compute the Mosquito Patch Emigration Rate
#'
#' @description This method dispatches on the type of `sigma_obj`. It should
#' set the values the patch emigration rate, \eqn{\sigma}
#'
#' @inheritParams ramp.xds::F_sigma
#'
#' @return a [numeric] vector osigma length `nPatches`
#'
#' @keywords internal
#' @export
F_sigma_q = function(t, xds_obj, s){
  UseMethod("F_sigma_q", xds_obj$MY_obj[[s]]$sigma_q_obj)
}

#' @title Static model patch emigration
#'
#' @description Implements [F_sigma_b] for a static model
#'
#' @inheritParams ramp.xds::F_sigma
#'
#' @return \eqn{sigma}, the patch emigration rate
#' @keywords internal
#' @export
F_sigma_b.static = function(t, xds_obj, s){
  return(xds_obj$MY_obj[[s]]$sigma_b_obj$sigma_b)
}

#' @title Static model patch emigration
#'
#' @description Implements [F_sigma_q] for a static model
#'
#' @inheritParams ramp.xds::F_sigma
#'
#' @return \eqn{sigma}, the patch emigration rate
#' @keywords internal
#' @export
F_sigma_q.static = function(t, xds_obj, s){
  return(xds_obj$MY_obj[[s]]$sigma_q_obj$sigma_q)
}
