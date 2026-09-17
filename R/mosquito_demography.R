

#' @title Change Omega
#'
#' @description Update the demographic matrix for
#' behavioral state models
#'
#' @inheritParams ramp.xds::change_Omega
#'
#' @return an **xds** model object
#' @keywords internal
#' @export
change_Omega_b = function(xds_obj, s){
  with(xds_obj$MY_obj[[s]],{
    xds_obj$MY_obj[[s]]$Omega_b <- F_Omega_xde(g, sigma_b, mu, Kb_matrix)
    return(xds_obj)
  })}

#' @title Change Omega
#'
#' @description Update the demographic matrix for
#' behavioral state models
#'
#' @inheritParams ramp.xds::change_Omega
#'
#' @return an **xds** model object
#' @keywords internal
#' @export
change_Omega_q = function(xds_obj, s){
  with(xds_obj$MY_obj[[s]],{
    xds_obj$MY_obj[[s]]$Omega_q <- F_Omega_xde(g, sigma_q, mu, Kq_matrix)
    return(xds_obj)
  })}

#' @title Change Omega
#'
#' @description Update the demographic matrix for
#' behavioral state models
#'
#' @inheritParams ramp.xds::change_Omega
#'
#' @return an **xds** model object
#' @keywords internal
#' @export
change_Omega_s = function(xds_obj, s){
  with(xds_obj$MY_obj[[s]],{
    xds_obj$MY_obj[[s]]$Omega_s <- F_Omega_xde(g, sigma_s, mu, Ks_matrix)
    return(xds_obj)
  })}
