#' @title Get the Mosquito Dispersal Matrix
#'
#' @inheritParams ramp.xds::get_K_matrix
#'
#' @return an **`xds`** object
#' @keywords internal
#'
#' @export
get_K_matrix.b = function(xds_obj, behave="b", s=1){
  return(xds_obj$MY_obj[[s]]$Kb_matrix)
}

#' @title Get the Mosquito Dispersal Matrix
#'
#' @inheritParams ramp.xds::get_K_matrix
#'
#' @return an **`xds`** object
#' @keywords internal
#'
#' @export
get_K_matrix.q = function(xds_obj, behave="q", s=1){
  return(xds_obj$MY_obj[[s]]$Kq_matrix)
}

#' @title Change Mosquito Dispersal Matrix
#' @description Change the dispersal matrix for
#' blood feeding
#'
#' @inheritParams ramp.xds::change_K_matrix
#'
#' @return an **`xds`** object
#' @keywords internal
#' @seealso [xds_info_mosquito_dispersal]; [setup_K_matrix]
#' @export
change_K_matrix.b = function(K_matrix, xds_obj, behave="b", s=1){
  stopifnot(is.matrix(K_matrix))
  stopifnot(dim(K_matrix) == rep(xds_obj$nPatches,2))
  stopifnot(colSums(K_matrix) < 1e-7)
  xds_obj$MY_obj[[s]]$Kb_matrix <- K_matrix
  xds_obj$MY_obj[[s]]$Omega_b <-
    with(xds_obj$MY_obj[[s]], compute_Omega_xde(g, sigma_b, mu, Kb_matrix))
  return(xds_obj)
}

#' @title Change Mosquito Dispersal Matrix
#' @description Change the dispersal matrix for
#' egg laying
#'
#' @inheritParams ramp.xds::change_K_matrix
#'
#' @return an **`xds`** object
#' @keywords internal
#' @seealso [xds_info_mosquito_dispersal]; [setup_K_matrix]
#' @export
change_K_matrix.q = function(K_matrix, xds_obj, behave="q", s=1){
  stopifnot(is.matrix(K_matrix))
  stopifnot(dim(K_matrix) == rep(xds_obj$nPatches,2))
  stopifnot(colSums(K_matrix) < 1e-7)
  xds_obj$MY_obj[[s]]$Kq_matrix <- K_matrix
  xds_obj$MY_obj[[s]]$Omega_q <-
    with(xds_obj$MY_obj[[s]], compute_Omega_xde(g, sigma_q, mu, Kq_matrix))
  return(xds_obj)
}
