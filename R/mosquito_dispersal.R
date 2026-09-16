#' @title Get the Mosquito Dispersal Matrix
#'
#' @inheritParams ramp.xds::get_K_matrix
#'
#' @return an **`xds`** object
#' @keywords internal
#'
#' @export
get_K_matrix.Kb = function(xds_obj, which_K="Kb", s=1){
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
get_K_matrix.Kq = function(xds_obj, which_K="Kq", s=1){
  return(xds_obj$MY_obj[[s]]$Kq_matrix)
}

#' @title Change mosquito dispersal matrix
#'
#' @description
#' Change the dispersal
#' matrix for blood feeding, called the `Kb_matrix`
#'
#' @inheritParams ramp.xds::change_K_matrix
#'
#' @return an **`xds`** object
#' @keywords internal
#' @seealso [xds_info_mosquito_dispersal]; [setup_K_matrix]
#' @export
change_K_matrix.Kb = function(K_matrix, xds_obj, which_K="Kb", s=1){
  check_K_matrix(K_matrix, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$Kb_matrix <- K_matrix
  xds_obj$MY_obj[[s]]$Kb_obj <- make_static_obj()
  xds_obj$MY_obj[[s]]$Omega_b_obj <- make_static_obj()
  xds_obj$MY_obj[[s]]$Omega_b <- F_Omega_xde(xds_obj, s)
  return(xds_obj)
}

#' @title Change mosquito dispersal matrix
#'
#' @description
#' Change the dispersal
#' matrix for egg laying, called the `Kq_matrix`
#'
#' @inheritParams ramp.xds::change_K_matrix
#'
#' @return an **`xds`** object
#' @keywords internal
#' @seealso [xds_info_mosquito_dispersal]; [setup_K_matrix]
#' @export
change_K_matrix.Kq = function(K_matrix, xds_obj, which_K="Kq", s=1){
  check_K_matrix(K_matrix, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$Kq_matrix <- K_matrix
  xds_obj$MY_obj[[s]]$Kq_obj <- make_static_obj()
  xds_obj$MY_obj[[s]]$Omega_q_obj <- make_static_obj()
  return(xds_obj)
}

#' @title Change mosquito dispersal matrix `Ks`
#' @description
#' Change the dispersal
#' matrix for sugar feeding, called the `Ks_matrix`
#'
#' @inheritParams ramp.xds::change_K_matrix
#'
#' @return an **`xds`** object
#' @keywords internal
#' @seealso [xds_info_mosquito_dispersal]; [setup_K_matrix]
#' @export
change_K_matrix.Ks = function(K_matrix, xds_obj, which_K="Ks", s=1){
  check_K_matrix(K_matrix, xds_obj$nPatches)
  xds_obj$MY_obj[[s]]$Ks_matrix <- K_matrix
  xds_obj$MY_obj[[s]]$Ks_obj <- make_static_obj()
  xds_obj$MY_obj[[s]]$Omega_s_obj <- make_static_obj()
  return(xds_obj)
}

