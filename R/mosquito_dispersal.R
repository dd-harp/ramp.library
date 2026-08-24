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

#' @title Change Mosquito Dispersal Matrix
#' @description
#' Run [check_K_matrix] then
#'
#' After passing checks, `xds_obj` is updated.
#'
#' In models with multiple species, use `s` to
#' specify the species to update.
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
  return(xds_obj)
}

#' @title Change Mosquito Dispersal Matrix
#' @description
#' Run [check_K_matrix] then
#'
#' After passing checks, `xds_obj` is updated.
#'
#' In models with multiple species, use `s` to
#' specify the species to update.
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
  return(xds_obj)
}
