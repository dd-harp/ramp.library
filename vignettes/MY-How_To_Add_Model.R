## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Make the model object for newMYname
#' @param nPatches number of patches
#' @param options named values that override defaults
#' @return a list
#' @export
make_MY_obj_newMYname <- function(nPatches, options = list(), p1 = 1) {
  with(options, {
    MY_obj <- list()
    MY_obj$p1 <- checkIt(p1, nPatches)
    MY_obj$nPatches <- nPatches
    class(MY_obj) <- "newMYname"
    MY_obj
  })
}

#' @title Set up the newMYname MY module
#' @description Implements [ramp.xds::setup_MY_obj] for newMYname.
#' @inheritParams ramp.xds::setup_MY_obj
#' @return an **`xds`** model object
#' @export
setup_MY_obj.newMYname <- function(MYname, xds_obj, options = list(), s = 1) {
  xds_obj$MY_obj[[s]] <- make_MY_obj_newMYname(xds_obj$nPatches, options)
  xds_obj <- setup_MY_inits(xds_obj, options = options, s = s)
  xds_obj
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Set up indices for newMYname
#' @description Implements [ramp.xds::setup_MY_ix] for newMYname.
#' @inheritParams ramp.xds::setup_MY_ix
#' @return an **`xds`** model object
#' @importFrom utils tail
#' @export
setup_MY_ix.newMYname <- function(xds_obj, s) {
  with(xds_obj, {
    M_ix <- seq(from = max_ix + 1, length.out = nPatches)
    max_ix <- tail(M_ix, 1)

    Y_ix <- seq(from = max_ix + 1, length.out = nPatches)
    max_ix <- tail(Y_ix, 1)

    xds_obj$max_ix <- max_ix
    xds_obj$MY_obj[[s]]$ix <- list(M_ix = M_ix, Y_ix = Y_ix)
    xds_obj
  })
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Get newMYname state variables
#' @description Implements [ramp.xds::get_MY_vars] for newMYname.
#' @inheritParams ramp.xds::get_MY_vars
#' @return a named list
#' @export
get_MY_vars.newMYname <- function(y, xds_obj, s) {
  with(xds_obj$MY_obj[[s]]$ix, {
    list(
      M = y[M_ix],
      Y = y[Y_ix]
    )
  })
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Parse newMYname model output
#' @description Implements [ramp.xds::parse_MY_orbits] for newMYname.
#' @inheritParams ramp.xds::parse_MY_orbits
#' @return a named list
#' @export
parse_MY_orbits.newMYname <- function(outputs, xds_obj, s) {
  with(xds_obj$MY_obj[[s]]$ix, {
    M <- outputs[, M_ix]
    Y <- outputs[, Y_ix]
    list(M = M, Y = Y)
  })
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Make newMYname initial values
#' @param nPatches number of patches
#' @param options named values that override defaults
#' @param M initial adult mosquito density
#' @param Y initial infected mosquito density
#' @return a named list
#' @export
make_MY_inits_newMYname <- function(nPatches, options = list(), M = 10, Y = 1) {
  with(options, {
    M <- checkIt(M, nPatches)
    Y <- checkIt(Y, nPatches)
    list(M = M, Y = Y)
  })
}

#' @title Set up newMYname initial values
#' @description Implements [ramp.xds::setup_MY_inits] for newMYname.
#' @inheritParams ramp.xds::setup_MY_inits
#' @return an **`xds`** model object
#' @export
setup_MY_inits.newMYname <- function(xds_obj, options = list(), s = 1) {
  xds_obj$MY_obj[[s]]$inits <-
    make_MY_inits_newMYname(xds_obj$nPatches, options)
  xds_obj
}

#' @title Change newMYname initial values
#' @description Implements [ramp.xds::change_MY_inits] for newMYname.
#' @inheritParams ramp.xds::change_MY_inits
#' @return an **`xds`** model object
#' @export
change_MY_inits.newMYname <- function(xds_obj, options = list(), s = 1) {
  xds_obj$MY_obj[[s]]$inits[names(options)] <- options
  xds_obj
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Compute newMYname derivatives
#' @description Implements [ramp.xds::dMYdt] for newMYname.
#' @inheritParams ramp.xds::dMYdt
#' @return a numeric vector of derivatives
#' @export
dMYdt.newMYname <- function(t, y, xds_obj, s) {
  # Retrieve the shared terms needed by this model.
  Lambda <- xds_obj$terms$Lambda[[s]]
  kappa <- xds_obj$terms$kappa[[s]]

  with(get_MY_vars(y, xds_obj, s), {
    with(xds_obj$MY_obj[[s]], {
      # Replace these placeholders with the model's equations.
      dM <- ...
      dY <- ...
      c(dM, dY)
    })
  })
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Update newMYname states
#' @description Implements [ramp.xds::Update_MYt] for newMYname.
#' @inheritParams ramp.xds::Update_MYt
#' @return a numeric vector of updated states
#' @export
Update_MYt.newMYname <- function(t, y, xds_obj, s) {
  with(get_MY_vars(y, xds_obj, s), {
    with(xds_obj$MY_obj[[s]], {
      M_next <- ...
      Y_next <- ...
      c(M_next, Y_next)
    })
  })
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Compute newMYname infective biting density
#' @description Implements [ramp.xds::F_fqZ] for newMYname.
#' @inheritParams ramp.xds::F_fqZ
#' @return a numeric vector with one value per patch
#' @export
F_fqZ.newMYname <- function(t, y, xds_obj, s) {
  with(get_MY_vars(y, xds_obj, s), {
    with(xds_obj$MY_obj[[s]], {
      fqZ <- ...
      fqZ
    })
  })
}

#' @title Compute newMYname human biting density
#' @description Implements [ramp.xds::F_fqM] for newMYname.
#' @inheritParams ramp.xds::F_fqM
#' @return a numeric vector with one value per patch
#' @export
F_fqM.newMYname <- function(t, y, xds_obj, s) {
  with(get_MY_vars(y, xds_obj, s), {
    with(xds_obj$MY_obj[[s]], {
      fqM <- ...
      fqM
    })
  })
}

#' @title Compute newMYname egg production
#' @description Implements [ramp.xds::F_eggs] for newMYname.
#' @inheritParams ramp.xds::F_eggs
#' @return a numeric vector with one value per patch
#' @export
F_eggs.newMYname <- function(t, y, xds_obj, s) {
  with(get_MY_vars(y, xds_obj, s), {
    with(xds_obj$MY_obj[[s]], {
      eggs <- ...
      eggs
    })
  })
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Set up newMYname skill set
#' @inheritParams ramp.xds::setup_skillset_MY
#' @return an **`xds`** model object
#' @export
setup_skillset_MY.newMYname <- function(xds_obj, s) {
  xds_obj$MY_obj[[s]]$skill_set <- list()
  xds_obj
}

#' @title Check newMYname setup
#' @inheritParams ramp.xds::check_MY
#' @return an **`xds`** model object
#' @export
check_MY.newMYname <- function(xds_obj, s) {
  xds_obj
}

