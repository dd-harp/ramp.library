## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title The `newMYname` module for the **MY** component
#' @section State Variables:
#' Describe the module's actual states and any derived quantities used by
#' its initial-density accessors.
#' @section Parameters:
#' Describe the parameters used by the implemented model, including their
#' biological meanings and units. Include `f`, `q`, `g`, and `sigma` if the
#' corresponding dynamical parameter accessors are implemented.
#' @section Dynamical Terms:
#' Describe the shared terms used by the implemented equations.
#' @section ODE:
#' Describe the continuous-time equations if the module implements `dMYdt`.
#' @section DTS:
#' Describe the updating equations if the module implements `Update_MYt`.
#' @section Egg Laying:
#' Describe the daily number of eggs laid per patch as the module's egg-laying
#' population times oviposition rate and eggs per batch.
#' [F_eggs.newMYname] returns this quantity before distribution among aquatic
#' habitats by [ramp.xds::EggLaying].
#' @name newMYname
NULL

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Compute the derivatives for `newMYname` (**MY**)
#' @description Compute the derivatives for the `newMYname` module:
#' Replace this placeholder with the equations for `dM` and `dY`,
#' in returned-vector order, when replacing the placeholder code below.
#' @inheritParams ramp.xds::dMYdt
#' @note Derivatives are returned in blocks for `M`, then `Y`, each in patch order.
#' @return a [numeric] vector of derivatives in the module's state-variable order
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
#' @title Update states for `newMYname` (**MY**)
#' @description Update the state variables for the `newMYname` module:
#' \deqn{\begin{array}{rl}
#' M_t &= \ldots \\\\
#' Y_t &= \ldots
#' \end{array}}
#' Replace these placeholders with the implemented updating equations.
#' @inheritSection newMYname State Variables
#' @inheritSection newMYname Dynamical Terms
#' @inheritSection newMYname Parameters
#' @inheritSection newMYname DTS
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
#' @title Get the state variables for `newMYname` (**MY**)
#' @description Extract the state variables for the `newMYname` module
#' from `y` and return them by name in a list for use during model solving.
#' @note Variable indices are stored at `xds_obj$MY_obj[[s]]$ix`.
#' @inheritParams ramp.xds::get_MY_vars
#' @return a named [list] of state variables and module-specific derived quantities
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
#' @title Compute infective mosquito biting density for `newMYname` (**MY**)
#' @description Return the daily infective mosquito biting density on humans,
#' `fqZ`, for each patch in the `newMYname` module.
#' @inheritParams ramp.xds::F_fqZ
#' @seealso [F_fqM.newMYname]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ.newMYname <- function(t, y, xds_obj, s) {
  with(get_MY_vars(y, xds_obj, s), {
    with(xds_obj$MY_obj[[s]], {
      fqZ <- ...
      fqZ
    })
  })
}

#' @title Compute mosquito biting density for `newMYname` (**MY**)
#' @description Return the daily mosquito biting density on humans,
#' `fqM`, for each patch in the `newMYname` module.
#' @inheritParams ramp.xds::F_fqM
#' @seealso [F_fqZ.newMYname]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqM.newMYname <- function(t, y, xds_obj, s) {
  with(get_MY_vars(y, xds_obj, s), {
    with(xds_obj$MY_obj[[s]], {
      fqM <- ...
      fqM
    })
  })
}

#' @title Compute egg production for `newMYname` (**MY**)
#' @description Return the number of eggs laid by the adult mosquito population, per patch
#' @inheritSection newMYname Egg Laying
#' @inheritParams ramp.xds::F_eggs
#' @seealso [newMYname], [ramp.xds::EggLaying]
#' @return a [numeric] vector of length `nPatches`
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
#' @title Get the feeding rate for `newMYname` (**MY**)
#' @description Return the stored realized feeding rate per patch.
#' @inheritSection newMYname Parameters
#' @inheritParams ramp.xds::get_f
#' @seealso [newMYname]
#' @return a [numeric] vector of length `nPatches`
#' @export
get_f.newMYname <- function(xds_obj, s = 1) {
  xds_obj$MY_obj[[s]]$f
}

#' @title Get the human blood fraction for `newMYname` (**MY**)
#' @description Return the stored realized human blood fraction per patch.
#' @inheritSection newMYname Parameters
#' @inheritParams ramp.xds::get_q
#' @seealso [newMYname]
#' @return a [numeric] vector of length `nPatches`
#' @export
get_q.newMYname <- function(xds_obj, s = 1) {
  xds_obj$MY_obj[[s]]$q
}

#' @title Get the adult mosquito mortality rate for `newMYname` (**MY**)
#' @description Return the stored realized adult mosquito mortality rate per patch.
#' @inheritSection newMYname Parameters
#' @inheritParams ramp.xds::get_g
#' @seealso [newMYname]
#' @return a [numeric] vector of length `nPatches`
#' @export
get_g.newMYname <- function(xds_obj, s = 1) {
  xds_obj$MY_obj[[s]]$g
}

#' @title Get the patch emigration rate for `newMYname` (**MY**)
#' @description Return the stored realized patch emigration rate per patch.
#' @inheritSection newMYname Parameters
#' @inheritParams ramp.xds::get_sigma
#' @seealso [newMYname]
#' @return a [numeric] vector of length `nPatches`
#' @export
get_sigma.newMYname <- function(xds_obj, s = 1) {
  xds_obj$MY_obj[[s]]$sigma
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Make a module object for `newMYname` (**MY**)
#' @description Create a `newMYname` module object and initialize its parameters.
#' Inspect parameters with [get_MY_pars.newMYname] and change supported values
#' with [change_MY_pars.newMYname]. Values in `options` override parameter arguments.
#' @seealso [get_MY_pars.newMYname], [change_MY_pars.newMYname]
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

#' @title Set up the module object for `newMYname` (**MY**)
#' @description Set up a `newMYname` module object for the **MY** component
#' and the \eqn{s^{th}} vector species.
#' + Create the module object with [make_MY_obj_newMYname], passing `nPatches` and `options`.
#' + Set initial values using [ramp.xds::setup_MY_inits], passing `options` to override defaults.
#' + Return the updated **`xds`** model object.
#' @note The module object is stored at `xds_obj$MY_obj[[s]]`.
#' @inheritParams ramp.xds::setup_MY_obj
#' @return an **`xds`** model object
#' @export
setup_MY_obj.newMYname <- function(MYname, xds_obj, options = list(), s = 1) {
  xds_obj$MY_obj[[s]] <- make_MY_obj_newMYname(xds_obj$nPatches, options)
  xds_obj <- setup_MY_inits(xds_obj, options = options, s = s)
  xds_obj
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Get parameters for `newMYname` (**MY**)
#' @description Return the parameters for the `newMYname` module as a named list.
#' @inheritParams ramp.xds::get_MY_pars
#' @seealso [make_MY_obj_newMYname], [change_MY_pars.newMYname]
#' @return a named [list] of parameter values
#' @export
get_MY_pars.newMYname <- function(xds_obj, s = 1) {
  with(xds_obj$MY_obj[[s]], list(p1 = p1))
}

#' @title Change parameters for `newMYname` (**MY**)
#' @description Change parameter values named in `options` for the `newMYname`
#' module, leaving other values unchanged.
#' This example supports changes to `p1`.
#' Call [ramp.xds::MBionomics] at time zero after changing parameters;
#' do not call [ramp.xds::MEffectSizes].
#' @inheritParams ramp.xds::change_MY_pars
#' @seealso [make_MY_obj_newMYname], [get_MY_pars.newMYname]
#' @return an **`xds`** model object
#' @export
change_MY_pars.newMYname <- function(xds_obj, options = list(), s = 1) {
  if ("p1" %in% names(options)) {
    xds_obj$MY_obj[[s]]$p1 <- checkIt(options$p1, xds_obj$nPatches)
  }
  y <- ramp.xds::get_inits(xds_obj, flatten = TRUE)
  xds_obj <- ramp.xds::MBionomics(0, y, xds_obj, s)
  xds_obj
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Make initial values for `newMYname` (**MY**)
#' @description Create a named list of initial values for the `newMYname` module.
#' Values in `options` override the defaults.
#' @details Initial values must be numeric scalars or have length `nPatches`.
#' Scalars are expanded to that length.
#' @param nPatches number of patches
#' @param options named values that override defaults
#' @param M initial adult mosquito density
#' @param Y initial infected mosquito density
#' @return a named [list] of initial values
#' @seealso [setup_MY_inits.newMYname], [change_MY_inits.newMYname]
#' @export
make_MY_inits_newMYname <- function(nPatches, options = list(), M = 10, Y = 1) {
  with(options, {
    M <- checkIt(M, nPatches)
    Y <- checkIt(Y, nPatches)
    list(M = M, Y = Y)
  })
}

#' @title Set up initial values for `newMYname` (**MY**)
#' @description Set up initial values for the `newMYname` module using
#' [make_MY_inits_newMYname]. Values in `options` override the defaults.
#' @note Initial values are stored at `xds_obj$MY_obj[[s]]$inits`.
#' @inheritParams ramp.xds::setup_MY_inits
#' @seealso [make_MY_inits_newMYname], [change_MY_inits.newMYname]
#' @return an **`xds`** model object
#' @export
setup_MY_inits.newMYname <- function(xds_obj, options = list(), s = 1) {
  xds_obj$MY_obj[[s]]$inits <-
    make_MY_inits_newMYname(xds_obj$nPatches, options)
  xds_obj
}

#' @title Change initial values for `newMYname` (**MY**)
#' @description Change initial values named in `options` for the `newMYname`
#' module, leaving other values unchanged.
#' @note Stored initial values for this module are at
#' `xds_obj$MY_obj[[s]]$inits`.
#' @inheritParams ramp.xds::change_MY_inits
#' @seealso [make_MY_inits_newMYname], [setup_MY_inits.newMYname]
#' @return an **`xds`** model object
#' @export
change_MY_inits.newMYname <- function(xds_obj, options = list(), s = 1) {
  xds_obj$MY_obj[[s]]$inits[names(options)] <- options
  xds_obj
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Get the initial adult mosquito density for `newMYname` (**MY**)
#' @description Return the initial adult mosquito density per patch.
#' @inheritSection newMYname State Variables
#' @inheritParams ramp.xds::get_M
#' @seealso [newMYname], [ramp.xds::get_MY_inits]
#' @return a [numeric] vector of length `nPatches`
#' @export
get_M.newMYname <- function(xds_obj, s = 1) {
  ramp.xds::get_MY_inits(xds_obj, s)$M
}

#' @title Get the initial infective mosquito density for `newMYname` (**MY**)
#' @description Return the initial infective mosquito density per patch.
#' @inheritSection newMYname State Variables
#' @inheritParams ramp.xds::get_Z
#' @seealso [newMYname], [ramp.xds::get_MY_inits]
#' @return a [numeric] vector of length `nPatches`
#' @export
get_Z.newMYname <- function(xds_obj, s = 1) {
  ramp.xds::get_MY_inits(xds_obj, s)$Z
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Set up variable indices for `newMYname` (**MY**)
#' @description Assign consecutive variable indices for the `newMYname` module,
#' starting after `xds_obj$max_ix`, and update the maximum index.
#' Allocate blocks in the order `M_ix`, `Y_ix`.
#' Each block has length `nPatches`.
#' @note Indices are stored at `xds_obj$MY_obj[[s]]$ix`.
#' The maximum index is stored at `xds_obj$max_ix`.
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
#' @title Parse outputs for `newMYname` (**MY**)
#' @description Extract the model state outputs for the `newMYname` module
#' from `outputs` and return them by name in a list.
#' @note Variable indices are stored at `xds_obj$MY_obj[[s]]$ix`.
#' @inheritParams ramp.xds::parse_MY_orbits
#' @return a [list]
#' @export
parse_MY_orbits.newMYname <- function(outputs, xds_obj, s) {
  with(xds_obj$MY_obj[[s]]$ix, {
    M <- outputs[, M_ix]
    Y <- outputs[, Y_ix]
    list(M = M, Y = Y)
  })
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Check `newMYname` (**MY**)
#' @description Return the model unchanged; no checks are implemented in this template.
#' @inheritParams ramp.xds::check_MY
#' @return an **`xds`** model object
#' @export
check_MY.newMYname <- function(xds_obj, s) {
  xds_obj
}

#' @title Get the skill set for `newMYname` (**MY**)
#' @description Return the capabilities and compatibilities of the `newMYname` module as a list.
#' @inheritParams ramp.xds::get_skillset_MY
#' @return the **MY** module skill set, as a list
#' @export
get_skillset_MY.newMYname <- function(xds_obj, s) {
  list(dts=FALSE, delay=FALSE)
}

