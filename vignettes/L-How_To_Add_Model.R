## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title The `newLname` module for the **L** component
#' @section State Variables:
#' Describe the module's actual states and their habitat-level interpretation.
#' @section Dynamical Terms:
#' Describe the shared inputs used by the implemented equations, including `eta`.
#' @section Parameters:
#' Describe the parameters used by the implemented model and their units.
#' @section ODE:
#' Replace this placeholder with the continuous-time model equations.
#' @section DTS:
#' Replace this placeholder with the discrete-time updating equations.
#' @section Emergence:
#' Replace this placeholder with the module's adult female emergence equation.
#' [F_emerge.newLname] returns habitat-level emergence before aggregation
#' to patches by [ramp.xds::Emergence].
#' @name newLname
NULL

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Compute the derivatives for `newLname` (**L**)
#' @description Compute the derivatives for the `newLname` module:
#' \deqn{dL/dt = 0}
#' Replace this equation and the derivative code together.
#' @inheritSection newLname State Variables
#' @inheritSection newLname Dynamical Terms
#' @inheritSection newLname Parameters
#' @inheritSection newLname ODE
#' @inheritParams ramp.xds::dLdt
#' @note Derivatives for `L` are returned in habitat order.
#' The egg input is read from `xds_obj$terms$eta[[s]]`.
#' @return a [numeric] vector of derivatives in the module's state-variable
#' order, with length `nHabitats`
#' @export
dLdt.newLname <- function(t, y, xds_obj, s) {
  eta <- xds_obj$terms$eta[[s]]

  with(ramp.xds::get_L_vars(y, xds_obj, s), {
    with(xds_obj$L_obj[[s]], {
      # Replace this expression with the model's derivative equation.
      dL <- rep(0, xds_obj$nHabitats)
      derivatives <- c(dL)
      return(derivatives)
    })
  })
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Update states for `newLname` (**L**)
#' @description Update the state variables for the `newLname` module:
#' \deqn{L_{t+1} = 0}
#' Replace this equation and the updating code together.
#' @inheritSection newLname State Variables
#' @inheritSection newLname Dynamical Terms
#' @inheritSection newLname Parameters
#' @inheritSection newLname DTS
#' @inheritParams ramp.xds::Update_Lt
#' @note Updated values for `L` are returned in habitat order.
#' The egg input is read from `xds_obj$terms$eta[[s]]`.
#' @return a [numeric] vector of updated states in the module's state-variable
#' order, with length `nHabitats`
#' @export
Update_Lt.newLname <- function(t, y, xds_obj, s) {
  eta <- xds_obj$terms$eta[[s]]

  with(ramp.xds::get_L_vars(y, xds_obj, s), {
    with(xds_obj$L_obj[[s]], {
      # Replace this expression with the model's discrete update.
      Lt <- rep(0, xds_obj$nHabitats)
      states <- c(Lt)
      return(states)
    })
  })
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Get the state variables for `newLname` (**L**)
#' @description Extract the state variables for the `newLname` module
#' from `y` and return them by name in a list for use during model solving.
#' @note Variable indices are stored at `xds_obj$L_obj[[s]]$ix`.
#' @inheritParams ramp.xds::get_L_vars
#' @return a named [list] of state variables and module-specific derived quantities
#' @export
get_L_vars.newLname <- function(y, xds_obj, s) {
  with(xds_obj$L_obj[[s]]$ix, {
    L <- y[L_ix]
    return(list(L = L))
  })
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Compute the baseline p1 parameter
#' @inheritParams ramp.xds::F_psi
#' @return a numeric vector of length `nHabitats`
#' @export
F_p1 <- function(t, xds_obj, s) {
  UseMethod("F_p1", xds_obj$L_obj[[s]]$p1_obj)
}

#' @rdname F_p1
#' @export
F_p1.static <- function(t, xds_obj, s) {
  xds_obj$L_obj[[s]]$p1_obj$p1
}

#' @title Mosquito bionomics for `newLname` (**L**)
#' @description Update baseline parameters and reset effect sizes.
#' @inheritSection newLname Parameters
#' @note Store the evaluated baseline as `p1_t` and reset `es_p1` to 1
#' per habitat at `xds_obj$L_obj[[s]]`. Leave the direct parameter `p2` unchanged.
#' @inheritParams ramp.xds::LBionomics
#' @seealso [newLname], [LEffectSizes.newLname]
#' @return an **`xds`** model object
#' @export
LBionomics.newLname <- function(t, y, xds_obj, s) {
  with(xds_obj$L_obj[[s]], {
    # Replace these assignments with the appropriate forcing functions.
    xds_obj$L_obj[[s]]$p1_t <- F_p1(t, xds_obj, s)
    xds_obj$L_obj[[s]]$p2 <- p2
    xds_obj$L_obj[[s]]$es_p1 <- rep(1, xds_obj$nHabitats)
    return(xds_obj)
  })
}

#' @title Apply effect sizes for `newLname` (**L**)
#' @description Apply vector-control effect sizes to baseline parameters.
#' @inheritSection newLname Parameters
#' @note Set `p1 = p1_t * es_p1` at `xds_obj$L_obj[[s]]`.
#' Do not evaluate ports, reset effect sizes, or change the direct parameter `p2`.
#' @inheritParams ramp.xds::LEffectSizes
#' @seealso [newLname], [LBionomics.newLname]
#' @return an **`xds`** model object
#' @export
LEffectSizes.newLname <- function(t, y, xds_obj, s) {
  with(xds_obj$L_obj[[s]], {
    xds_obj$L_obj[[s]]$p1 <- p1_t * es_p1
    return(xds_obj)
  })
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Compute adult mosquito emergence for `newLname` (**L**)
#' @description Return adult female mosquito emergence per habitat.
#' @inheritSection newLname Emergence
#' @inheritParams ramp.xds::F_emerge
#' @seealso [newLname], [ramp.xds::Emergence]
#' @return a [numeric] vector of length `nHabitats`
#' @export
F_emerge.newLname <- function(t, y, xds_obj, s) {
  with(ramp.xds::get_L_vars(y, xds_obj, s), {
    with(xds_obj$L_obj[[s]], {
      # Replace this expression with the model's emergence equation.
      alpha <- rep(0, xds_obj$nHabitats)
      return(alpha)
    })
  })
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Get the maturation rate for `newLname` (**L**)
#' @description Return the stored realized maturation rate per habitat.
#' @inheritSection newLname Parameters
#' @inheritParams ramp.xds::get_psi
#' @seealso [newLname], [ramp.xds::change_psi]
#' @return a [numeric] vector of length `nHabitats`
#' @export
get_psi.newLname <- function(xds_obj, s = 1) {
  xds_obj$L_obj[[s]]$psi
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Make a module object for `newLname` (**L**)
#' @description Create a `newLname` module object and initialize its parameters.
#' Inspect parameters with [get_L_pars.newLname] and change supported values
#' with [change_L_pars.newLname]. Values in `options` override parameter arguments.
#' @seealso [get_L_pars.newLname], [change_L_pars.newLname]
#' @param nHabitats the number of habitats in the model
#' @param options a named list of parameter values
#' @return an L module object
#' @export
make_L_obj_newLname <- function(nHabitats, options = list(), p1 = 1, p2 = 1) {
  with(options, {
    L_obj <- list(
      p1 = ramp.xds::checkIt(p1, nHabitats),
      p2 = ramp.xds::checkIt(p2, nHabitats),
      p1_t = ramp.xds::checkIt(p1, nHabitats),
      es_p1 = rep(1, nHabitats),
      p1_obj = structure(list(p1 = ramp.xds::checkIt(p1, nHabitats)),
                         class = "static"),
      nHabitats = nHabitats
    )
    class(L_obj) <- "newLname"
    return(L_obj)
  })
}

#' @title Set up the module object for `newLname` (**L**)
#' @description Set up a `newLname` module object for the **L** component
#' and the \eqn{s^{th}} vector species.
#' + Create the module object with [make_L_obj_newLname], passing `nHabitats` and `options`.
#' + Set the skill set using [ramp.xds::get_skillset_L].
#' + Set initial values using [ramp.xds::setup_L_inits], passing `options` to override defaults.
#' + Initialize habitat search weights using [ramp.xds::setup_habitat_search_weights.setup].
#' + Return the updated **`xds`** model object.
#' @note The module object is stored at `xds_obj$L_obj[[s]]`.
#' Its classes are `newLname` and `newLname_` followed by `xds_obj$xds`.
#' @inheritParams ramp.xds::setup_L_obj
#' @return an **`xds`** model object
#' @export
setup_L_obj.newLname <- function(Lname, membership, xds_obj,
                                 options = list(), s = 1) {
  L_obj <- make_L_obj_newLname(xds_obj$nHabitats, options)
  class(L_obj) <- c("newLname", paste0("newLname_", xds_obj$xds))
  xds_obj$L_obj[[s]] <- L_obj
  xds_obj$L_obj[[s]]$skill_set <- ramp.xds::get_skillset_L(xds_obj, s)
  xds_obj <- ramp.xds::setup_L_inits(xds_obj, options = options, s = s)
  xds_obj <- ramp.xds::setup_habitat_search_weights(
    "setup", xds_obj, list(membership = membership), s = s
  )
  return(xds_obj)
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Get parameters for `newLname` (**L**)
#' @description Return the parameters for the `newLname` module as a named list.
#' @seealso [make_L_obj_newLname], [change_L_pars.newLname]
#' @inheritParams ramp.xds::get_L_pars
#' @return a named [list] of parameter values
#' @export
get_L_pars.newLname <- function(xds_obj, s = 1) {
  with(xds_obj$L_obj[[s]], {
    return(list(p1 = p1, p2 = p2))
  })
}

#' @title Change parameters for `newLname` (**L**)
#' @description Change parameter values named in `options` for the `newLname`
#' module, leaving other values unchanged.
#' Update the `p1` port and the direct parameter `p2`.
#' The effective `p1` value is set by [LEffectSizes.newLname].
#' Call [ramp.xds::LBionomics] at time zero after the updates,
#' without calling [ramp.xds::LEffectSizes].
#' @seealso [make_L_obj_newLname], [get_L_pars.newLname]
#' @inheritParams ramp.xds::change_L_pars
#' @return an **`xds`** model object
#' @export
change_L_pars.newLname <- function(xds_obj, options = list(), s = 1) {
  nHabitats <- xds_obj$nHabitats
  if ("p1" %in% names(options)) {
    xds_obj$L_obj[[s]]$p1_obj$p1 <- ramp.xds::checkIt(options$p1, nHabitats)
  }
  if ("p2" %in% names(options)) {
    xds_obj$L_obj[[s]]$p2 <- ramp.xds::checkIt(options$p2, nHabitats)
  }
  y <- ramp.xds::get_inits(xds_obj, flatten = TRUE)
  xds_obj <- ramp.xds::LBionomics(0, y, xds_obj, s)
  return(xds_obj)
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Make initial values for `newLname` (**L**)
#' @description Create a named list of initial values for the `newLname` module.
#' Values in `options` override the defaults.
#' @details Initial values must be numeric scalars or have length `nHabitats`.
#' Scalars are expanded to that length.
#' @param nHabitats the number of habitats in the model
#' @param options a named list of initial values
#' @param L the initial value for state `L`
#' @return a named [list] of initial values
#' @seealso [setup_L_inits.newLname], [change_L_inits.newLname]
#' @export
make_L_inits_newLname <- function(nHabitats, options = list(), L = 1) {
  with(options, {
    L <- ramp.xds::checkIt(L, nHabitats)
    return(list(L = L))
  })
}

#' @title Set up initial values for `newLname` (**L**)
#' @description Set up initial values for the `newLname` module using
#' [make_L_inits_newLname]. Values in `options` override the defaults.
#' @note Initial values are stored at `xds_obj$L_obj[[s]]$inits`.
#' @inheritParams ramp.xds::setup_L_inits
#' @seealso [make_L_inits_newLname], [change_L_inits.newLname]
#' @return an **`xds`** model object
#' @export
setup_L_inits.newLname <- function(xds_obj, options = list(), s = 1) {
  xds_obj$L_obj[[s]]$inits <- make_L_inits_newLname(
    xds_obj$nHabitats, options
  )
  return(xds_obj)
}

#' @title Change initial values for `newLname` (**L**)
#' @description Change initial values named in `options` for the `newLname`
#' module, leaving other values unchanged.
#' @note Stored initial values for this module are at
#' `xds_obj$L_obj[[s]]$inits`.
#' @inheritParams ramp.xds::change_L_inits
#' @seealso [make_L_inits_newLname], [setup_L_inits.newLname]
#' @return an **`xds`** model object
#' @export
change_L_inits.newLname <- function(xds_obj, options = list(), s = 1) {
  with(options, {
    if (!is.null(L)) {
      xds_obj$L_obj[[s]]$inits$L <- ramp.xds::checkIt(
        L, xds_obj$nHabitats
      )
    }
  })
  return(xds_obj)
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Set up variable indices for `newLname` (**L**)
#' @description Assign consecutive variable indices for the `newLname` module,
#' starting after `xds_obj$max_ix`, and update the maximum index.
#' Allocate the `L_ix` block with length `nHabitats`.
#' @note Indices are stored at `xds_obj$L_obj[[s]]$ix`.
#' The maximum index is stored at `xds_obj$max_ix`.
#' @inheritParams ramp.xds::setup_L_ix
#' @return an **`xds`** model object
#' @importFrom utils tail
#' @export
setup_L_ix.newLname <- function(xds_obj, s) {
  L_ix <- seq(from = xds_obj$max_ix + 1, length.out = xds_obj$nHabitats)
  xds_obj$max_ix <- utils::tail(L_ix, 1)
  xds_obj$L_obj[[s]]$ix <- list(L_ix = L_ix)
  return(xds_obj)
}


## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Parse outputs for `newLname` (**L**)
#' @description Extract the model state outputs for the `newLname` module
#' from `outputs` and return them by name in a list.
#' @note Variable indices are stored at `xds_obj$L_obj[[s]]$ix`.
#' @inheritParams ramp.xds::parse_L_orbits
#' @return a [list]
#' @export
parse_L_orbits.newLname <- function(outputs, xds_obj, s) {
  with(xds_obj$L_obj[[s]]$ix, {
    L <- outputs[, L_ix]
    return(list(L = L))
  })
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Check `newLname` (**L**)
#' @description Return the model unchanged; no checks are implemented in this template.
#' @inheritParams ramp.xds::check_L
#' @return an **`xds`** model object
#' @export
check_L.newLname <- function(xds_obj, s) {
  return(xds_obj)
}

#' @title Get the skill set for `newLname` (**L**)
#' @description Return the capabilities and compatibilities of the `newLname` module as a list.
#' @inheritParams ramp.xds::get_skillset_L
#' @return the **L** module skill set, as a list
#' @export
get_skillset_L.newLname <- function(xds_obj, s) {
  list(dts = FALSE, xde = TRUE, delay = FALSE)
}

