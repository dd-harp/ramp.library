## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Compute derivatives for `newLname` (**L**)
#' @description Implements [ramp.xds::dLdt] for `newLname`.
#' @inheritParams ramp.xds::dLdt
#' @return a numeric vector of derivatives, in state-variable order
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
#' @description Implements [ramp.xds::Update_Lt] for `newLname`.
#' @inheritParams ramp.xds::Update_Lt
#' @return a numeric vector of updated states, in state-variable order
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
#' @title Mosquito bionomics for `newLname` (**L**)
#' @description Update baseline parameters and reset effect sizes.
#' @inheritParams ramp.xds::LBionomics
#' @return an **`xds`** model object
#' @export
LBionomics.newLname <- function(t, y, xds_obj, s) {
  with(xds_obj$L_obj[[s]], {
    # Replace these assignments with the appropriate forcing functions.
    xds_obj$L_obj[[s]]$p1_t <- p1
    xds_obj$L_obj[[s]]$p2_t <- p2
    xds_obj$L_obj[[s]]$es_p1 <- rep(1, xds_obj$nHabitats)
    xds_obj$L_obj[[s]]$es_p2 <- rep(1, xds_obj$nHabitats)
    return(xds_obj)
  })
}

#' @title Apply effect sizes for `newLname` (**L**)
#' @description Apply vector-control effect sizes to baseline parameters.
#' @inheritParams ramp.xds::LEffectSizes
#' @return an **`xds`** model object
#' @export
LEffectSizes.newLname <- function(t, y, xds_obj, s) {
  with(xds_obj$L_obj[[s]], {
    xds_obj$L_obj[[s]]$p1 <- p1_t * es_p1
    xds_obj$L_obj[[s]]$p2 <- p2_t * es_p2
    return(xds_obj)
  })
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Make the `newLname` module object
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
      p2_t = ramp.xds::checkIt(p2, nHabitats),
      es_p1 = rep(1, nHabitats),
      es_p2 = rep(1, nHabitats),
      nHabitats = nHabitats
    )
    class(L_obj) <- "newLname"
    return(L_obj)
  })
}

#' @title Set up `newLname` (**L**)
#' @description Implements [ramp.xds::setup_L_obj] for `newLname`.
#' @inheritParams ramp.xds::setup_L_obj
#' @return an **`xds`** model object
#' @export
setup_L_obj.newLname <- function(Lname, membership, xds_obj,
                                 options = list(), s = 1) {
  L_obj <- make_L_obj_newLname(xds_obj$nHabitats, options)
  class(L_obj) <- c("newLname", paste0("newLname_", xds_obj$xds))
  xds_obj$L_obj[[s]] <- L_obj
  xds_obj <- ramp.xds::setup_L_inits(xds_obj, options = options, s = s)
  xds_obj <- ramp.xds::setup_habitat_search_weights(
    "setup", xds_obj, list(membership = membership), s = s
  )
  return(xds_obj)
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Get parameters for `newLname` (**L**)
#' @inheritParams ramp.xds::get_L_pars
#' @return a named list of parameters
#' @export
get_L_pars.newLname <- function(xds_obj, s = 1) {
  with(xds_obj$L_obj[[s]], {
    return(list(p1 = p1, p2 = p2))
  })
}

#' @title Change parameters for `newLname` (**L**)
#' @inheritParams ramp.xds::change_L_pars
#' @return an **`xds`** model object
#' @export
change_L_pars.newLname <- function(xds_obj, options = list(), s = 1) {
  nHabitats <- xds_obj$nHabitats
  if ("p1" %in% names(options)) {
    xds_obj$L_obj[[s]]$p1_t <- ramp.xds::checkIt(options$p1, nHabitats)
  }
  if ("p2" %in% names(options)) {
    xds_obj$L_obj[[s]]$p2_t <- ramp.xds::checkIt(options$p2, nHabitats)
  }
  return(xds_obj)
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Set up state indices for `newLname` (**L**)
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

#' @title Extract states for `newLname` (**L**)
#' @inheritParams ramp.xds::get_L_vars
#' @return a named list of state variables
#' @export
get_L_vars.newLname <- function(y, xds_obj, s) {
  with(xds_obj$L_obj[[s]]$ix, {
    L <- y[L_ix]
    return(list(L = L))
  })
}

## ----eval=TRUE, purl=TRUE-----------------------------------------------------
#' @title Make initial values for `newLname` (**L**)
#' @param nHabitats the number of habitats in the model
#' @param options a named list of initial values
#' @param L the initial value for state `L`
#' @return a named list of initial values
#' @export
make_L_inits_newLname <- function(nHabitats, options = list(), L = 1) {
  with(options, {
    L <- ramp.xds::checkIt(L, nHabitats)
    return(list(L = L))
  })
}

#' @title Set up initial values for `newLname` (**L**)
#' @inheritParams ramp.xds::setup_L_inits
#' @return an **`xds`** model object
#' @export
setup_L_inits.newLname <- function(xds_obj, options = list(), s = 1) {
  xds_obj$L_obj[[s]]$inits <- make_L_inits_newLname(
    xds_obj$nHabitats, options
  )
  return(xds_obj)
}

#' @title Change initial values for `newLname` (**L**)
#' @inheritParams ramp.xds::change_L_inits
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
#' @title Compute emergence for `newLname` (**L**)
#' @inheritParams ramp.xds::F_emerge
#' @return a numeric vector of length `nHabitats`
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
#' @title Parse outputs for `newLname` (**L**)
#' @inheritParams ramp.xds::parse_L_orbits
#' @return a named list of state-variable orbits
#' @export
parse_L_orbits.newLname <- function(outputs, xds_obj, s) {
  with(xds_obj$L_obj[[s]]$ix, {
    L <- outputs[, L_ix]
    return(list(L = L))
  })
}

