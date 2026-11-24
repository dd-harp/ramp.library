testthat::test_that("SEIS discrete formulas use ar rather than an undefined nr", {
  env <- new.env(parent = baseenv())
  for (definition in parse(testthat::test_path("..", "..", "R", "human-XH-SEIS.R"))) {
    if (is.call(definition) && as.character(definition[[1]])[1] %in% c("<-", "=") &&
        as.character(definition[[2]])[1] %in% c("steady_state_X_SEIS.dts", "Update_XHt.SEIS")) {
      eval(definition, env)
    }
  }
  object <- list(nu = 0.3, Hpar = list())
  model <- list(XH_obj = list(list(), object), terms = list(AR = list(0, 0.2)))
  testthat::expect_false("nr" %in% all.names(body(env$Update_XHt.SEIS)))
  testthat::expect_true(any(grepl("It <- ar * I", deparse(body(env$Update_XHt.SEIS)), fixed = TRUE)))
  foi <- c(0.1, 0.3)
  H <- c(100, 200)
  ar <- exp(-foi)
  denominator <- (1-ar)+ar-object$nu+ar*object$nu*(1-ar)
  I <- ar*H*(1-object$nu)/denominator
  E <- ar*H*(1-ar)/denominator
  testthat::expect_equal(env$steady_state_X_SEIS.dts(foi, H, model, 2),
                         c(S = H-I-E, E = E, I = I))
})

testthat::test_that("XH steady states dispatch on module and then equation system", {
  files <- c(SEIR = "host-xde-SEIR.R", SEIRV = "host-xde-SEIRV.R",
             SIR = "host-xde-SIR.R", SIRS = "host-xde-SIRS.R",
             SEIS = "human-XH-SEIS.R", SIP = "malaria-XH-SIP.R",
             SIPd = "malaria-XH-SIPd.R")
  env <- new.env(parent = baseenv())
  env$steady_state_X <- function(foi, H, xds_obj, i = 1) {
    UseMethod("steady_state_X", xds_obj$XH_obj[[i]])
  }
  environment(env$steady_state_X) <- env
  for (module in names(files)) {
    for (definition in parse(testthat::test_path("..", "..", "R", files[[module]]))) {
      if (is.call(definition) && as.character(definition[[1]])[1] %in% c("<-", "=") &&
          is.symbol(definition[[2]]) &&
          startsWith(as.character(definition[[2]]), "steady_state_X")) {
        eval(definition, env)
      }
    }
    router <- paste0("steady_state_X_", module)
    registerS3method("steady_state_X", module, env[[paste0("steady_state_X.", module)]],
                     envir = env)
    testthat::expect_false(exists(paste0("steady_state_X.", module, "_dts"), env,
                                 inherits = FALSE))
    for (system in c("xde", "dts")) {
      target <- paste0(router, ".", system)
      original <- env[[target]]
      env[[target]] <- function(foi, H, xds_obj, i = 1) {
        list(foi = foi, H = H, i = i, object = xds_obj$XH_obj[[i]])
      }
      registerS3method(router, system, env[[target]], envir = env)
      for (i in 1:2) {
        model <- list(XH_obj = list(structure(list(id = 1), class = module),
                                    structure(list(id = 2), class = module)),
                      xds = structure(system, class = system))
        env$model <- model
        env$i <- i
        result <- evalq(steady_state_X(c(0.1, 0.2), c(100, 200), model, i), env)
        testthat::expect_identical(result, list(foi = c(0.1, 0.2), H = c(100, 200),
                                               i = i, object = model$XH_obj[[i]]))
      }
      env[[target]] <- original
    }
    model$xds <- structure("unsupported", class = "unsupported")
    env$model <- model
    testthat::expect_error(evalq(steady_state_X(0.1, 100, model), env), "no applicable method")
  }
})

testthat::test_that("SIPd discrete steady state accepts FoI and selected host species", {
  env <- new.env(parent = baseenv())
  for (definition in parse(testthat::test_path("..", "..", "R", "malaria-XH-SIPd.R"))) {
    if (is.call(definition) && as.character(definition[[1]])[1] %in% c("<-", "=") &&
        identical(as.character(definition[[2]]), "steady_state_X_SIPd.dts")) {
      eval(definition, env)
    }
  }
  object <- list(eta = 0.2, xi = 0.03, rho = 0.1, r = 0.05)
  model <- list(XH_obj = list(list(), object))
  foi <- c(0.1, 0.3)
  H <- c(100, 200)
  result <- env$steady_state_X_SIPd.dts(foi, H, model, 2)
  expected <- with(object, {
    ar <- exp(-foi)
    I <- (ar*H*eta*(1-rho))/((eta+xi)*(ar*(1-r)+(r+xi)+ar*(eta*(r-1)+r+xi)))
    P <- (H*(xi*(ar+r*(1-ar)+xi)+ar*(1+xi)*rho))/
      ((eta+xi)*(ar*(1-r)+(r+xi)+ar*(eta*(r-1)+r+xi)))
    c(S = H-I-P, I = I, P = P)
  })
  testthat::expect_identical(result, expected)
})
