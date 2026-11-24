testthat::test_that("library changers refresh bionomics without applying effect sizes", {
  files <- c("adult-MY-BQ.R", "adult-MY-BQm.R", "adult-MY-SEI.R",
             "adult-MY-BRQS.R", "aquatic-L-stages.R")
  for (file in files) {
    definitions <- parse(testthat::test_path("..", "..", "R", file))
    definition <- Filter(function(e) {
      is.call(e) && as.character(e[[1]])[1] %in% c("<-", "=") &&
        grepl("^change_(MY|L)_pars\\.", as.character(e[[2]])[1])
    }, definitions)
    testthat::expect_length(definition, 1L)
    name <- as.character(definition[[1]][[2]])
    aquatic <- startsWith(name, "change_L")
    component <- if (aquatic) "L_obj" else "MY_obj"
    env <- new.env(parent = baseenv())
    eval(definition[[1]], env)
    parameters <- c("f", "q", "g", "sigma", "sigma_b", "sigma_q", "mu", "nu",
                    "eip", "eggsPerBatch", "psi_E", "phi_E", "phi_P", "psi_P",
                    outer(c("psi", "xi", "phi", "theta"), 1:4, paste0))
    object <- setNames(as.list(rep(1, length(parameters))), parameters)
    for (par in parameters) object[[paste0(par, "_obj")]] <- setNames(list(1), par)
    object$nPatches <- 2L
    model <- list(nHabitats = 2L, nPatches = 2L)
    model[[component]] <- list(object, object)
    env$checkIt <- ramp.xds::checkIt
    env$get_inits <- function(xds_obj, flatten) {
      testthat::expect_true(flatten)
      c(state = 7)
    }
    calls <- 0L
    env[[if (aquatic) "LBionomics" else "MBionomics"]] <- function(t, y, xds_obj, s) {
      calls <<- calls + 1L
      testthat::expect_identical(t, 0)
      testthat::expect_identical(y, c(state = 7))
      testthat::expect_identical(s, 2L)
      updated <- xds_obj[[component]][[s]]
      value <- if (aquatic) updated$psi_E else if (name == "change_MY_pars.BRQS") updated$f else updated$f_obj$f
      testthat::expect_equal(value, if (name == "change_MY_pars.BRQS" || aquatic) c(2, 2) else 2)
      xds_obj$refreshed <- TRUE
      xds_obj
    }
    env$MEffectSizes <- env$LEffectSizes <- function(...) stop("Unexpected effect-size call")
    result <- env[[name]](model, if (aquatic) list(psi_E = 2) else list(f = 2), 2L)
    testthat::expect_true(result$refreshed)
    testthat::expect_identical(calls, 1L)
    testthat::expect_identical(result[[component]][[1]], object)
  }
})

testthat::test_that("BRQS changes its direct parameter fields and returns the model", {
  env <- new.env(parent = baseenv())
  definitions <- parse(testthat::test_path("..", "..", "R", "adult-MY-BRQS.R"))
  for (definition in definitions) {
    if (is.call(definition) && as.character(definition[[1]])[1] %in% c("<-", "=") &&
        as.character(definition[[2]])[1] %in% c("change_MY_pars.BRQS", "MBionomics.BRQS")) {
      eval(definition, env)
    }
  }
  env$checkIt <- ramp.xds::checkIt
  env$get_inits <- function(...) numeric()
  env$MBionomics <- env$MBionomics.BRQS
  baseline <- c("f", "q", "nu", "g_s", "g_b", "g_r", "g_q",
                "mu_s", "mu_b", "mu_q", "sigma_b", "sigma_q", "sigma_s",
                "theta", "xi", "eta")
  object <- c(list(nPatches = 2L, eip = 12, eggsPerBatch = 60),
              setNames(rep(list(c(1, 1)), length(baseline)), baseline))
  model <- list(nPatches = 2L, MY_obj = list(object, object))
  options <- setNames(rep(list(c(2, 3)), length(baseline)), baseline)
  options$eip <- 8
  options$eggsPerBatch <- 40
  result <- env$change_MY_pars.BRQS(model, options, 2L)
  testthat::expect_identical(result$MY_obj[[1]], object)
  for (par in baseline) testthat::expect_identical(result$MY_obj[[2]][[par]], c(2, 3))
  testthat::expect_identical(result$MY_obj[[2]]$eip, 8)
  testthat::expect_identical(result$MY_obj[[2]]$eggsPerBatch, 40)
  testthat::expect_false(any(grepl("_t$|^es_", names(result$MY_obj[[2]]))))
  unchanged <- env$change_MY_pars.BRQS(model)
  testthat::expect_identical(unchanged$MY_obj[[1]]$f, object$f)
})
