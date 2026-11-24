testthat::test_that("BRQS setup stores initial values on the selected mosquito module", {
  expressions <- parse(testthat::test_path("..", "..", "R", "adult-MY-BRQS.R"))
  env <- new.env(parent = baseenv())
  for (name in c("make_MY_inits_BRQS", "setup_MY_inits.BRQS")) {
    definition <- Filter(function(expr) {
      is.call(expr) &&
        as.character(expr[[1]])[1] %in% c("<-", "=") &&
        identical(expr[[2]], as.name(name))
    }, expressions)
    testthat::expect_length(definition, 1L)
    eval(definition[[1]], env)
  }
  env$checkIt <- ramp.xds::checkIt
  states <- c("Su", "Bu", "Ru", "Qu", "Sy", "By", "Ry", "Qy",
              "Sz", "Bz", "Rz", "Qz")
  defaults <- c(100, 100, 1, 10, rep(0, 8))

  for (nPatches in c(2L, 5L)) {
    model <- list(
      nPatches = nPatches,
      MY_obj = list(list(parameter = 1), list(parameter = 2))
    )
    for (s in seq_along(model$MY_obj)) {
      expected <- model
      expected$MY_obj[[s]]$inits <- setNames(
        lapply(defaults, rep, times = nPatches), states
      )
      testthat::expect_identical(env$setup_MY_inits.BRQS(model, s = s), expected)

      options <- list(Su = seq_len(nPatches) + 10, Qz = 3)
      expected$MY_obj[[s]]$inits$Su <- options$Su
      expected$MY_obj[[s]]$inits$Qz <- rep(3, nPatches)
      result <- env$setup_MY_inits.BRQS(model, options, s)
      testthat::expect_identical(result, expected)
      testthat::expect_false("MYinits" %in% names(result))

      options <- setNames(lapply(seq_along(states), function(value) {
        rep(as.numeric(value), nPatches)
      }), states)
      expected$MY_obj[[s]]$inits <- options
      testthat::expect_identical(env$setup_MY_inits.BRQS(model, options, s), expected)
    }
  }
})
