testthat::test_that("garki initial-value changes update only the selected host's stored values", {
  expressions <- parse(testthat::test_path("..", "..", "R", "malaria-garki.R"))
  env <- new.env(parent = baseenv())
  for (name in c("make_XH_inits_garki", "change_XH_inits.garki")) {
    definition <- Filter(function(expr) {
      is.call(expr) &&
        as.character(expr[[1]])[1] %in% c("<-", "=") &&
        identical(expr[[2]], as.name(name))
    }, expressions)
    testthat::expect_length(definition, 1L)
    eval(definition[[1]], env)
  }
  env$checkIt <- ramp.xds::checkIt
  model <- list(
    nStrata = c(2L, 5L),
    XH_obj = list(
      list(inits = env$make_XH_inits_garki(2L), parameter = 1),
      list(inits = env$make_XH_inits_garki(5L), parameter = 2)
    )
  )

  for (i in seq_along(model$XH_obj)) {
    testthat::expect_identical(env$change_XH_inits.garki(model, i = i), model)
    original <- model$XH_obj[[i]]$inits
    for (state in names(original)) {
      replacement <- seq_len(model$nStrata[i]) + 10
      expected <- model
      expected$XH_obj[[i]]$inits[[state]] <- replacement
      result <- env$change_XH_inits.garki(
        model, options = setNames(list(replacement), state), i = i
      )
      testthat::expect_identical(result, expected, info = paste(i, state))
      testthat::expect_false("XH_inits" %in% names(result))
    }
    options <- lapply(original, function(value) value + 1)
    expected <- model
    expected$XH_obj[[i]]$inits <- options
    testthat::expect_identical(env$change_XH_inits.garki(model, options, i), expected)
  }
})
