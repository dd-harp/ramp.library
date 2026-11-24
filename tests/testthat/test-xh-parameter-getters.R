testthat::test_that("XH parameter getters expose constructor parameters for the selected host", {
  modules <- c(
    SEIR = "host-xde-SEIR.R",
    SEIRV = "host-xde-SEIRV.R",
    SIR = "host-xde-SIR.R",
    SIRS = "host-xde-SIRS.R",
    garki = "malaria-garki.R"
  )

  for (module in names(modules)) {
    constructor <- paste0("make_XH_obj_", module)
    getter <- paste0("get_XH_pars.", module)
    changer <- paste0("change_XH_pars.", module)
    expressions <- parse(testthat::test_path("..", "..", "R", modules[[module]]))
    env <- new.env(parent = baseenv())
    for (name in c(constructor, getter, changer)) {
      definition <- Filter(function(expr) {
        is.call(expr) &&
          as.character(expr[[1]])[1] %in% c("<-", "=") &&
          identical(expr[[2]], as.name(name))
      }, expressions)
      testthat::expect_length(definition, 1L)
      eval(definition[[1]], env)
    }
    env$checkIt <- ramp.xds::checkIt
    parameter_names <- setdiff(names(formals(env[[constructor]])), c("nStrata", "options"))
    objects <- list(env[[constructor]](2L), env[[constructor]](5L, list(b = 0.3)))
    model <- list(XH_obj = objects)

    for (i in seq_along(objects)) {
      expected <- objects[[i]][parameter_names]
      result <- env[[getter]](model, i)
      testthat::expect_setequal(names(result), parameter_names)
      testthat::expect_identical(result, expected[names(result)])
      testthat::expect_true(all(lengths(result) == c(2L, 5L)[i]))
      testthat::expect_identical(model$XH_obj, objects)

      changed <- env[[changer]](model, options = list(b = 0.6), i = i)
      testthat::expect_identical(env[[getter]](changed, i)$b, 0.6)
      testthat::expect_identical(changed$XH_obj[[3L - i]], objects[[3L - i]])
    }

    namespace <- readLines(testthat::test_path("..", "..", "NAMESPACE"))
    testthat::expect_true(paste0("S3method(get_XH_pars,", module, ")") %in% namespace)
  }
})
