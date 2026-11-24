testthat::test_that("initial-value changers preserve stored state structure and other species", {
  modules <- c(
    BQ = "adult-MY-BQ.R", BQm = "adult-MY-BQm.R",
    SEIR = "host-xde-SEIR.R", SEIRV = "host-xde-SEIRV.R",
    SIR = "host-xde-SIR.R", SIRS = "host-xde-SIRS.R",
    SEIS = "human-XH-SEIS.R", SEISd = "human-XH-SEISd.R",
    SIPd = "malaria-XH-SIPd.R"
  )
  for (module in names(modules)) {
    component <- if (module %in% c("BQ", "BQm")) "MY" else "XH"
    constructor <- paste0("make_", component, "_inits_", module)
    changer <- paste0("change_", component, "_inits.", module)
    env <- new.env(parent = baseenv())
    expressions <- parse(testthat::test_path("..", "..", "R", modules[[module]]))
    for (name in c(constructor, changer)) {
      definition <- Filter(function(expr) {
        is.call(expr) && as.character(expr[[1]])[1] %in% c("<-", "=") &&
          identical(expr[[2]], as.name(name))
      }, expressions)
      testthat::expect_length(definition, 1L)
      eval(definition[[1]], env)
    }
    env$checkIt <- ramp.xds::checkIt
    env$get_XH_inits <- function(xds_obj, i) xds_obj$XH_obj[[i]]$inits
    states <- setdiff(names(formals(env[[constructor]])),
                      c("nStrata", "nPatches", "options"))
    field <- paste0(component, "_obj")
    for (dimension in c(2L, 5L)) {
      inits <- env[[constructor]](dimension)
      model <- list(nPatches = dimension, nStrata = rep(dimension, 2))
      model[[field]] <- list(
        list(inits = inits, parameter = 1),
        list(inits = inits, parameter = 2)
      )
      for (i in seq_along(model[[field]])) {
        testthat::expect_identical(env[[changer]](model, list(), i), model,
                                  info = module)
        for (state in states) {
          options <- setNames(list(seq_len(dimension) + 10), state)
          expected <- model
          expected[[field]][[i]]$inits[[state]] <- options[[state]]
          result <- env[[changer]](model, options, i)
          testthat::expect_identical(result, expected, info = paste(module, i, state))
        }
        options <- setNames(lapply(seq_along(states), function(value) {
          rep(as.numeric(value + 10), dimension)
        }), states)
        expected <- model
        expected[[field]][[i]]$inits[states] <- options
        testthat::expect_identical(env[[changer]](model, options, i), expected,
                                  info = module)
      }
    }
  }
})
