testthat::test_that("XH setup passes the selected host's stratum count to its constructor", {
  modules <- c("SEIR", "SEIRV", "SIR", "SIRS", "SEIS", "SEISd",
               "SIP", "SIPd", "garki")

  for (module in modules) {
    method_name <- paste0("setup_XH_obj.", module)
    method <- getFromNamespace(method_name, "ramp.library")
    env <- new.env(parent = baseenv())
    env[[method_name]] <- eval(
      call("function", formals(method), body(method)),
      envir = env
    )

    # Isolate constructor arguments from unrelated interface setup.
    env[[paste0("make_XH_obj_", module)]] <- function(nStrata, options) {
      list(nStrata = nStrata, options = options)
    }
    env$get_skillset_XH <- function(...) stop("XH setup must not call get_skillset_XH")
    env$setup_XH_inits <- function(xds_obj, ...) xds_obj
    env$setup_mass_treatment <- function(xds_obj, ...) xds_obj
    env$setup_demog_obj <- function(xds_obj, ...) xds_obj
    for (helper in c(
      "setup_timespent", "setup_blood_search_weights",
      "setup_time_away", "setup_travel_eir"
    )) {
      env[[helper]] <- function(name, xds_obj, ...) xds_obj
    }

    model <- list(nStrata = c(2L, 5L), XH_obj = list(list(), list()))
    options <- list(b = 0.3)
    for (i in seq_along(model$nStrata)) {
      result <- env[[method_name]](
        module, residence = seq_len(model$nStrata[i]),
        xds_obj = model, options = options, i = i
      )
      testthat::expect_identical(result$XH_obj[[i]]$nStrata, model$nStrata[i],
                                info = paste(module, i))
      testthat::expect_identical(result$XH_obj[[i]]$options, options,
                                info = paste(module, i))
      testthat::expect_false("skill_set" %in% names(result$XH_obj[[i]]),
                            info = paste(module, i))
      testthat::expect_identical(result$XH_obj[[3L - i]], model$XH_obj[[3L - i]],
                                info = paste(module, i))
    }
  }
})
