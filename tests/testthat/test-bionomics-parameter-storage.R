testthat::test_that("MY ports apply all effects once and getters read realized values", {
  for (module in c("BQ", "BQm", "SEI")) {
    env <- new.env(parent = asNamespace("ramp.library"))
    for (file in c(paste0("adult-MY-", module, ".R"), "mosquito_demography.R")) {
      for (definition in parse(testthat::test_path("..", "..", "R", file))) {
        if (is.call(definition) && as.character(definition[[1]])[1] %in% c("<-", "=")) {
          eval(definition, env)
        }
      }
    }
    ports <- c("f", "q", "g", if (module == "SEI") "sigma" else c("sigma_b", "sigma_q"),
               "mu", "nu", "eip")
    if (module == "SEI") ports <- c(ports, "K_matrix")
    for (npatches in c(1L, 3L)) {
      for (s in 1:2) {
        object <- list(unported = rep(9, npatches))
        for (par in ports) {
          object[[par]] <- rep(7, npatches)
          object[[paste0(par, "_obj")]] <- setNames(list(rep(2, npatches)), par)
          env[[paste0("F_", par)]] <- local({
            parameter <- par
            function(t, xds_obj, s) {
              xds_obj$MY_obj[[s]][[paste0(parameter, "_obj")]][[parameter]] + t
            }
          })
        }
        object$Kb_matrix <- object$Kq_matrix <- matrix(0, npatches, npatches)
        if (module == "SEI") {
          object$K_matrix <- matrix(0, npatches, npatches)
          object$K_matrix_obj <- list(K_matrix = matrix(0, npatches, npatches))
          object$Upsilon_obj <- structure(list(), class = "na")
        }
        model <- list(nPatches = npatches, xds = structure("xde", class = "xde"),
                      MY_obj = list(object, object))
        baseline <- env[[paste0("MBionomics.", module)]](0, numeric(), model, s)
        neutral <- env[[paste0("MEffectSizes.", module)]](0, numeric(), baseline, s)
        for (par in ports) {
          testthat::expect_equal(neutral$MY_obj[[s]][[par]], baseline$MY_obj[[s]][[paste0(par, "_t")]])
        }
        for (par in ports) {
          testthat::expect_identical(baseline$MY_obj[[s]][[par]], object[[par]])
          testthat::expect_equal(baseline$MY_obj[[s]][[paste0(par, "_t")]],
                                object[[paste0(par, "_obj")]][[par]])
          testthat::expect_true(all(baseline$MY_obj[[s]][[paste0("es_", par)]] == 1))
          baseline$MY_obj[[s]][[paste0("es_", par)]] <- if (par == "K_matrix") matrix(1, npatches, npatches) else rep(0.5, npatches)
        }
        result <- env[[paste0("MEffectSizes.", module)]](0, numeric(), baseline, s)
        for (par in ports) {
          testthat::expect_equal(result$MY_obj[[s]][[par]],
                                baseline$MY_obj[[s]][[paste0(par, "_t")]] *
                                  baseline$MY_obj[[s]][[paste0("es_", par)]])
        }
        for (par in c("f", "q", "g", if (module == "SEI") "sigma" else "sigma_b")) {
          # Changing a multiplier after application must not change getter results.
          result$MY_obj[[s]][[paste0("es_", par)]] <- rep(0.25, npatches)
          getter <- if (par == "sigma_b") "sigma" else par
          testthat::expect_identical(env[[paste0("get_", getter, ".", module)]](result, s),
                                    result$MY_obj[[s]][[par]])
          result$MY_obj[[s]][[paste0("es_", par)]] <- rep(0.5, npatches)
        }
        testthat::expect_identical(result$MY_obj[[3L - s]], object)
        testthat::expect_identical(result$MY_obj[[s]]$unported, object$unported)
        testthat::expect_identical(env[[paste0("MEffectSizes.", module)]](0, numeric(), result, s), result)
        realized <- result$MY_obj[[s]]
        if (module == "SEI") {
          testthat::expect_equal(realized$Omega, ramp.xds::F_Omega_xde(
            realized$g, realized$sigma, realized$mu, realized$K_matrix))
        } else {
          for (state in c("b", "q")) {
            testthat::expect_equal(realized[[paste0("Omega_", state)]], ramp.xds::F_Omega_xde(
              realized$g, realized[[paste0("sigma_", state)]], realized$mu,
              realized[[paste0("K", state, "_matrix")]]))
          }
          if (module == "BQ") testthat::expect_equal(realized$phi, 1 / realized$eip)
        }
        reset <- env[[paste0("MBionomics.", module)]](1, numeric(), result, s)
        for (par in ports) testthat::expect_true(all(reset$MY_obj[[s]][[paste0("es_", par)]] == 1))
      }
    }
  }
})

testthat::test_that("non-port stages and BRQS parameters remain direct", {
  for (module in c("stages", "BRQS")) {
    env <- new.env(parent = asNamespace("ramp.library"))
    aquatic <- module == "stages"
    file <- if (aquatic) "aquatic-L-stages.R" else "adult-MY-BRQS.R"
    for (definition in parse(testthat::test_path("..", "..", "R", file))) {
      if (is.call(definition) && as.character(definition[[1]])[1] %in% c("<-", "=")) {
        eval(definition, env)
      }
    }
    component <- if (aquatic) "L_obj" else "MY_obj"
    prefix <- if (aquatic) "L" else "M"
    constructor <- env[[if (aquatic) "make_L_obj_stages" else "make_MY_obj_BRQS"]]
    testthat::expect_false(any(grepl("_t$|^es_", all.names(body(constructor)))))
    object <- if (aquatic) env$make_L_obj_stages(3L) else list(f = c(1, 2, 3), q = c(0.5, 0.6, 0.7))
    testthat::expect_false(any(grepl("_t$|^es_", names(object))))
    model <- list(nHabitats = 3L, nPatches = 3L)
    model[[component]] <- list(object, object)
    for (s in 1:2) {
      testthat::expect_identical(env[[paste0(prefix, "Bionomics.", module)]](0, numeric(), model, s), model)
      testthat::expect_identical(env[[paste0(prefix, "EffectSizes.", module)]](0, numeric(), model, s), model)
      if (!aquatic) {
        testthat::expect_identical(env$get_f.BRQS(model, s), object$f)
        testthat::expect_identical(env$get_q.BRQS(model, s), object$q)
      }
    }
  }
})
