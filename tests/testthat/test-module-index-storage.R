testthat::test_that("BRQS and SIPd index setup and consumers share module-level storage", {
  for (module in c("BRQS", "SIPd")) {
    component <- if (module == "BRQS") "MY" else "XH"
    file <- if (module == "BRQS") "adult-MY-BRQS.R" else "malaria-XH-SIPd.R"
    env <- new.env(parent = baseenv())
    env$tail <- utils::tail
    env$get_f <- function(...) 2
    for (definition in parse(testthat::test_path("..", "..", "R", file))) {
      if (is.call(definition) && as.character(definition[[1]])[1] %in% c("<-", "=")) {
        eval(definition, env)
      }
    }
    env$F_ni <- function(vars, XH_obj) env$F_ni.SIPd(vars, XH_obj)
    env$F_prevalence <- function(vars, XH_obj) env$F_prevalence.SIPd(vars, XH_obj)
    for (species in 1:2) {
      object <- list(c = 0.5, nu = 3, eggsPerBatch = 10)
      field <- paste0(component, "_obj")
      model <- list(max_ix = 7L, nPatches = 3L, nStrata = c(1L, 5L))
      model[[field]] <- list(object, object)
      result <- env[[paste0("setup_", component, "_ix.", module)]](model, species)
      ix <- result[[field]][[species]]$ix
      dimension <- if (module == "BRQS") 3L else model$nStrata[species]
      states <- if (module == "BRQS") {
        c("Su", "Bu", "Ru", "Qu", "Sy", "By", "Ry", "Qy", "Sz", "Bz", "Rz", "Qz")
      } else c("S", "I", "P", "treated")
      testthat::expect_identical(names(ix), paste0(states, "_ix"))
      testthat::expect_true(all(lengths(ix) == dimension))
      testthat::expect_equal(unlist(ix, use.names = FALSE),
                             seq(8, length.out = length(states) * dimension))
      testthat::expect_equal(result$max_ix, 7 + length(states) * dimension)
      testthat::expect_null(result$ix)
      testthat::expect_false("MY" %in% names(result))
      testthat::expect_identical(result[[field]][[3L - species]], object)
      y <- seq_len(result$max_ix)
      if (module == "BRQS") {
        testthat::expect_equal(env$F_fqZ.BRQS(0, y, result, species), 4 * y[ix$Bz_ix])
        testthat::expect_equal(env$F_fqM.BRQS(0, y, result, species),
                               4 * (y[ix$Bu_ix] + y[ix$By_ix] + y[ix$Bz_ix]))
        testthat::expect_equal(env$F_eggs.BRQS(0, y, result, species),
                               30 * (y[ix$Qu_ix] + y[ix$Qy_ix] + y[ix$Qz_ix]))
        outputs <- rbind(y, y + 1)
        parsed <- env$parse_MY_orbits.BRQS(outputs, result, species)
        for (state in states) {
          testthat::expect_identical(parsed[[state]], outputs[, ix[[paste0(state, "_ix")]]])
        }
      } else {
        vars <- env$get_XH_vars.SIPd(y, result, species)
        testthat::expect_identical(vars$I, y[ix$I_ix])
        testthat::expect_identical(vars$H, y[ix$S_ix] + y[ix$I_ix] + y[ix$P_ix])
        testthat::expect_equal(env$F_I.SIPd(0, y, result, species), 0.5 * y[ix$I_ix])
        for (ntimes in c(1L, 3L)) {
          outputs <- matrix(seq_len(ntimes * result$max_ix), nrow = ntimes)
          parsed <- env$parse_XH_orbits.SIPd(outputs, result, species)
          testthat::expect_identical(names(parsed), c("S", "I", "P", "H", "ni", "true_pr"))
          for (state in c("S", "I", "P")) {
            testthat::expect_identical(parsed[[state]], outputs[, ix[[paste0(state, "_ix")]]])
          }
          testthat::expect_identical(parsed$H, parsed$S + parsed$I + parsed$P)
          testthat::expect_equal(parsed$ni, 0.5 * parsed$I / parsed$H)
          testthat::expect_equal(parsed$true_pr, parsed$I / parsed$H)
        }
      }
    }
  }
})

testthat::test_that("index consumers no longer refer to legacy storage", {
  for (file in c("adult-MY-BRQS.R", "malaria-XH-SIPd.R", "human-XH-SEISd.R")) {
    code <- paste(deparse(parse(testthat::test_path("..", "..", "R", file))), collapse = "\n")
    testthat::expect_false(grepl("xds_obj$ix$X", code, fixed = TRUE))
    testthat::expect_false(grepl("xds_obj$ix$MY", code, fixed = TRUE))
    testthat::expect_false(grepl("xds_obj$MY[[", code, fixed = TRUE))
    testthat::expect_false(grepl("xds_obj$X_obj[[", code, fixed = TRUE))
  }
})
