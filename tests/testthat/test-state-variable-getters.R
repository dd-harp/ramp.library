testthat::test_that("BRQS and garki getters return uniquely named states and derived quantities", {
  for (module in c("BRQS", "garki")) {
    component <- if (module == "BRQS") "MY" else "XH"
    file <- if (module == "BRQS") "adult-MY-BRQS.R" else "malaria-garki.R"
    getter <- paste0("get_", component, "_vars.", module)
    env <- new.env(parent = baseenv())
    for (definition in parse(testthat::test_path("..", "..", "R", file))) {
      if (is.call(definition) &&
          as.character(definition[[1]])[1] %in% c("<-", "=") &&
          identical(as.character(definition[[2]]), getter)) {
        eval(definition, env)
      }
    }
    states <- if (module == "BRQS") {
      c("Su", "Bu", "Ru", "Qu", "Sy", "By", "Ry", "Qy", "Sz", "Bz", "Rz", "Qz")
    } else c("H", "x2", "x3", "x4", "y1", "y2", "y3")
    expected_names <- if (module == "BRQS") {
      c(states, "B", "R", "Q", "S", "M")
    } else c("H", "x1", states[-1])
    for (species in 1:2) {
      for (dimension in c(1L, 3L)) {
        ix <- setNames(lapply(seq_along(states), function(j) {
          5L + (j - 1L) * dimension + seq_len(dimension)
        }), paste0(states, "_ix"))
        model <- list()
        objects <- list(list(ix = list()), list(ix = list()))
        objects[[species]]$ix <- ix
        model[[paste0(component, "_obj")]] <- objects
        y <- seq_len(5L + length(states) * dimension)
        vars <- env[[getter]](y, model, species)
        testthat::expect_identical(names(vars), expected_names)
        testthat::expect_true(all(lengths(vars) == dimension))
        for (state in states) {
          testthat::expect_identical(vars[[state]], y[ix[[paste0(state, "_ix")]]])
        }
        if (module == "BRQS") {
          for (behavior in c("B", "R", "Q", "S")) {
            expected <- Reduce(`+`, lapply(paste0(behavior, c("u", "y", "z")), function(state) {
              y[ix[[paste0(state, "_ix")]]]
            }))
            testthat::expect_identical(vars[[behavior]], expected)
          }
          testthat::expect_identical(vars$M, Reduce(`+`, vars[states]))
          testthat::expect_identical(with(vars, Bu + By + Bz), vars$B)
        } else {
          testthat::expect_identical(vars$x1, vars$H - vars$x2 - vars$x3 -
                                      vars$x4 - vars$y1 - vars$y2 - vars$y3)
        }
      }
    }
  }
})
