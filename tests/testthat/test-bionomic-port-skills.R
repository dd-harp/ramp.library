test_that("library MY skillsets identify the implemented bionomic ports", {
  parameters <- c("f", "q", "g", "sigma", "sigma_b", "sigma_q",
                  "mu", "nu", "eip", "K_matrix")
  ports <- list(
    SEI = c("f", "q", "g", "sigma", "mu", "nu", "eip", "K_matrix"),
    BQ = c("f", "q", "g", "sigma_b", "sigma_q", "mu", "nu", "eip"),
    BQm = c("f", "q", "g", "sigma_b", "sigma_q", "mu", "nu", "eip"),
    BRQS = character()
  )
  for (module in names(ports)) {
    skills <- get_skillset_MY(module)
    flags <- skills[grepl("_is_port$", names(skills))]
    expect_named(flags, paste0(parameters, "_is_port"))
    expect_identical(unname(unlist(flags)), parameters %in% ports[[module]],
                     info = module)
    calls <- all.names(body(getS3method("MBionomics", module)), functions = TRUE)
    expect_setequal(parameters[parameters %in% ports[[module]]],
                    parameters[paste0("F_", parameters) %in% calls])
    model <- list(MY_obj = list(structure(list(), class = "trivialMY"),
                                structure(list(), class = module)))
    expect_identical(get_skillset_MY(model, 2), skills)
    expect_false(skills$delay)
    expect_false("not_implemented" %in% names(skills))
  }
})

test_that("stages reports no bionomic parameter ports", {
  skills <- get_skillset_L("stages")
  flags <- skills[grepl("_is_port$", names(skills))]
  expect_identical(flags, list(psi_is_port = FALSE, phi_is_port = FALSE,
                               xi_is_port = FALSE, theta_is_port = FALSE))
  calls <- all.names(body(getS3method("LBionomics", "stages")), functions = TRUE)
  expect_false(any(paste0("F_", c("psi", "phi", "xi", "theta")) %in% calls))
  model <- list(L_obj = list(structure(list(), class = "trivialL"),
                             structure(list(), class = "stages")))
  expect_identical(get_skillset_L(model, 2), skills)
  expect_identical(skills[c("dts", "delay")],
                   list(dts = FALSE, delay = FALSE))
  expect_false("not_implemented" %in% names(skills))
})
