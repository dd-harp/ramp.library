test_that("library XH modules expose independent test positivity", {
  modules <- c("SEIR", "SEIRV", "SIR", "SIRS", "SEIS", "SEISd",
               "SIP", "SIPd", "garki")
  for (module in modules) {
    constructor <- get(paste0("make_XH_obj_", module))
    default <- constructor(2)
    expect_identical(default$q, c(0.8, 0.8), info = module)
    expect_identical(constructor(2, q = 0.6)$q, c(0.6, 0.6), info = module)
    object <- constructor(2, options = list(q = c(0.4, 0.7)))
    expect_identical(object$q, c(0.4, 0.7), info = module)
    expect_error(constructor(2, q = c(0.1, 0.2, 0.3)))
    expect_error(constructor(2, q = "invalid"))
    if (module == "SIP") {
      object$mda <- object$msat <- ramp.xds::F_zero
    }
    model <- list(XH_obj = list(default, object))
    original <- model
    pars <- get_XH_pars(model, 2)
    expect_identical(pars$q, object$q, info = module)
    expect_identical(change_XH_pars(model, list(), 2), model, info = module)
    changed <- change_XH_pars(model, list(q = c(0.3, 0.9)), 2)
    expect_identical(changed$XH_obj[[1]], default, info = module)
    expect_identical(get_XH_pars(changed, 2)$q, c(0.3, 0.9), info = module)
    expect_identical(changed$XH_obj[[2]][setdiff(names(object), "q")],
                     object[setdiff(names(object), "q")], info = module)

    states <- list(I = c(10, 40), H = c(100, 200),
                   y1 = c(1, 2), y2 = c(3, 4), y3 = c(6, 8))
    expected <- object$q * F_prevalence(states, object)
    expect_equal(F_pr_obs(states, object), expected, info = module)
    expect_true(get_skillset_XH(module)$pr_obs, info = module)
    trajectory <- lapply(states, function(x) rbind(x, x * 2, x * 3))
    model$outputs <- list(time = 0:2, orbits = list(XH = list(NULL, trajectory)))
    expected <- object$q * F_prevalence(trajectory, object)
    expect_equal(get_PR(model, "obs", 2),
                 list(time = 0:2, pr = expected, method = "obs"), info = module)
    expect_identical(model$XH_obj, original$XH_obj, info = module)
  }
})
