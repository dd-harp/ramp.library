test_that("SEIS modules configure and change treatment parameters without changing other parameters", {
  for (module in c("SEIS", "SEISd")) {
    constructor <- get(paste0("make_XH_obj_", module))
    expect_equal(constructor(2)$rho, c(0, 0))
    expect_equal(constructor(2)$sigma, c(0, 0))
    expect_equal(constructor(2, rho=0.3, sigma=0.2)$rho, c(0.3, 0.3))
    expect_equal(constructor(2, list(sigma=0.4), sigma=0.2)$sigma, c(0.4, 0.4))

    model <- make_xds_object_template()
    model$nStrata <- 2
    model <- setup_XH_obj(module, c(1, 1), model,
                          list(rho=c(0.2, 0.5), sigma=c(0.1, 0.3)))
    pars <- get_XH_pars(model)
    expect_equal(pars$rho, c(0.2, 0.5))
    expect_equal(pars$sigma, c(0.1, 0.3))
    changed <- change_XH_pars(model, list(rho=c(0.4, 0.8), sigma=c(0.2, 0.6)))
    expect_equal(get_XH_pars(changed)$rho, c(0.4, 0.8))
    expect_equal(get_XH_pars(changed)$sigma, c(0.2, 0.6))
    expect_equal(get_XH_pars(changed)[setdiff(names(pars), c("rho", "sigma"))],
                 pars[setdiff(names(pars), c("rho", "sigma"))])
    expect_true(model$XH_obj[[1]]$skill_set$mda)
    expect_true(model$XH_obj[[1]]$skill_set$msat)
    expect_equal(model$XH_obj[[1]]$mda(0), 0)
    expect_equal(model$XH_obj[[1]]$msat(0), 0)
  }
})

test_that("treatment clears I and diverts incident infections without altering E", {
  for (module in c("SEIS", "SEISd")) {
    model <- make_xds_object_template()
    model <- setup_XH_obj(module, 1, model, list(H=100, E=20, I=10))
    model$max_ix <- 0
    model <- setup_XH_ix(model, 1)
    model$terms$FoI[[1]] <- 0.1
    y <- unlist(model$XH_obj[[1]]$inits, use.names=FALSE)
    baseline <- dXHdt(0, y, model, 1)
    changed <- change_XH_pars(model, list(rho=0.4, sigma=0.2))
    changed <- setup_mass_treatment(changed, mda=function(t) 0.03*t,
                                    msat=function(t) 0.05)
    actual <- dXHdt(2, y, changed, 1)
    incident <- if (module == "SEIS") model$XH_obj[[1]]$nu*20 else 0
    expect_equal(actual[1:2], baseline[1:2])
    expect_equal(actual[3], baseline[3] - 0.4*incident - (0.2+0.06+0.05)*10)
    if (module == "SEISd") expect_equal(actual[4], baseline[4])
    fully_treated <- change_XH_pars(model, list(rho=1))
    expect_equal(dXHdt(0, y, fully_treated, 1)[3], -model$XH_obj[[1]]$r*10)
  }
})

test_that("SEISd applies rho to delayed incidence rather than current exposures", {
  for (rho in c(0, 0.4, 1)) {
    model <- make_xds_object_template()
    model <- setup_XH_obj("SEISd", 1, model,
                          list(H=100, E=0, I=0, nu=2, r=0.1, rho=rho, sigma=0.2))
    model <- setup_mass_treatment(model, mda=function(t) 0.05, msat=function(t) 0.05)
    model$max_ix <- 0
    model <- setup_XH_ix(model, 1)
    derivative <- function(t, y, parms) {
      vars <- get_XH_vars(y, parms, 1)
      parms$terms$FoI[[1]] <- 2/vars$S
      list(dXHdt(t, y, parms, 1))
    }
    output <- deSolve::dede(
      y=unlist(model$XH_obj[[1]]$inits, use.names=FALSE),
      times=seq(0, 6, by=0.1), func=derivative, parms=model,
      rtol=1e-8, atol=1e-8, hmax=0.1
    )
    duration <- pmax(output[, 1]-2, 0)
    expect_equal(unname(output[, 4]), (1-rho)*2/0.4*(1-exp(-0.4*duration)),
                 tolerance=1e-5)
    expect_equal(unname(output[, 3]), 2*pmin(output[, 1], 2), tolerance=1e-5)
  }
})

test_that("SEIS equilibrium helpers include treatment parameters", {
  for (module in c("SEIS", "SEISd")) {
    model <- make_xds_object_template()
    model <- setup_XH_obj(module, 1, model, list(r=0.1, rho=0.4, sigma=0.2))
    pars <- get_XH_pars(model)
    state <- steady_state_X(0.05, 100, model)
    expect_equal(unname(state["S"] + state["E"] + state["I"]), 100)
    expect_equal(unname(0.6*0.05*state["S"] - 0.3*state["I"]), 0)
    if (module == "SEIS") {
      expect_equal(unname(pars$nu*state["E"]), unname(0.05*state["S"]))
    } else {
      expect_equal(unname(state["E"]), unname(0.05*state["S"]*pars$nu))
    }
    fully_treated <- change_XH_pars(model, list(rho=1))
    expect_equal(unname(steady_state_X(0.05, 100, fully_treated)["I"]), 0)
  }
})
