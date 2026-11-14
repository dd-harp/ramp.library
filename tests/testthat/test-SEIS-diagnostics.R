test_that("SEIS modules configure, return, and change SIS-style diagnostic parameters", {
  for (module in c("SEIS", "SEISd")) {
    constructor <- get(paste0("make_XH_obj_", module))
    defaults <- constructor(2)
    expect_equal(defaults$q_lm, c(0.8, 0.8))
    expect_equal(defaults$q_rdt, c(0.8, 0.8))
    expect_equal(defaults$q_pcr, c(0.9, 0.9))

    options <- list(q_lm=c(0.2, 0.4), q_rdt=c(0.5, 0.6), q_pcr=c(0.7, 0.9))
    expect_equal(constructor(2, q_lm=0.3)$q_lm, c(0.3, 0.3))
    expect_equal(constructor(2, options, q_lm=0.3)$q_lm, options$q_lm)

    model <- make_xds_object_template()
    model$nStrata <- 2
    model <- setup_XH_obj(module, c(1, 1), model, options)
    pars <- get_XH_pars(model)
    expect_named(pars, c("b", "c", "r", "nu", "q_lm", "q_rdt", "q_pcr", "rho", "sigma"))
    expect_equal(pars[names(options)], options)

    changed <- change_XH_pars(model, list(q_rdt=c(0.1, 0.3)))
    expect_equal(get_XH_pars(changed)$q_rdt, c(0.1, 0.3))
    expect_equal(get_XH_pars(changed)[c("b", "c", "r", "nu", "q_lm", "q_pcr")],
                 pars[c("b", "c", "r", "nu", "q_lm", "q_pcr")])

    vars <- list(H=c(100, 200), E=c(20, 40), I=c(10, 40))
    expect_equal(F_prevalence(vars, model$XH_obj[[1]]), c(0.1, 0.2))
    for (test in c("lm", "rdt", "pcr")) {
      diagnostic <- get(paste0("F_pfpr_by_", test))
      expect_equal(diagnostic(vars, model$XH_obj[[1]]),
                   options[[paste0("q_", test)]] * vars$I/vars$H)
      expect_null(dim(diagnostic(vars, model$XH_obj[[1]])))
      model$outputs$orbits$XH[[1]] <- vars
      model$outputs$time <- c(0, 1)
      expect_equal(get_PR(model, method=test)$pr,
                   diagnostic(vars, model$XH_obj[[1]]))
      zero_obj <- constructor(1, setNames(list(0), paste0("q_", test)))
      expect_equal(diagnostic(list(H=100, I=10), zero_obj), 0)
    }
  }
})

test_that("diagnostic parameters do not change SEIS infection dynamics", {
  for (module in c("SEIS", "SEISd")) {
    model <- make_xds_object_template()
    model <- setup_XH_obj(module, 1, model, list(H=100, E=20, I=10))
    model$max_ix <- 0
    model <- setup_XH_ix(model, 1)
    model$terms$FoI[[1]] <- 0.1
    y <- unlist(model$XH_obj[[1]]$inits, use.names=FALSE)
    baseline <- dXHdt(0, y, model, 1)
    changed <- change_XH_pars(model, list(q_lm=0, q_rdt=1, q_pcr=0.3))
    expect_equal(dXHdt(0, y, changed, 1), baseline)
    pars <- get_XH_pars(model)
    if (module == "SEIS") {
      expect_equal(baseline, c(0, 7-pars$nu*20, pars$nu*20-pars$r*10))
    } else {
      expect_equal(baseline, c(0, 7, -pars$r*10, 7))
    }
  }
})
