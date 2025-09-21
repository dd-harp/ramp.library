library(deSolve)
library(ramp.xds)

numeric_tol <- 1e-5

test_that("human SIP_xde model remains at equilibrium", {
  HPop <- c(100, 500, 250)
  nPatches=3
  residence = 1:3

  I <- c(20, 120, 80)
  b <- 0.55
  c <- 0.15
  r <- 1/200
  eta <- c(1/30, 1/40, 1/35)
  rho <- c(0.05, 0.1, 0.15)
  xi <- rep(0, 3)
  Xo = list(b=b, c=c, r=r, eta=eta, rho=rho, xi=xi)
  foi = stats::rnorm(3, 0.5, .1)
  params <- xds_setup_human(Xname ="SIP", XHoptions=Xo, nPatches=nPatches, HPop=HPop, residence = residence)
  params$terms$FoI[[1]] <- foi

  steady_state_X(foi, HPop, params, 1) ->ss
  params <- change_XH_inits(params, 1, ss)


  # set initial conditions
  y0 <- as.vector(unlist(get_inits(params)))

  out <- deSolve::ode(y = y0, times = c(0, 1000), func = function(t, y, pars, s) {
    list(dXHdt(t, y, pars, s))
  }, parms = params, method = 'lsoda', s=1)

  expect_equal(as.vector(out[2L, params$XH_obj[[1]]$ix$I_ix+1]), ss$I, tolerance = numeric_tol)
  expect_equal(as.vector(out[2L, params$XH_obj[[1]]$ix$P_ix+1]), ss$P, tolerance = numeric_tol)
})
