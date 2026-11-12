test_that("library MY initial-value interfaces match the generic defaults", {
  for (module in c("BQ", "BQm", "BRQS", "SEI")) {
    expect_identical(formals(getS3method("setup_MY_inits", module)),
                     formals(setup_MY_inits))
  }
})

test_that("library MY initial values accept positional, named and default arguments", {
  model <- xds_setup(MYname="BQ")
  options <- list(Bu0=20)
  positional <- setup_MY_inits(model, options)
  expect_equal(positional, setup_MY_inits(model, options=options, s=1))
  expect_equal(positional$MY_obj[[1]]$inits$Bu, 20)
  expect_equal(setup_MY_inits(model), setup_MY_inits(model, list(), 1))

  model <- setup_MY_obj("BQ", model, s=2)
  second <- setup_MY_inits(model, options, 2)
  expect_equal(second$MY_obj[[1]], model$MY_obj[[1]])
  expect_equal(second$MY_obj[[2]]$inits$Bu, 20)
})
