test_that("library MY setup accepts options before the species index", {
  model <- xds_setup()
  options <- list(f=0.4)

  positional <- setup_MY_obj("BQ", model, options)
  named <- setup_MY_obj("BQ", model, options=options, s=1)
  expect_equal(positional, named)
  expect_equal(positional$MY_obj[[1]]$f, 0.4)
  expect_equal(setup_MY_obj("BQ", model), setup_MY_obj("BQ", model, list(), 1))

  second <- setup_MY_obj("BQ", model, options, 2)
  expect_equal(second$MY_obj[[1]], model$MY_obj[[1]])
  expect_equal(second$MY_obj[[2]]$f, 0.4)
})
