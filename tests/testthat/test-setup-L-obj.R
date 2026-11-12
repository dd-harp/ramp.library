test_that("library L setup accepts options before the species index", {
  model <- make_xds_object_template()
  options <- list(psi_E=0.2)

  positional <- setup_L_obj("stages", 1, model, options)
  named <- setup_L_obj("stages", 1, model, options=options, s=1)
  expect_equal(positional, named)
  expect_equal(positional$L_obj[[1]]$psi_E, 0.2)
  expect_equal(
    setup_L_obj("stages", 1, model),
    setup_L_obj("stages", 1, model, list(), 1)
  )

  second <- setup_L_obj("stages", 1, positional, options, 2)
  expect_equal(second$L_obj[[1]], positional$L_obj[[1]])
  expect_equal(second$L_obj[[2]]$psi_E, 0.2)
})
