test_that("library XH setup accepts options before the species index", {
  model <- make_xds_object_template()
  options <- list(r=0.02)

  positional <- setup_XH_obj("SIP", 1, model, options)
  named <- setup_XH_obj("SIP", 1, model, options=options, i=1)
  expect_equal(positional, named)
  expect_equal(positional$XH_obj[[1]]$r, 0.02)
  expect_equal(
    setup_XH_obj("SIP", 1, model),
    setup_XH_obj("SIP", 1, model, list(), 1)
  )
})
