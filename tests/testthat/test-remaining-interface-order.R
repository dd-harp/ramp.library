test_that("library initial-value methods match the reordered generic defaults", {
  families <- list(
    change_MY_inits=c("BQ", "BQm", "BRQS", "SEI"),
    setup_L_inits="stages",
    change_L_inits="stages"
  )
  for (name in names(families)) {
    for (module in families[[name]]) {
      expect_identical(formals(getS3method(name, module)), formals(get(name)))
    }
  }
})

test_that("library initial-value methods accept options before the species index", {
  model <- xds_setup(MYname="SEI", Lname="stages")
  for (name in c("change_MY_inits", "setup_L_inits", "change_L_inits")) {
    fun <- get(name)
    is_MY <- name == "change_MY_inits"
    component <- if (is_MY) "MY_obj" else "L_obj"
    field <- if (is_MY) "M" else "E"
    options <- setNames(list(20), field)
    model2 <- if (is_MY) setup_MY_obj("SEI", model, s=2) else
      setup_L_obj("stages", 1, model, s=2)
    changed <- fun(model2, options)
    expect_equal(changed, fun(model2, options=options, s=1))
    expect_equal(changed[[component]][[1]]$inits[[field]], 20)
    expect_equal(fun(model2), fun(model2, list(), 1))
    second <- fun(model2, options, 2)
    expect_equal(second[[component]][[1]], model2[[component]][[1]])
    expect_equal(second[[component]][[2]]$inits[[field]], 20)
  }
})
