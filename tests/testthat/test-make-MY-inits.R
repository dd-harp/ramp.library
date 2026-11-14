test_that("BQ initial-value constructors use variable names without trailing zeros", {
  variables <- c("Bu", "Qu", "By", "Qy", "Bz", "Qz")
  values <- setNames(as.list(seq_along(variables) + 10), variables)
  expected <- lapply(values, rep, times=2)
  defaults <- list(Bu=c(5, 5), Qu=c(1, 1), By=c(1, 1),
                   Qy=c(1, 1), Bz=c(1, 1), Qz=c(1, 1))

  for (module in c("BQ", "BQm")) {
    constructor <- get(paste0("make_MY_inits_", module))
    expect_identical(names(formals(constructor)), c("nPatches", "options", variables))
    expect_equal(constructor(2), defaults)
    expect_equal(do.call(constructor, c(list(nPatches=2), values)), expected)
    expect_equal(constructor(2, values), expected)
    expect_equal(constructor(2, list(Bu=30), Bu=20)$Bu, c(30, 30))

    model <- xds_setup(MYname=module, MYoptions=values)
    expect_equal(model$MY_obj[[1]]$inits, values)
  }
})
