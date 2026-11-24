skillset_modules <- list(
  XH=c("SEIR", "SEIRV", "SEIS", "SEISd", "SIP", "SIPd", "SIR", "SIRS", "garki"),
  MY=c("BQ", "BQm", "BRQS", "SEI"),
  L="stages"
)

test_that("library skill-set methods return lists by name and model", {
  for (component in names(skillset_modules)) {
    getter <- get(paste0("get_skillset_", component))
    for (module in skillset_modules[[component]]) {
      skills <- getter(module)
      expect_type(skills, "list")
      expect_identical(skills$delay, module %in% c("SEISd", "SIPd"))
      model <- list()
      field <- paste0(component, "_obj")
      model[[field]] <- list(structure(list(), class=module))
      original <- model
      expect_identical(getter(model), skills)
      expect_identical(model, original)
      model[[field]] <- c(list(NULL), model[[field]])
      expect_identical(getter(model, 2), skills)
    }
  }
})

test_that("setup preserves skill-set lookup and selects the delay solver", {
  for (component in names(skillset_modules)) {
    # BRQS construction has an existing undefined `theta` parameter.
    for (module in setdiff(skillset_modules[[component]], "BRQS")) {
      args <- list()
      args[[paste0(component, "name")]] <- module
      model <- do.call(xds_setup, args)
      getter <- get(paste0("get_skillset_", component))
      object <- model[[paste0(component, "_obj")]][[1]]
      expect_false("skill_set" %in% names(object))
      expect_identical(getter(model), getter(module))
      expect_identical(as.character(model$xde),
                       if (module %in% c("SEISd", "SIPd")) "dde" else "ode")
    }
  }
})
