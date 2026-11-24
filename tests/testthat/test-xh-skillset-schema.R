test_that("library XH skill sets report diagnostics only for malaria models", {
  modules <- c("SEIR", "SEIRV", "SIR", "SIRS", "SEIS", "SEISd",
               "SIP", "SIPd", "garki")
  methods <- new.env(parent = baseenv())
  for (file in list.files(test_path("..", "..", "R"), "[.]R$", full.names = TRUE)) {
    for (definition in parse(file)) {
      if (is.call(definition) && length(definition) >= 3L &&
          is.symbol(definition[[2]]) &&
          grepl("^get_skillset_XH[.]", as.character(definition[[2]]))) {
        eval(definition, methods)
      }
    }
  }
  fields <- c("dts", "xde", "delay", "demography", "malaria",
              "pr_obs", "pf_lm", "pf_rdt", "pf_pcr", "mda", "msat")
  for (module in modules) {
    getter <- get(paste0("get_skillset_XH.", module), methods)
    skills <- getter(NULL, 1)
    malaria <- !module %in% c("SIR", "SIRS", "SEIR", "SEIRV")
    expected_fields <- if (malaria) fields else
      setdiff(fields, c("pf_lm", "pf_rdt", "pf_pcr"))
    expect_named(skills, expected_fields)
    expect_true(all(vapply(skills, function(x) {
      is.logical(x) && length(x) == 1L && !is.na(x)
    }, logical(1))))
    expect_false(skills$dts)
    expect_true(skills$xde)
    expect_identical(skills$delay, module %in% c("SEISd", "SIPd"))
    expect_true(skills$demography)
    expect_true(skills$pr_obs)
    expect_identical(skills$malaria,
                     !module %in% c("SIR", "SIRS", "SEIR", "SEIRV"))
    for (field in c("pf_lm", "pf_rdt", "pf_pcr")) {
      if (malaria) {
        expect_identical(skills[[field]],
                         module %in% c("SEIS", "SEISd", "SIP", "SIPd"),
                         info = paste(module, field))
      } else {
        expect_false(field %in% names(skills), info = paste(module, field))
      }
    }
    for (field in c("mda", "msat")) {
      expect_identical(skills[[field]],
                       module %in% c("SEIS", "SEISd", "SIP", "garki"),
                       info = paste(module, field))
    }
  }
})
