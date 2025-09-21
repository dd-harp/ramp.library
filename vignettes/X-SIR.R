## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(ramp.xds)
library(ramp.library)
library(deSolve)
library(data.table)
library(ggplot2)

## ----echo=FALSE---------------------------------------------------------------
#devtools::load_all()

## -----------------------------------------------------------------------------
test_SIR <- xds_setup(MYname ="macdonald", Xname="SIR")

## -----------------------------------------------------------------------------
xds_solve(test_SIR, 365*10)-> test_SIR
unlist(get_XH_vars(test_SIR$outputs$last_y, test_SIR, 1)) -> out1

## -----------------------------------------------------------------------------
foi_eq = tail(test_SIR$outputs$orbits$XH[[1]]$foi, 1)
steady_state_X(1/365, 1000, test_SIR, 1) -> out2

## -----------------------------------------------------------------------------
sum(abs(out2 - out1))


