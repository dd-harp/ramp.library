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
test_SIR <- xds_setup(MYZname ="RM", Xname="SIR")

## -----------------------------------------------------------------------------
xds_solve(test_SIR, 365*10)-> test_SIR
unlist(list_Xvars(test_SIR$outputs$last_y, test_SIR, 1)) -> out1

## -----------------------------------------------------------------------------
foi_eq = test_SIR$Xpar[[1]]$b*tail(test_SIR$outputs$terms$EIR,1)
xde_steady_state_X(1/365, 1000, test_SIR$Xpar[[1]]) -> out2

## -----------------------------------------------------------------------------
out2
out1

