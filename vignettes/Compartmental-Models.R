## ----message=FALSE, warning=FALSE---------------------------------------------
library(knitr)
library(ramp.xds)
library(ramp.library)
library(deSolve)
library(ramp.library)

## ----echo=F-------------------------------------------------------------------
#devtools::load_all()

## -----------------------------------------------------------------------------
library(ramp.xds)
library(ramp.library)
library(deSolve)
library(data.table)
library(ggplot2)

## -----------------------------------------------------------------------------
test_SIR<- xds_setup(MYZname="macdonald", Xname="SIR")
xds_solve(test_SIR, 365*5) -> test_SIR
foi_eq = test_SIR$Xpar[[1]]$b*tail(test_SIR$outputs$terms$EIR,1)
unlist(list_Xvars(test_SIR$outputs$last_y, test_SIR, 1)) -> inf 
xde_steady_state_X(foi_eq, 1000, test_SIR$Xpar[[1]]) -> ss
xds_plot_X(test_SIR)

## -----------------------------------------------------------------------------
test_SIRS<- xds_setup(MYZname="macdonald", Xname="SIRS")
xds_solve(test_SIRS, 365*3) -> test_SIRS
foi_eq = test_SIRS$Xpar[[1]]$b*tail(test_SIRS$outputs$terms$EIR,1)
unlist(list_Xvars(test_SIRS$outputs$last_y, test_SIRS, 1)) -> inf_SIRS
xde_steady_state_X(foi_eq, 1000, test_SIRS$Xpar[[1]]) -> ss_SIRS
xds_plot_X(test_SIRS)

## -----------------------------------------------------------------------------
test_SEIR<- xds_setup(MYZname="macdonald", Xname="SEIR")
xds_solve(test_SEIR, 365*5) -> test_SEIR
foi_eq = test_SEIR$Xpar[[1]]$b*tail(test_SEIR$outputs$terms$EIR,1)
unlist(list_Xvars(test_SEIR$outputs$last_y, test_SEIR, 1)) -> inf_SEIR
xde_steady_state_X(foi_eq, 1000, test_SEIR$Xpar[[1]]) -> ss_SEIR
xds_plot_X(test_SEIR)

## -----------------------------------------------------------------------------
test_SEIRV<- xds_setup(MYZname="macdonald", Xname="SEIRV")
xds_solve(test_SEIRV, 365*5) -> test_SEIRV
foi_eq = test_SEIRV$Xpar[[1]]$b*tail(test_SEIRV$outputs$terms$EIR,1)
unlist(list_Xvars(test_SEIRV$outputs$last_y, test_SEIRV, 1)) -> inf_SEIRV
xde_steady_state_X(foi_eq, 1000, test_SEIRV$Xpar[[1]]) -> ss_SEIRV
xds_plot_X(test_SEIRV)

