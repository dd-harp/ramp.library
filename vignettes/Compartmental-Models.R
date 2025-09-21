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
#devtools::load_all()

## -----------------------------------------------------------------------------
test_SIR<- xds_setup(MYname="macdonald", Xname="SIR")
xds_solve(test_SIR, 365*5) -> test_SIR
foi_eq = tail(test_SIR$outputs$orbits$XH[[1]]$foi,1)
unlist(get_XH_vars(test_SIR$outputs$last_y, test_SIR, 1)) -> inf 
steady_state_X(foi_eq, 1000, test_SIR, 1) -> ss
xds_plot_X(test_SIR)

## -----------------------------------------------------------------------------
test_SIRS<- xds_setup(MYname="macdonald", Xname="SIRS")
xds_solve(test_SIRS, 365*3) -> test_SIRS
foi_eq = tail(test_SIRS$outputs$orbits$XH[[1]]$foi,1)
unlist(get_XH_vars(test_SIRS$outputs$last_y, test_SIRS, 1)) -> inf_SIRS
steady_state_X(foi_eq, 1000, test_SIRS, 1) -> ss_SIRS
xds_plot_X(test_SIRS)

## -----------------------------------------------------------------------------
test_SEIR<- xds_setup(MYname="macdonald", Xname="SEIR")
xds_solve(test_SEIR, 365*5) -> test_SEIR
foi_eq = tail(test_SEIR$outputs$orbits$XH[[1]]$foi,1)
unlist(get_XH_vars(test_SEIR$outputs$last_y, test_SEIR, 1)) -> inf_SEIR
steady_state_X(foi_eq, 1000, test_SEIR, 1) -> ss_SEIR
xds_plot_X(test_SEIR)

## -----------------------------------------------------------------------------
test_SEIRV<- xds_setup(MYname="macdonald", Xname="SEIRV")
xds_solve(test_SEIRV, 365*5) -> test_SEIRV
foi_eq = tail(test_SEIRV$outputs$orbits$XH[[1]]$foi,1)
unlist(get_XH_vars(test_SEIRV$outputs$last_y, test_SEIRV, 1)) -> inf_SEIRV
steady_state_X(foi_eq, 1000, test_SEIRV, 1) -> ss_SEIRV
xds_plot_X(test_SEIRV)

