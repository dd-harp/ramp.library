## ----message=FALSE, warning=FALSE---------------------------------------------
library(knitr)
library(ramp.xds)
library(deSolve) 
library(ramp.library)

## ----echo=F-------------------------------------------------------------------
#devtools::load_all()

## -----------------------------------------------------------------------------
HPop = rep(1000, 3)
residence = c(1:3) 
model <- xds_setup(MYname="RMG", Lname="trivial", Xname = "trivial",  residence=residence, HPop =HPop, nPatches=3)

## -----------------------------------------------------------------------------
model <- xds_solve(model)

## ----fig.height=4, fig.width=5.5----------------------------------------------
xds_plot_M(model)
xds_plot_Y(model, add = T)
xds_plot_Z(model, add = T)

## ----fig.height=4, fig.width=5.5----------------------------------------------
xds_plot_Y_fracs(model)
xds_plot_Z_fracs(model, add=T)

