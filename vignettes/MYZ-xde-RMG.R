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
model <- xds_setup(MYZname="RMG", Lname="trivial", Xname = "trivial",  residence=residence, HPop =HPop, nPatches=3)

## -----------------------------------------------------------------------------
model <- xds_solve(model, dt=5)

## ----fig.height=7.5, fig.width=5.5, eval=F------------------------------------
#  par(mfrow = c(2,1))
#  xds_plot_M(model)
#  xds_plot_YZ(model, add_axes = F)
#  xds_plot_YZ_fracs(model)

## ----fig.height=4, fig.width=5.5, eval=F--------------------------------------
#  xds_plot_YZ_fracs(model)

