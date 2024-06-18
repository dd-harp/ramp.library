## -----------------------------------------------------------------------------
EIPo = list(EIPname = "static_xde", eip=12) 
model <- xde_setup(MYZname="RMG_xde", Lname="trace", Xname = "trace", EIPopts = EIPo, nPatches=3)

## -----------------------------------------------------------------------------
model <- xde_solve(model, dt=5)

## ----fig.height=7.5, fig.width=5.5--------------------------------------------
par(mfrow = c(2,1))
xds_plot_M(model)
xds_plot_YZ(model, add_axes = F)
xds_plot_YZ_fracs(model)

## ----fig.height=4, fig.width=5.5----------------------------------------------
xds_plot_YZ_fracs(model)

