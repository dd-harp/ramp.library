## ----fig.width=6, fig.height=3.5----------------------------------------------
garki <- xds_setup(Xname = "garki", MYZname = "si")
garki <- xds_solve(garki, Tmax=3650, dt = 15)
xds_plot_X(garki, 1)

