## -----------------------------------------------------------------------------
suppressMessages(library(knitr))
suppressMessages(library(ramp.xds))
#suppressMessages(library(ramp.library))
library(deSolve)
library(ramp.library)
#devtools::load_all()

## -----------------------------------------------------------------------------
#devtools::load_all()

## ----eval=F-------------------------------------------------------------------
#  wh <- xds_setup(dlay = 'dde',
#                  Xname = "workhorse",
#                  HPop=1000,
#                  Lname = "trivial",
#                  Lopts = list(
#                    Lambda = 3000,
#                    season = function(t){1+sin(2*pi*t/365)}
#                    )
#                  )

## ----eval=F-------------------------------------------------------------------
#  wh <- xde_solve(wh, Tmax = 5*365, dt=10)

## ----eval=F-------------------------------------------------------------------
#  wh1 = wh
#  wh1$Xpar[[1]]$zeta_1 = 0.02
#  wh1 <- xde_solve(wh1, Tmax = 5*365, dt=10)

## ----eval=F-------------------------------------------------------------------
#  xds_plot_X(wh)
#  xds_plot_X(wh1, llty=2, add=FALSE)

## ----eval=F-------------------------------------------------------------------
#  time = wh$outputs$time
#  with(wh$outputs$orbits$XH[[1]],{
#    A = A0+A1+A2+A3+A4
#    I = I1+I2+I3+I4
#    plot(time, U + P + G, type = "l", ylim = range(0, H), col = "darkblue")
#    lines(time, A+I, col = "darkred")
#    lines(time, A, col = "purple")
#    lines(time, I-I1, col = "orange")
#    lines(time, P+G, col = "darkgreen")
#  })
#  
#  with(wh1$outputs$orbits$XH[[1]],{
#    A = A0+A1+A2+A3+A4
#    I = I1+I2+I3+I4
#    lines(time, U + P + G, col = "darkblue", lty=2)
#    lines(time, A+I, col = "darkred", lty=2)
#    lines(time, A, col = "purple", lty=2)
#    lines(time, I-I1, col = "orange", lty=2)
#    lines(time, P+G, col = "darkgreen", lty=2)
#  })

