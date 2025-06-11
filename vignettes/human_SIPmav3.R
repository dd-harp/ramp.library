## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(ramp.xds)
library(viridisLite)
#library(ramp.library)
library(deSolve)
devtools::load_all()

## ----echo=FALSE---------------------------------------------------------------
mav3_opts = list()
mav3_opts$mod_b = "im" 
mav3_opts$b_opts = list(b0 = .55, Nb=10) 
mav3_opts$mod_severe_h = "protect" 
mav3_opts$severe_h_opts = list(p0 = .55, Np=1) 
mav3_opts$mod_severe_x = "protect" 
mav3_opts$severe_x_opts = list(rate0 = 1/300, Np=1) 
mav3_opts$mod_moderate_h = "protect" 
mav3_opts$moderate_h_opts = list(p0 = 1, Np=20) 
mav3_opts$mod_moderate_x = "protect" 
mav3_opts$moderate_x_opts = list(rate0 = 1/300, Np=20) 
mav3_opts$mod_c = "dim" 
mav3_opts$c_opts = list(c0=.8, pda=80, Nc=4) 

mod <- xds_setup_cohort(Xname = "SIPmav3", Xopts = mav3_opts, age_par = makepar_F_type2())
mod <- setup_exposure_nb(mod, .1)

## -----------------------------------------------------------------------------

mod$EIRpar$eir = 1/365
mod1 <- xds_solve_cohort(mod, A=10, da=1)
XH1 <- get_XH(mod1)

mod$EIRpar$eir = 5/365
mod5 <- xds_solve_cohort(mod, A=10, da=1)
XH5 <- get_XH(mod5)

mod$EIRpar$eir = 25/365
mod25 <- xds_solve_cohort(mod, A=10, da=1)
XH25 <- get_XH(mod25)

## ----fig.height=5, fig.width=8------------------------------------------------
omega = make_function(makepar_F_type2())
aa = seq(0, 7300, by=10)
plot(aa/365, omega(aa), type = "l", ylim = c(0, 1.5), xlab = "Age (Years)", ylab = expression(omega(a)))
segments(0,1, 20, 1, lty =2)

## ----fig.height=5, fig.width=8------------------------------------------------
clrs <- turbo(7)
age <- XH1$time

b_a_1 = with(mod$Xpar[[1]]$b_mod, F_infect_im(XH1$vh, b0, Nb))
b_a_3 = with(mod$Xpar[[1]]$b_mod, F_infect_im(XH5$vh, b0, Nb))
b_a_10 = with(mod$Xpar[[1]]$b_mod, F_infect_im(XH25$vh, b0, Nb))
plot(age/365, b_a_1, type="l", xlab = "Age (Years)", ylab = expression(F[b](a)), ylim = c(0,1), col = clrs[2])
lines(age/365, b_a_3, col = clrs[4])
lines(age/365, b_a_10, col = clrs[6])
segments(0, 0.55, 20, 0.55, lty =2)

## ----fig.height=5, fig.width=8------------------------------------------------
b_a_1 = with(mod$Xpar[[1]]$b_mod, F_infect_im(XH1$vh, b0, Nb))
b_a_3 = with(mod$Xpar[[1]]$b_mod, F_infect_im(XH5$vh, b0, Nb))
b_a_10 = with(mod$Xpar[[1]]$b_mod, F_infect_im(XH25$vh, b0, Nb))
w_a <- omega(age)
plot(age/365, 0.55*w_a, type="l", xlab = "Age (Years)", ylab=expression(list(omega(a), omega(a)*b(a))), lwd=2)
lines(age/365, b_a_1*w_a, col = clrs[2])
lines(age/365, b_a_3*w_a, col = clrs[4])
lines(age/365, b_a_10*w_a, col = clrs[6])

## ----fig.height=5, fig.width=8------------------------------------------------
a1 = age[-1]
foi1 <- diff(XH1$vh)
foi3 <- diff(XH5$vh)
foi10 <- diff(XH25$vh)
plot(a1/365, foi10, type="l", xlab = "Age (Years)", ylab=expression(list(omega(a), omega(a)*b(a))), col = clrs[6], ylim = c(0, 0.014))
lines(a1/365, foi1, col = clrs[2])
lines(a1/365, foi3, col = clrs[4])

## ----fig.height=5, fig.width=8------------------------------------------------
plot(a1/365, foi10/foi1, type="l", xlab = "Age (Years)", ylab=expression(list(omega(a), omega(a)*b(a))), col = clrs[7], ylim = c(0,25))
lines(a1/365, foi10/foi3, type="l", col = clrs[5])
lines(a1/365, foi3/foi1, type="l", col = clrs[3])
segments(0,1,10,1, col = grey(0.5))
#segments(0,3,10,3, col = grey(0.5))
#segments(0,2,10,2, col = grey(0.5))

## ----fig.height=5, fig.width=8------------------------------------------------
severe1 = diff(XH1$severe)
severe5 = diff(XH5$severe)
severe25 = diff(XH25$severe)
plot(XH1$time[-1]/365, diff(XH25$severe)*365, type = "l", col = clrs[6], ylab = "Severe Malaria, per Year", xlab = "Age (in Years)")
lines(XH1$time[-1]/365, diff(XH5$severe)*365, col = clrs[4])
lines(XH1$time[-1]/365, diff(XH1$severe)*365, col = clrs[2]) 


## ----fig.height=5, fig.width=8------------------------------------------------
plot(XH1$time[-1]/365, diff(XH25$moderate)*365, type = "l", col = clrs[6], ylab = "moderate Malaria, per Year", xlab = "Age (in Years)", ylim = c(0,3.7))
lines(XH1$time[-1]/365, diff(XH5$moderate)*365, col = clrs[4])
lines(XH1$time[-1]/365, diff(XH1$moderate)*365, col = clrs[2]) 

## ----fig.height=5, fig.width=8------------------------------------------------
plot(XH1$time[-1]/365, diff(XH25$mild)*365, type = "l", col = clrs[6], ylab = "mild Malaria, per Year", xlab = "Age (in Years)", ylim = c(0,5))
lines(XH1$time[-1]/365, diff(XH5$mild)*365, col = clrs[4])
lines(XH1$time[-1]/365, diff(XH1$mild)*365, col = clrs[2]) 

## ----fig.height=5, fig.width=8------------------------------------------------
plot(XH1$time[-1]/365, diff(XH25$treat)*365, type = "l", col = clrs[6], ylab = "Yearly Treatment Rate", xlab = "Age (in Years)", ylim = c(0, 2000))
lines(XH1$time[-1]/365, diff(XH5$treat)*365, col = clrs[4])
lines(XH1$time[-1]/365, diff(XH1$treat)*365, col = clrs[2]) 

