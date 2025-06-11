## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----fig.height=5, fig.width=8------------------------------------------------
aoi = seq(20, 400, by = 5)
test_pos = function(aoi, mx=.95, rt=1/100){
  mx*exp(-rt*(aoi-20)) 
}
plot(aoi, test_pos(aoi), type = "l", 
     xlab = expression(list(alpha, " Age of Infection")), 
     ylab = expression(list(xi(alpha), " Test Positivity")), 
     ylim = c(0,1), xlim = c(0,400))
segments(20, 0, 20, 1, lty = 2)
lines(aoi, test_pos(aoi, .8, 1/90))
lines(aoi, test_pos(aoi, .7, 1/80))

## -----------------------------------------------------------------------------
pr_by_aoi = function(aoi, m, mx=0.95, rt=1/100){
  xi = test_pos(aoi, mx, rt)  
  1 - exp(-xi*m)
}

## ----fig.height=5, fig.width=8------------------------------------------------
plot(aoi, pr_by_aoi(aoi, 1, .8, 1/90), type = "l", xlab = expression(list(bar(alpha), " Mean AoI")), ylab = expression(list(x(bar (alpha)), " Apparent PR")), ylim = c(0,1), xlim = c(0,400))
segments(20, 0, 20, 1, lty = 2)
lines(aoi, pr_by_aoi(aoi, 2,  .8, 1/90))
lines(aoi, pr_by_aoi(aoi, 4, .8, 1/90))

## ----fig.height=5, fig.width=8------------------------------------------------
foi_age = function(a, foi=1, shft=30, A=1.5, B=2){
  aa = (a+shft)/365
  foi*A*aa/(B+aa)
}

age = seq(0, 10*365, by=365/12)
plot(age/365, foi_age(age), type = "l", ylim = c(0,1.3), xlab = "Age (in Years)", ylab = "FoI by Age")
segments(0, 1, 10, 1, col = grey(0.5))

## ----fig.height=5, fig.width=8------------------------------------------------
cum_exposure = function(a, foi=1/365, shft=30, A=1.5, B=2){
  fa <- function(a, foi, shft, A, B) 
    integrate(foi_age, 0, a, foi=foi, shft=shft, A=A, B=B)$value
  sapply(a, fa, foi=foi, shft=shft, A=A,B=B)
}


plot(age/365, cum_exposure(age, foi = 10/365), type = "l", 
     xlab = "Age (in Years)", ylab = "Cumulative Exposure")
lines(age/365, cum_exposure(age, foi = 3/365)) 
lines(age/365, cum_exposure(age, foi = 1/365)) 

## -----------------------------------------------------------------------------
incidence = function(a, foi=3/365, z=1/10, N=100, shft=30, A=1.5, B=2){
  w <- sapply(a, cum_exposure, foi=foi, shft=shft, A=A, B=B)
  foi_age(a, foi)*z*exp(-w*z) 
}

## ----fig.height=5, fig.width=8------------------------------------------------
base = incidence(age, 30/365)
plot(age/365, incidence(age, 30/365), type = "l", xlab = "Age (in Years)", ylab = "Incidence", ylim = range(0, base))
lines(age/365, incidence(age, 10/365), col = "violet")
lines(age/365, incidence(age, 3/365), col = "darkred")
lines(age/365, incidence(age, 1/365), col = "darkblue")
lines(age/365, incidence(age, 1/365/3), col = "darkgreen")


## ----fig.height=5, fig.width=8------------------------------------------------
cum_incidence = function(a, foi=3/365, z=1/10, N=100, shft=30, A=1.5, B=2){
  ci <- function(a, foi, z, N, shft, A, B)
    integrate(incidence, 0, a, foi=foi, z=z, shft=shft, A=A, B=B)$value
  sapply(a, ci, foi=foi, z=z, shft=shft, A=A, B=B)
}

base = cum_incidence(age, 10/365)

plot(age/365, cum_incidence(age, 30/365), type = "l", xlab = "Age (in Years)", ylab = "Incidence", ylim = range(0, base))
lines(age/365, cum_incidence(age, 10/365), col = "violet")
lines(age/365, cum_incidence(age, 3/365), col = "darkred")
lines(age/365, cum_incidence(age, 1/365), col = "darkblue")
lines(age/365, cum_incidence(age, 1/365/3), col = "darkgreen")


## ----message=FALSE, warning=FALSE---------------------------------------------
library(ramp.xds)
library(ramp.library)
library(deSolve)
library(data.table)
library(ggplot2)

## ----echo=FALSE---------------------------------------------------------------
devtools::load_all()

