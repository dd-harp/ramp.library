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
nStrata <- 3
H <- c(100, 500, 250)
nPatches <- 3
residence <- 1:3 
params <- make_xds_object_template("ode", "human", nPatches, 1, residence) 

## -----------------------------------------------------------------------------
b <- 0.5
c <- 0.15
r <- 1/200
eta <- c(1/30, 1/40, 1/35)
rho <- c(0.05, 0.1, 0.15)
xi <- rep(0, 3)
sigma <-rep(0.5,3)
Xo = list(b=b,c=c,r=r,eta=eta,rho=rho,sigma=sigma,xi=xi)
params = setup_XH_obj("SIP", params, 1, Xo) 

## -----------------------------------------------------------------------------
eir <- c(1,2,3)/365
steady_state_X(eir*b, H, params) -> ss
params = setup_XH_inits(params, H, 1, ss)

## -----------------------------------------------------------------------------
Xo <- c(Xo, ss) 

## -----------------------------------------------------------------------------
MYo = list(
  Z = eir*H, f=1, q=1
)

## -----------------------------------------------------------------------------
params = setup_MY_obj("trivial", params, 1, MYo)
params = setup_MY_inits(params, 1)
params = setup_L_obj("trivial", params, 1)
params = setup_L_inits(params, 1)

## -----------------------------------------------------------------------------
params = make_indices(params)

## -----------------------------------------------------------------------------
steady_state_X(eir*b, H, params, 1)

## -----------------------------------------------------------------------------
y0 <- as.vector(unlist(get_inits(params)))

## -----------------------------------------------------------------------------
out <- deSolve::ode(y = y0, times = c(0, 730), xde_derivatives, parms= params, method = 'lsoda') 
get_XH_vars(out[2,-1], params, 1)

## ----out.width = "100%"-------------------------------------------------------
colnames(out)[params$XH_obj[[1]]$ix$H_ix+1] <- paste0('H_', 1:params$nStrata)
colnames(out)[params$XH_obj[[1]]$ix$I_ix+1] <- paste0('I_', 1:params$nStrata)
colnames(out)[params$XH_obj[[1]]$ix$P_ix+1] <- paste0('P_', 1:params$nStrata)

out <- as.data.table(out)
out <- melt(out, id.vars = 'time')
out[, c("Component", "Strata") := tstrsplit(variable, '_', fixed = TRUE)]
out[, variable := NULL]

ggplot(data = out, mapping = aes(x = time, y = value, color = Strata)) +
  geom_line() +
  facet_wrap(. ~ Component, scales = 'free') +
  theme_bw()

## -----------------------------------------------------------------------------
xds_setup_human(Xname="SIP", nPatches=3, residence = 1:3, HPop=H, XHoptions= Xo, MYoptions = MYo) -> test_SIP_xde

## -----------------------------------------------------------------------------
steady_state_X(b*eir, H, test_SIP_xde, 1) -> out1
out1 <- unlist(out1)

## -----------------------------------------------------------------------------
xds_solve(test_SIP_xde, 365, 365)$outputs$last_y -> out2
approx_equal(out2,out1) 

