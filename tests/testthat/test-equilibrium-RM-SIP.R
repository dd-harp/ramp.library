library(expm)
library(ramp.xds)
library(MASS)
library(deSolve)

numeric_tol <- 1e-5

test_that("test equilibrium with macdonald adults (DDE), SIP_xde humans, trivial", {


  # set number of patches and strata
  nPatches <- 2
  residence <- c(1:2)
  nStrata <- length(residence)
  membership <- c(1:2)

  # parameters
  b <- 0.55
  c <- 0.15
  r <- 1/200
  eta <- c(1/30, 1/40)
  xi <- c(0, 0)
  rho <- c(0.05, 0.1)
  Xo <- list(b=b, c=c, r=r, eta=eta, xi=xi, rho=rho)

  f <- rep(0.3, nPatches)
  q <- rep(0.9, nPatches)
  g <- rep(1/10, nPatches)
  sigma <- rep(1/100, nPatches)
  mu <- rep(0, nPatches)
  nu <- rep(1/2, nPatches)
  eggsPerBatch <- 30
  eip <- 11

  # mosquito movement calK
  calK <- matrix(0, nPatches, nPatches)
  calK[upper.tri(calK)] <- rexp(sum(1:(nPatches-1)))
  calK[lower.tri(calK)] <- rexp(sum(1:(nPatches-1)))
  calK <- calK/rowSums(calK)
  calK <- t(calK)

  # omega matrix
  Omega <- compute_Omega_xde(g, sigma, mu, calK)
  Upsilon <- expm::expm(-Omega * eip)

  MYZo <- list(nPatches=nPatches,
               f=f, q=q, g=g, sigma=sigma, mu=mu, nu=nu, eggsPerBatch=eggsPerBatch,
               eip=eip, Omega=Omega, Upsilon=Upsilon, calK=calK)


  # human PfPR and H
  foi <- rnorm(2, 1/500, .002)
  H <- rpois(n = nStrata, lambda = 1000)
  residence = c(1,2)
  searchWtsH = c(1,1)
  class(Xo) <- "SIP"
  xde_steady_state_X(foi,H,Xo) -> ssI
  I <- ssI$I
  Px <- ssI$P

  eir <- foi/b
  Xo$I=I
  Xo$P=Px

  # TaR
  dg <- rbeta(nStrata, 80,20)
  TaR <- matrix(c(dg[1], 1-dg, dg[2]), 2, 2)

  # ambient pop
  W <- compute_W(searchWtsH, H, TaR)
  beta <- compute_beta(H, W, searchWtsH, TaR)

  # biting distribution matrix
  fqZ <- eir2fqZ(eir, beta)

  # kappa
  kappa <- t(beta) %*% (I*c)

  # equilibrium solutions for adults
  Z <- fqZ/f/q
  MY <- diag(1/as.vector(f*q*kappa), nPatches, nPatches) %*% solve(Upsilon) %*% Omega %*% Z
  Y <- solve(Omega) %*% (diag(as.vector(f*q*kappa), nPatches, nPatches) %*% MY)
  M <- MY + Y
  P <- solve(diag(f, nPatches) + Omega) %*% diag(f, nPatches) %*% M
  Lambda <- Omega %*% M
  class(MYZo) <- "macdonald"
  xde_steady_state_MYZ(Lambda, kappa, MYZo) -> ss

  MYZo$M=M
  MYZo$P=P
  MYZo$Y=Y
  MYZo$Z=Z


  Lo = list(Lambda=Lambda)

  # parameters for exDE
  params <- xds_setup(MYZname = "macdonald", MYZopts=MYZo, Lname = "trivial", Lopts=Lo, TimeSpent = TaR, calK=calK,
                      Xname = "SIP", Xopts=Xo, HPop=H, membership=membership, nPatches=nPatches, residence=residence)


  params <- xds_solve(params, 730, 1)

  out <- params$outputs$last_y

  M_sim <- out[params$ix$MYZ[[1]]$M_ix]
  P_sim <- out[params$ix$MYZ[[1]]$P_ix]
  Y_sim <- out[params$ix$MYZ[[1]]$Y_ix]
  Z_sim <- out[params$ix$MYZ[[1]]$Z_ix]
  I_sim <- out[params$ix$X[[1]]$I_ix]
  Px_sim <- out[params$ix$X[[1]]$P_ix]


  expect_equal(as.vector(M_sim), as.vector(M), tolerance = numeric_tol)
  expect_equal(as.vector(P_sim), as.vector(P), tolerance = numeric_tol)
  expect_equal(as.vector(Y_sim), as.vector(Y), tolerance = numeric_tol)
  expect_equal(as.vector(Z_sim), as.vector(Z), tolerance = numeric_tol)
  expect_equal(as.vector(I_sim), as.vector(I), tolerance = numeric_tol)
  expect_equal(as.vector(Px_sim), as.vector(Px), tolerance = numeric_tol)
})
