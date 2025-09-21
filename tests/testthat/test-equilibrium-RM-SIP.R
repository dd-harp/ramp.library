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

  # mosquito movement K_matrix
  K_matrix <- matrix(0, nPatches, nPatches)
  K_matrix[upper.tri(K_matrix)] <- rexp(sum(1:(nPatches-1)))
  K_matrix[lower.tri(K_matrix)] <- rexp(sum(1:(nPatches-1)))
  K_matrix <- K_matrix/rowSums(K_matrix)
  K_matrix <- t(K_matrix)

  # omega matrix
  Omega <- make_Omega_xde(g, sigma, mu, K_matrix)
  Upsilon <- expm::expm(-Omega * eip)

  MYo <- list(nPatches=nPatches,
               f=f, q=q, g=g, sigma=sigma, mu=mu, nu=nu, eggsPerBatch=eggsPerBatch,
               eip=eip, K_matrix=K_matrix)


  # human PfPR and H
  foi <- rnorm(2, 1/500, .002)
  H <- rpois(n = nStrata, lambda = 1000)
  residence = c(1,2)
  searchWtsH = c(1,1)

  dg <- rbeta(nStrata, 80,20)
  TaR <- matrix(c(dg[1], 1-dg, dg[2]), 2, 2)

  params <- xds_setup(MYname = "macdonald", MYoptions = MYo,
                      Lname = "trivial", XHoptions = Xo,
                      Xname = "SIP",
                      TimeSpent = TaR, K_matrix=K_matrix,
                      HPop=H, membership=membership,
                      nPatches=nPatches, residence=residence)
  params <- check_MY(params, 1)
  steady_state_X(foi, H, params, 1) -> ssI
  params <- change_XH_inits(params, 1, ssI)

  I <- ssI$I
  Px <- ssI$P

  eir <- foi/b

  # ambient pop
  W <- F_W_available(searchWtsH, H, TaR)
  beta <- F_beta(H, W, searchWtsH, TaR)

  # biting distribution matrix
  fqZ <- solve(beta) %*% eir

  # kappa
  kappa <- t(beta) %*% (I*c)

  # equilibrium solutions for adults
  Z <- fqZ/f/q
  fqk <- as.vector(f*q*kappa)
  MY <- solve(Upsilon) %*% Omega %*% Z
  Y <- solve(Omega) %*% MY
  M <- diag(1/fqk)%*%(diag(fqk)%*%Y + Omega %*%Y)
  P <- solve(diag(f, nPatches) + Omega) %*% diag(f, nPatches) %*% M
  Lambda <- Omega %*% M

  Lo = list(Lambda=Lambda)
  params <- change_L_pars(params, 1, Lo)

  steady_state_MY(Lambda, kappa, params, 1) -> ss
  params <- change_MY_inits(params, 1, ss)

  params <- xds_solve(params, 730, 1)

  out <- params$outputs$last_y

  M_sim <- get_MY_vars(out, params, 1)$M
  P_sim <- get_MY_vars(out, params, 1)$P
  Y_sim <- get_MY_vars(out, params, 1)$Y
  Z_sim <- get_MY_vars(out, params, 1)$Z
  I_sim <- get_XH_vars(out, params, 1)$I
  Px_sim <- get_XH_vars(out, params, 1)$P


  expect_equal(as.vector(M_sim), as.vector(M), tolerance = numeric_tol)
  expect_equal(as.vector(P_sim), as.vector(P), tolerance = numeric_tol)
  expect_equal(as.vector(Y_sim), as.vector(Y), tolerance = numeric_tol)
  expect_equal(as.vector(Z_sim), as.vector(Z), tolerance = numeric_tol)
  expect_equal(as.vector(I_sim), as.vector(I), tolerance = numeric_tol)
  expect_equal(as.vector(Px_sim), as.vector(Px), tolerance = numeric_tol)
})
