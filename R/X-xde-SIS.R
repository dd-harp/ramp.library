#' Split a stratum into two strata, assigning a fraction `p` a new biting weight that is multiplied by a factor `fac` compared with the old one. The biting weight for the remaining `1-p` gets a new factor `1/fac`
#'
#' @inheritParams ramp.work::split_stratum_by_biting
#'
#' @return pars a list
#' @param i the host species index
#' @param j the stratum to split
#' @param p the fraction that gets multiplied by `fac`
#' @param fac a factor
#' @return a [list]
#' @export
split_stratum_by_biting.SIS = function(pars, i, j, p, fac){
  stopifnot(i <= pars$Hpar[[i]]$nStrata)
  nj_ix = pars$Hpar[[i]]$nStrata + 1
  pars$Hpar[[i]]$nStrata = nj_ix


  Xpar = pars$Xpar[[1]]
  Xpar$b <- c(Xpar$b, Xpar$b[j])
  Xpar$r <- c(Xpar$r, Xpar$r[j])
  Xpar$c <- c(Xpar$c, Xpar$c[j])
  pars$Xpar[[1]] = Xpar

  Xinits = pars$Xinits[[1]]
  Xinits$S <- c(Xinits$S, p*Xinits$S[j])
  Xinits$S[j] <- Xinits$S[j]*(1-p)
  Xinits$I <- c(Xinits$I, p*Xinits$I[j])
  Xinits$I[j] <- Xinits$I[j]*(1-p)
  pars$Xinits[[1]] = Xinits

  residence = pars$BFpar$residence[[i]]
  pars$BFpar$residence[[i]] = c(residence, residence[j])
  for(s in 1:pars$nVectors){
    wts_f = pars$BFpar$searchWts[[i]][[s]]
    pars$BFpar$searchWts[[i]][[s]] = c(wts_f, fac*wts_f[j])
  }

  H =  pars$Hpar[[i]]$H
  pars$Hpar[[i]]$H = c(H, p*H[j])
  pars$Hpar[[i]]$H[j] = H[j]*(1-p)

  TimeSpent = pars$BFpar$TimeSpent[[i]]
  TimeSpent = cbind(TimeSpent, TimeSpent[,j])

  for(s in 1:pars$nVectors) pars = make_TaR(0, pars, i, 1)

  pars <- ramp.xde::make_indices(pars)

  return(pars)
}
