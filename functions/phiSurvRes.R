phiSurvRes <- function(surv, cens, treat, newtimes, newstatus, newX, newA, predtimes){
  require(data.table)
  #browser()
  nobs <- dim(newX)[1]
  
  # predict nuisance at jump times of hat(Lambda)
  time.jump1 <- unique(sort(surv$time[surv$status == 1]))
  Shats <- predict(surv, newtimes = c(0, time.jump1), newX = newX, newA = newA)
  Ghats <- predict(cens, newtimes = c(0, time.jump1), newX = newX, newA = newA)
  pmhat <- predict(treat, newX = newX)
  
  newA1 <- rep(1, length(newA))
  newA0 <- rep(0, length(newA))
  Shats1 <- predict(surv, newtimes = c(0, time.jump1), newX = newX, newA = newA1)
  Shats0 <- predict(surv, newtimes = c(0, time.jump1), newX = newX, newA = newA0)
  
  
  Ht <- Htu(Shats$surv, c(0, time.jump1), predtimes)
  tau1 <- Htu(Shats1$surv, c(0, time.jump1), predtimes)[,1]
  tau0 <- Htu(Shats0$surv, c(0, time.jump1), predtimes)[,1]
  
  # martingales
  
  dN <- do.call(cbind, lapply(predtimes, function(t0){
    newtimes <= t0 & newstatus == 1
  }))
  
  #dL <- cbind(Shats$cumhaz[,1], t(diff(t(Shats$cumhaz))))
  dL <- cbind(t(diff(t(Shats$cumhaz))))
  
  mgEvalIndex <- do.call(cbind,lapply(predtimes, function(tau){
    prodlim::sindex(jump.times = time.jump1, eval.times = pmin(newtimes,tau))
  }))
  
  ### martingale integrand
  integrand <- Ht/(Shats$surv[,1:(1 + max(mgEvalIndex)), drop = FALSE] * Ghats$surv[,1:(1 + max(mgEvalIndex)), drop = FALSE])
  
  ### dLambda integral
  L_tmp <- rowCumSum(integrand * dL[,1:(1 + max(mgEvalIndex))])
  L <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L_tmp[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  ### dN integral
  N <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- integrand[cbind(1:nobs, mgEvalIndex[, x])] * dN[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  ### martingale
  M <- N - L
  
  ## pseudo outcome; phi
  # survival as outcome
  est1_tmp <- tau1 - (1*(newA == 1) / pmhat) * M
  est0_tmp <- tau0 - (1*(newA == 0) / (1 - pmhat)) * M
  
  
  # 1-S as outcome
  # est1_tmp <- 1 - Stau1 + (1*(data[[treatvar]] == 1) / pmhat) * M * StauObs
  # est0_tmp <- 1 - Stau0 + 1*(data[[treatvar]] == 0) / (1 - pmhat) * M * StauObs
  
  ## ate estimates
  est1 <- (1/nobs) * colSums(est1_tmp)
  est0 <- (1/nobs) * colSums(est0_tmp)
  ate <- (1/nobs) * colSums(est1_tmp - est0_tmp)
  
  res <- cbind(est1, est0, ate)
  rownames(res) <- predtimes
  
  phi1 <- est1_tmp
  colnames(phi1) <- predtimes
  phi0 <- est0_tmp
  colnames(phi0) <- predtimes
  
  return(list(phi = phi1 - phi0, phi1 = phi1, phi0 = phi0))
}