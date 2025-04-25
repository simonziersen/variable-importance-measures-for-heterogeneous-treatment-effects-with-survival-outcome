Htu2 <- function(surv, data, times, cause = 1){
  
  dat <- as.data.table(data)
  evalTimes_tmp <- c(0, sort(surv$y[, "time"][surv$y[, "status"] == cause]))
  evalTimes <- evalTimes_tmp[evalTimes_tmp < times]
  sHat <- predictRisk(S, type = "survival", newdata = dat, times = evalTimes[which(evalTimes <= times)])
  tmp <- t(t(sHat) * diff(c(evalTimes, times)))
  res <- rowCumSum(tmp[, dim(tmp)[2]:1, drop = FALSE])[, dim(tmp)[2]:1]
  res
}
Htu <- function(surv, times, predtimes){
  evalTimes_tmp <- sort(times)
  evalTimes <- evalTimes_tmp[evalTimes_tmp < predtimes]
  sHat <- surv[, which(evalTimes <= predtimes)]
  tmp <- t(t(sHat) * diff(c(evalTimes, predtimes)))
  res <- rowCumSum(tmp[, dim(tmp)[2]:1])[, dim(tmp)[2]:1]
  res
}
