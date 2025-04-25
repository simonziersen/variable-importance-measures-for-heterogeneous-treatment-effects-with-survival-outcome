resSurv2 <- function(surv, data, times, cause = 1, diff = TRUE){
  dat <- as.data.table(data)
  dat1 <- copy(dat)[, A := factor(1)]
  dat0 <- copy(dat)[, A := factor(0)]
  evalTimes_tmp <- c(0, sort(surv$y[, "time"][surv$y[, "status"] == cause]))
  evalTimes <- evalTimes_tmp[evalTimes_tmp < times]
  if(diff){
    sHat <- predictRisk(S, type = "survival", newdata = dat1, times = evalTimes[which(evalTimes <= times)]) -
      predictRisk(S, type = "survival", newdata = dat0, times = evalTimes[which(evalTimes <= times)])
  } else{
    sHat <- predictRisk(S, type = "survival", newdata = dat, times = evalTimes[which(evalTimes <= times)])
  }
  tau <- sHat %*% diff(c(evalTimes, times))
  tau
}

resSurv <- function(surv, newtimes, newX, newA){
  evalTimes_tmp <- c(0, sort(surv$time[surv$status == 1]))
  evalTimes <- evalTimes_tmp[evalTimes_tmp < newtimes]
  sHat <- predict(surv, newtimes = evalTimes, newX = newX, newA = newA)
  res <- sHat$surv %*% diff(c(evalTimes, newtimes))
  res
}

