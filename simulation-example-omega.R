library(survival)
library(mgcv)
library(riskRegression)

f1_tmp <- list.files("functions", full.names = TRUE)
f1 <- f1_tmp[grepl(".R", f1_tmp)]
lapply(f1, source)



# Cox-Weibull distribution function and its invers
inv_weibull <- function(u, theta, scale, shape){
  t <- (-log(u) / (scale * exp(theta))) ^ (1 / shape)
  return(t)
}  

weibull <- function(t, theta, scale, shape){
  p <- exp(-scale*exp(theta)*t^shape)
  return(p)
}

expit <- function(x) exp(x)/(1+exp(x))

# Set parameters 
times <- 10                      # event horizon
shape <- 2                      # shape parameter for weibull distributed event times 
scale <- 0.0025                    # scale parameter for weibull distributed event times 
n <- 500

  

  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  x4 <- rnorm(n)
  xx <- data.frame(x1,x2,x3,x4)
  A <- rbinom(n, 1, expit(0.3*x1 + 0.3*x2))
  A_fac <- factor(A)
  eventTime <- inv_weibull(runif(n), theta = -x1 - x2 - 0.2*x3 + 0.1*x4 + A*(-2 + 0.5*x1 - 0.3*x2), scale = scale, shape = shape)
  censTime <- inv_weibull(runif(n), theta = -0.3*x1, scale = 0.1*scale, shape = shape)
  obsTime <- pmin(eventTime, censTime)
  delta <- 1*(eventTime <= censTime)

  # prob scale covTim
  # correct cox and glm
  
  correctGAMprob <- covVimSurv(time = obsTime,
                              status = delta,
                              treatment = A_fac,
                              confounders = xx,
                              l = 1,
                              nuisance = list(surv = list(method = "coxph", form = "x1*A+x2*A+x3+x4"),
                                              cens = list(method = "coxph"),
                                              prop = list(method = "glm", interactions = FALSE),
                                              pseudoRegTau = list(method = "gam", interactions = TRUE),
                                              pseudoRegX = list(method = "gam", interactions = TRUE)),
                              evaltimes = times,
                              tauLearner = "S",
                              folds = 1)
  
  correctGAMCFprob <- covVimSurv(time = obsTime,
                              status = delta,
                              treatment = A_fac,
                              confounders = xx,
                              l = 1,
                              nuisance = list(surv = list(method = "coxph", form = "x1*A+x2*A+x3+x4"),
                                              cens = list(method = "coxph"),
                                              prop = list(method = "glm", interactions = FALSE),
                                              pseudoRegTau = list(method = "gam", interactions = TRUE),
                                              pseudoRegX = list(method = "gam", interactions = TRUE)),
                              evaltimes = times,
                              tauLearner = "S",
                              folds = 10)
  
  correctRFprob <- covVimSurv(time = obsTime,
                              status = delta,
                              treatment = A_fac,
                              confounders = xx,
                              l = 1,
                              nuisance = list(surv = list(method = "coxph", form = "x1*A+x2*A+x3+x4"),
                                              cens = list(method = "coxph"),
                                              prop = list(method = "glm", interactions = FALSE),
                                              pseudoRegTau = list(method = "rfsrc"),
                                              pseudoRegX = list(method = "rfsrc")),
                              evaltimes = times,
                              tauLearner = "S",
                              folds = 1)
  
  correctRFCFprob <- covVimSurv(time = obsTime,
                              status = delta,
                              treatment = A_fac,
                              confounders = xx,
                              l = 1,
                              nuisance = list(surv = list(method = "coxph", form = "x1*A+x2*A+x3+x4"),
                                              cens = list(method = "coxph"),
                                              prop = list(method = "glm", interactions = FALSE),
                                              pseudoRegTau = list(method = "rfsrc"),
                                              pseudoRegX = list(method = "rfsrc")),
                              evaltimes = times,
                              tauLearner = "S",
                              folds = 10)
  
  
  # RF
  RFGAMprob <- covVimSurv(time = obsTime,
                               status = delta,
                               treatment = A_fac,
                               confounders = xx,
                               l = 1,
                               nuisance = list(surv = list(method = "rfsrc"),
                                               cens = list(method = "rfsrc"),
                                               prop = list(method = "rfsrc"),
                                               pseudoRegTau = list(method = "gam", interactions = TRUE),
                                               pseudoRegX = list(method = "gam", interactions = TRUE)),
                               evaltimes = times,
                               tauLearner = "S",
                               folds = 1)
  
  RFGAMCFprob <- covVimSurv(time = obsTime,
                                 status = delta,
                                 treatment = A_fac,
                                 confounders = xx,
                                 l = 1,
                                 nuisance = list(surv = list(method = "rfsrc"),
                                                 cens = list(method = "rfsrc"),
                                                 prop = list(method = "rfsrc"),
                                                 pseudoRegTau = list(method = "gam", interactions = TRUE),
                                                 pseudoRegX = list(method = "gam", interactions = TRUE)),
                                 evaltimes = times,
                                 tauLearner = "S",
                                 folds = 10)
  
  RFRFprob <- covVimSurv(time = obsTime,
                         status = delta,
                         treatment = A_fac,
                         confounders = xx,
                         l = 1,
                         nuisance = list(surv = list(method = "rfsrc"),
                                         cens = list(method = "rfsrc"),
                                         prop = list(method = "rfsrc"),
                                         pseudoRegTau = list(method = "rfsrc"),
                                         pseudoRegX = list(method = "rfsrc")),
                         evaltimes = times,
                         tauLearner = "S",
                         folds = 1)
  
  RFRFCFprob <- covVimSurv(time = obsTime,
                                status = delta,
                                treatment = A_fac,
                                confounders = xx,
                                l = 1,
                           nuisance = list(surv = list(method = "rfsrc"),
                                           cens = list(method = "rfsrc"),
                                           prop = list(method = "rfsrc"),
                                           pseudoRegTau = list(method = "rfsrc"),
                                           pseudoRegX = list(method = "rfsrc")),
                                evaltimes = times,
                                tauLearner = "S",
                                folds = 10)
  
  
  rescorrectprob <- rbind(c(unlist(correctGAMprob), method = "correctGAM"),
                          c(unlist(correctGAMCFprob), method = "correctGAMCF"),
                          c(unlist(correctRFprob), method = "correctRF"),
                          c(unlist(correctRFCFprob), method = "correctRFCF"))
  
  
  resrfprob <- rbind(c(unlist(RFRFprob), method = "RFRF"),
                     c(unlist(RFGAMprob), method = "RFGAM"),
                     c(unlist(RFRFCFprob), method = "RFRFCF"),
                     c(unlist(RFGAMCFprob), method = "RFGAMCF"))
  
  resprob <- rbind(rescorrectprob, resrfprob)
  
  
  








