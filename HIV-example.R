library(speff2trial)
library(riskRegression)
library(data.table)
library(randomForestSRC)
f1_tmp <- list.files("functions", full.names = TRUE)
f1 <- f1_tmp[grepl(".R", f1_tmp)]
lapply(f1, source)
data("ACTG175")
armstmp <- "02"
if(armstmp == "13"){
  dat <- data.table(ACTG175)[arms %in% c(1, 3)]
  dat[,A := as.factor(ifelse(arms == 1, 1, 0))]
} else if(armstmp == "02"){
  dat <- data.table(ACTG175)[arms %in% c(0, 2)]
  dat[,A := as.factor(ifelse(arms == 2, 1, 0))]
}

dat[,pidnum := NULL]
dat[,zprior := NULL]
dat[, strat := as.factor(strat)]
dat[, treat := NULL]
dat[, cd420 := NULL]
dat[, cd496 := NULL]
dat[, cd820 := NULL]
dat[, arms := NULL]
dat[, r := NULL]
dat[, strat := NULL]
dat[, oprior := NULL]
dat[, z30 := NULL]
dat[, preanti := NULL]
dat[, offtrt := NULL]

dat[, hemo := as.factor(hemo)]
dat[, homo := as.factor(homo)]
dat[, drugs := as.factor(drugs)]
dat[, race := as.factor(race)]
dat[, gender := as.factor(gender)]
dat[, str2 := as.factor(str2)]
dat[, symptom := as.factor(symptom)]
dat[karnof == 70, karnof := 80]

dat[, karnof := as.numeric(karnof)]
dat[, age := as.numeric(age)]
dat[, cd40 := as.numeric(cd40)]
dat[, cd80 := as.numeric(cd80)]

RFCF1 <- teVimSurvRes_tmp(time = dat$days,
                     status = dat$cens,
                     treatment = dat$A,
                     confounders = dat[, 1:12],
                     s = 1,
                     nuisance = list(surv = list(method = "rfsrc", interactions = FALSE),
                                     cens = list(method = "rfsrc"),
                                     prop = list(method = "rfsrc"),
                                     pseudoReg = list(method = "rfsrc", interactions = FALSE)),
                     evaltimes = 1000,
                     tauLearner = "S",
                     folds = 10)
RFCF2 <- teVimSurvRes_tmp(time = dat$days,
                     status = dat$cens,
                     treatment = dat$A,
                     confounders = dat[, 1:12],
                     s = 2,
                     nuisance = list(surv = list(method = "rfsrc", interactions = FALSE),
                                     cens = list(method = "rfsrc"),
                                     prop = list(method = "rfsrc"),
                                     pseudoReg = list(method = "rfsrc", interactions = FALSE)),
                     evaltimes = 1000,
                     tauLearner = "S",
                     folds = 10)
RFCF3 <- teVimSurvRes_tmp(time = dat$days,
                     status = dat$cens,
                     treatment = dat$A,
                     confounders = dat[, 1:12],
                     s = 3,
                     nuisance = list(surv = list(method = "rfsrc", interactions = FALSE),
                                     cens = list(method = "rfsrc"),
                                     prop = list(method = "rfsrc"),
                                     pseudoReg = list(method = "rfsrc", interactions = FALSE)),
                     evaltimes = 1000,
                     tauLearner = "S",
                     folds = 10)
RFCF4 <- teVimSurvRes_tmp(time = dat$days,
                     status = dat$cens,
                     treatment = dat$A,
                     confounders = dat[, 1:12],
                     s = 4,
                     nuisance = list(surv = list(method = "rfsrc", interactions = FALSE),
                                     cens = list(method = "rfsrc"),
                                     prop = list(method = "rfsrc"),
                                     pseudoReg = list(method = "rfsrc", interactions = FALSE)),
                     evaltimes = 1000,
                     tauLearner = "S",
                     folds = 10)
RFCF5 <- teVimSurvRes_tmp(time = dat$days,
                     status = dat$cens,
                     treatment = dat$A,
                     confounders = dat[, 1:12],
                     s = 5,
                     nuisance = list(surv = list(method = "rfsrc", interactions = FALSE),
                                     cens = list(method = "rfsrc"),
                                     prop = list(method = "rfsrc"),
                                     pseudoReg = list(method = "rfsrc", interactions = FALSE)),
                     evaltimes = 1000,
                     tauLearner = "S",
                     folds = 10)

RFCF6 <- teVimSurvRes_tmp(time = dat$days,
              status = dat$cens,
              treatment = dat$A,
              confounders = dat[, 1:12],
              s = 6,
              nuisance = list(surv = list(method = "rfsrc"),
                              cens = list(method = "rfsrc"),
                              prop = list(method = "rfsrc"),
                              pseudoReg = list(method = "rfsrc")),
              evaltimes = 1000,
              tauLearner = "S",
              folds = 10)
RFCF7 <- teVimSurvRes_tmp(time = dat$days,
              status = dat$cens,
              treatment = dat$A,
              confounders = dat[, 1:12],
              s = 7,
              nuisance = list(surv = list(method = "rfsrc"),
                              cens = list(method = "rfsrc"),
                              prop = list(method = "rfsrc"),
                              pseudoReg = list(method = "rfsrc")),
              evaltimes = 1000,
              tauLearner = "S",
              folds = 10)
RFCF8 <- teVimSurvRes_tmp(time = dat$days,
              status = dat$cens,
              treatment = dat$A,
              confounders = dat[, 1:12],
              s = 8,
              nuisance = list(surv = list(method = "rfsrc"),
                              cens = list(method = "rfsrc"),
                              prop = list(method = "rfsrc"),
                              pseudoReg = list(method = "rfsrc")),
              evaltimes = 1000,
              tauLearner = "S",
              folds = 10)
RFCF9 <- teVimSurvRes_tmp(time = dat$days,
              status = dat$cens,
              treatment = dat$A,
              confounders = dat[, 1:12],
              s = 9,
              nuisance = list(surv = list(method = "rfsrc"),
                              cens = list(method = "rfsrc"),
                              prop = list(method = "rfsrc"),
                              pseudoReg = list(method = "rfsrc")),
              evaltimes = 1000,
              tauLearner = "S",
              folds = 10)
RFCF10 <- teVimSurvRes_tmp(time = dat$days,
              status = dat$cens,
              treatment = dat$A,
              confounders = dat[, 1:12],
              s = 10,
              nuisance = list(surv = list(method = "rfsrc"),
                              cens = list(method = "rfsrc"),
                              prop = list(method = "rfsrc"),
                              pseudoReg = list(method = "rfsrc")),
              evaltimes = 1000,
              tauLearner = "S",
              folds = 10)
RFCF11 <- teVimSurvRes_tmp(time = dat$days,
              status = dat$cens,
              treatment = dat$A,
              confounders = dat[, 1:12],
              s = 11,
              nuisance = list(surv = list(method = "rfsrc"),
                              cens = list(method = "rfsrc"),
                              prop = list(method = "rfsrc"),
                              pseudoReg = list(method = "rfsrc")),
              evaltimes = 1000,
              tauLearner = "S",
              folds = 10)
RFCF12 <- teVimSurvRes_tmp(time = dat$days,
              status = dat$cens,
              treatment = dat$A,
              confounders = dat[, 1:12],
              s = 12,
              nuisance = list(surv = list(method = "rfsrc"),
                              cens = list(method = "rfsrc"),
                              prop = list(method = "rfsrc"),
                              pseudoReg = list(method = "rfsrc")),
              evaltimes = 1000,
              tauLearner = "S",
              folds = 10)

teVimRes <- cbind(covariate = names(dat[,1:12]),rbind(unlist(RFCF1),
                                          unlist(RFCF2),
                                          unlist(RFCF3),
                                          unlist(RFCF4),
                                          unlist(RFCF5),
                                          unlist(RFCF6),
                                          unlist(RFCF7),
                                          unlist(RFCF8),
                                          unlist(RFCF9),
                                          unlist(RFCF10),
                                          unlist(RFCF11),
                                          unlist(RFCF12)),
                  estimand = "teVim")

fwrite(teVimRes, file = paste0("analysis-results/teVimRes-", armstmp, ".txt"))



# covVIm

cov1 <- covVimSurvRes(time = dat$days,
                   status = dat$cens,
                   treatment = dat$A,
                   confounders = dat[, 1:12],
                   s = 1,
                   nuisance = list(surv = list(method = "rfsrc"),
                                   cens = list(method = "rfsrc"),
                                   prop = list(method = "rfsrc"),
                                   pseudoReg = list(method = "rfsrc")),
                   tauLearner = "S",
                   evaltimes = 1000,
                   folds = 10)

cov1 <- covVimSurvRes(time = dat$days,
                   status = dat$cens,
                   treatment = dat$A,
                   confounders = dat[, 1:12],
                   s = 1,
                   nuisance = list(surv = list(method = "coxph", form = "age + gender"),
                                   cens = list(method = "coxph", form = "age + gender"),
                                   prop = list(method = "rfsrc"),
                                   pseudoReg = list(method = "rfsrc")),
                   tauLearner = "S",
                   evaltimes = 1000,
                   folds = 2)

cov2 <- covVimSurvRes(time = dat$days,
                   status = dat$cens,
                   treatment = dat$A,
                   confounders = dat[, 1:12],
                   s = 2,
                   nuisance = list(surv = list(method = "rfsrc"),
                                   cens = list(method = "rfsrc"),
                                   prop = list(method = "rfsrc"),
                                   pseudoReg = list(method = "rfsrc")),
                   tauLearner = "S",
                   evaltimes = 1000,
                   folds = 10)

cov3 <- covVimSurvRes(time = dat$days,
                   status = dat$cens,
                   treatment = dat$A,
                   confounders = dat[, 1:12],
                   s = 3,
                   nuisance = list(surv = list(method = "rfsrc"),
                                   cens = list(method = "rfsrc"),
                                   prop = list(method = "rfsrc"),
                                   pseudoReg = list(method = "rfsrc")),
                   tauLearner = "S",
                   evaltimes = 1000,
                   folds = 10)

cov4 <- covVimSurvRes(time = dat$days,
                   status = dat$cens,
                   treatment = dat$A,
                   confounders = dat[, 1:12],
                   s = 4,
                   nuisance = list(surv = list(method = "rfsrc"),
                                   cens = list(method = "rfsrc"),
                                   prop = list(method = "rfsrc"),
                                   pseudoReg = list(method = "rfsrc")),
                   tauLearner = "S",
                   evaltimes = 1000,
                   folds = 10)

cov5 <- covVimSurvRes(time = dat$days,
                   status = dat$cens,
                   treatment = dat$A,
                   confounders = dat[, 1:12],
                   s = 5,
                   nuisance = list(surv = list(method = "rfsrc"),
                                   cens = list(method = "rfsrc"),
                                   prop = list(method = "rfsrc"),
                                   pseudoReg = list(method = "rfsrc")),
                   tauLearner = "S",
                   evaltimes = 1000,
                   folds = 10)

cov6 <- covVimSurvRes(time = dat$days,
                   status = dat$cens,
                   treatment = dat$A,
                   confounders = dat[, 1:12],
                   s = 6,
                   nuisance = list(surv = list(method = "rfsrc"),
                                   cens = list(method = "rfsrc"),
                                   prop = list(method = "rfsrc"),
                                   pseudoReg = list(method = "rfsrc")),
                   tauLearner = "S",
                   evaltimes = 1000,
                   folds = 10)

cov7 <- covVimSurvRes(time = dat$days,
                   status = dat$cens,
                   treatment = dat$A,
                   confounders = dat[, 1:12],
                   s = 7,
                   nuisance = list(surv = list(method = "rfsrc"),
                                   cens = list(method = "rfsrc"),
                                   prop = list(method = "rfsrc"),
                                   pseudoReg = list(method = "rfsrc")),
                   tauLearner = "S",
                   evaltimes = 1000,
                   folds = 10)
cov8 <- covVimSurvRes(time = dat$days,
                   status = dat$cens,
                   treatment = dat$A,
                   confounders = dat[, 1:12],
                   s = 8,
                   nuisance = list(surv = list(method = "rfsrc"),
                                   cens = list(method = "rfsrc"),
                                   prop = list(method = "rfsrc"),
                                   pseudoReg = list(method = "rfsrc")),
                   tauLearner = "S",
                   evaltimes = 1000,
                   folds = 10)
cov9 <- covVimSurvRes(time = dat$days,
                   status = dat$cens,
                   treatment = dat$A,
                   confounders = dat[, 1:12],
                   s = 9,
                   nuisance = list(surv = list(method = "rfsrc"),
                                   cens = list(method = "rfsrc"),
                                   prop = list(method = "rfsrc"),
                                   pseudoReg = list(method = "rfsrc")),
                   tauLearner = "S",
                   evaltimes = 1000,
                   folds = 10)
cov10 <- covVimSurvRes(time = dat$days,
                   status = dat$cens,
                   treatment = dat$A,
                   confounders = dat[, 1:12],
                   s = 10,
                   nuisance = list(surv = list(method = "rfsrc"),
                                   cens = list(method = "rfsrc"),
                                   prop = list(method = "rfsrc"),
                                   pseudoReg = list(method = "rfsrc")),
                   tauLearner = "S",
                   evaltimes = 1000,
                   folds = 10)

cov11 <- covVimSurvRes(time = dat$days,
                   status = dat$cens,
                   treatment = dat$A,
                   confounders = dat[, 1:12],
                   s = 11,
                   nuisance = list(surv = list(method = "rfsrc"),
                                   cens = list(method = "rfsrc"),
                                   prop = list(method = "rfsrc"),
                                   pseudoReg = list(method = "rfsrc")),
                   tauLearner = "S",
                   evaltimes = 1000,
                   folds = 10)

cov12 <- covVimSurvRes(time = dat$days,
                   status = dat$cens,
                   treatment = dat$A,
                   confounders = dat[, 1:12],
                   s = 12,
                   nuisance = list(surv = list(method = "rfsrc"),
                                   cens = list(method = "rfsrc"),
                                   prop = list(method = "rfsrc"),
                                   pseudoReg = list(method = "rfsrc")),
                   tauLearner = "S",
                   evaltimes = 1000,
                   folds = 10)

covVimRes <- cbind(names(dat[,1:12]),rbind(unlist(cov1),
                                          unlist(cov2),
                                          unlist(cov3),
                                          unlist(cov4),
                                          unlist(cov5),
                                          unlist(cov6),
                                          unlist(cov7),
                                          unlist(cov8),
                                          unlist(cov9),
                                          unlist(cov10),
                                          unlist(cov11),
                                          unlist(cov12)),
                  estimand = "teVim")

fwrite(covVimRes, file = paste0("analysis-results/covVimRes-", armstmp, ".txt"))