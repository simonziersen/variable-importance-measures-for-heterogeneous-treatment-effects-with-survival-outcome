teVimSurv <- function(time,
                      status,
                      treatment,
                      confounders,
                      s,
                      evaltimes,
                      nuisance = list(),
                      tauLearner = "S",
                      dr = FALSE,
                      drs = FALSE,
                      folds = 1){
  #browser()
  A <- treatment
  treatLabs <- levels(treatment)
  treatment <- 1*(treatment == treatLabs[2])
  confounders <- as.data.frame(confounders)
  N <- length(time)
  sminus <- setdiff(1:dim(confounders)[2], s)
  k <- sample(1:folds, length(status), replace = TRUE)
  if(folds == 1){
    pi <- do.call(propML, c(nuisance$prop, list(X = confounders, A = A)))
    if(tauLearner == "S"){
      S <- do.call(survML, c(nuisance$surv, list(time = time, status = status, X = confounders, A = treatment)))
      C <- do.call(survML, c(nuisance$cens, list(time = time, status = 1*(status==0), X = confounders, A = treatment)))
      phi_tmp <- phiSurv(surv = S,
                            cens = C,
                            treat = pi,
                            newtimes = time,
                            newstatus = status,
                            newX = confounders,
                            newA = treatment,
                            predtimes = evaltimes)
      phi <- phi_tmp$phi
      phi1 <- phi_tmp$phi1
      phi0 <- phi_tmp$phi0
      tau <- predict(S, newX = confounders, newtimes = evaltimes, newA = 1)$surv -
        predict(S, newX = confounders, newtimes = evaltimes, newA = 0)$surv
      #taup <- mean(tau)
      taup <- mean(phi)
      taus_tmp <- do.call(pseudoReg, c(nuisance$pseudoReg, list(pseudo = tau, X = confounders[, sminus, drop = FALSE])))
      taus <- predict(taus_tmp, newX = confounders[, sminus, drop = FALSE])
      varTauPlug <- var(tau)
      varTausPlug <- var(taus)
      thetasPlug <- varTauPlug - varTausPlug
      psisPlug <- thetasPlug/varTauPlug
    } else if(tauLearner == "T"){
      S1 <- do.call(survML, c(nuisance$surv, list(time = time[treatment == 1], status = status[treatment == 1], X = confounders[treatment == 1, ])))
      S0 <- do.call(survML, c(nuisance$surv, list(time = time[treatment == 0], status = status[treatment == 0], X = confounders[treatment == 0, ])))
      C1 <- do.call(survML, c(nuisance$cens, list(time = time[treatment == 1], status = 1*(status[treatment == 1]==0), X = confounders[treatment == 1, ])))
      C0 <- do.call(survML, c(nuisance$cens, list(time = time[treatment == 0], status = 1*(status[treatment == 0]==0), X = confounders[treatment == 0, ])))
      phi1 <- phiSurv(surv = S1,
                     cens = C1,
                     treat = pi,
                     newtimes = time,
                     newstatus = status,
                     newX = confounders,
                     newA = treatment,
                     predtimes = evaltimes)$phi1
      phi0 <- phiSurv(surv = S0,
                     cens = C0,
                     treat = pi,
                     newtimes = time,
                     newstatus = status,
                     newX = confounders,
                     newA = treatment,
                     predtimes = evaltimes)$phi0
      phi <- phi1 - phi0
      tau <- predict(S1, newX = confounders, newtimes = evaltimes)$surv -
        predict(S0, newX = confounders, newtimes = evaltimes)$surv
      #taup <- mean(tau)
      taup <- mean(phi)
      taus_tmp <- do.call(pseudoReg, c(nuisance$pseudoReg, list(pseudo = tau, X = confounders[, sminus, drop = FALSE])))
      taus <- predict(taus_tmp, newX = confounders[, sminus, drop = FALSE])
      varTauPlug <- c(var(tau))
      varTausPlug <- c(var(taus))
      thetasPlug <- varTauPlug - varTausPlug
      psisPlug <- thetasPlug/varTauPlug
    } else if(tauLearner == "DR"){
      S <- do.call(survML, c(nuisance$surv, list(time = time, status = status, X = confounders, A = treatment)))
      C <- do.call(survML, c(nuisance$cens, list(time = time, status = 1*(status==0), X = confounders, A = treatment)))
      phi <- phiSurv(surv = S,
                     cens = C,
                     treat = pi,
                     newtimes = time,
                     newstatus = status,
                     newX = confounders,
                     newA = treatment,
                     predtimes = evaltimes)$phi
      tau_tmps <- do.call(drSurv, list(time = time,
                                       status = status,
                                       treatment = treatment,
                                       confounders = confounders,
                                       evaltimes = evaltimes,
                                       nuisance = nuisance,
                                       s = 0))
      tau_tmp <- predict(tau_tmps, newX = confounders)
      tau <- tau_tmp$tau
      taup <- tau_tmps$taup
      taus_tmp <- do.call(pseudoReg, c(nuisance$pseudoReg, list(pseudo = tau, X = confounders[, sminus, drop = FALSE])))
      taus <- predict(taus_tmp, newX = confounders[, sminus, drop = FALSE])
      varTauPlug <- var(tau)
      varTausPlug <- var(taus)
      thetasPlug <- varTauPlug - varTausPlug
      psisPlug <- thetasPlug/varTauPlug
    }
    
    #taup <- mean(phi)
    
    eifVarTaus_uncen <- (taup - phi)^2 - (taus - phi)^2
    eifVarTau_uncen <- (taup - phi)^2 - (tau - phi)^2
    eifThetas_uncen <- (taus - phi)^2 - (tau - phi)^2
    eifTaup_uncen <- phi
    
    varTau <- mean(eifVarTau_uncen)
    varTaus <- mean(eifVarTaus_uncen)
    thetas <- mean(eifThetas_uncen)
    psis <- thetas/varTau
    taup <- mean(eifTaup_uncen)
    
    eifVarTau <- eifVarTau_uncen - varTau
    eifVarTaus <- eifVarTaus_uncen - varTaus
    eifThetas <- eifThetas_uncen - thetas
    eifPsis <- (eifThetas - psis*eifVarTau)/varTau
    eifTaup <- eifTaup_uncen - taup
    
    var_varTau <- (1/N)*mean(eifVarTau^2)
    var_varTaus <- (1/N)*mean(eifVarTaus^2)
    var_thetas <- (1/N)*mean(eifThetas^2)
    var_psis <- (1/N)*mean(eifPsis^2)
    var_taup <- (1/N)*mean(eifTaup^2)
    
    thetasLog <- log(thetasPlug) + (thetas - thetasPlug)/thetasPlug
    varTauLog <- log(varTauPlug) + (varTau - varTauPlug)/varTauPlug
    psisLogit <- logit(psisPlug) + (varTau/varTauPlug) * ((psis - psisPlug)/(psisPlug - psisPlug^2))
    var_thetasLog <- var_thetas/exp(thetasLog)^2
    var_varTauLog <- var_varTau/exp(varTauLog)^2
    var_psisLogit <- var_psis/(expit(psisLogit) - expit(psisLogit)^2)^2
    
  } else{
    vals <- lapply(1:folds, function(i){
      pi <- do.call(propML, c(nuisance$prop, list(X = confounders[k != i,], A = A[k != i])))
      if(tauLearner == "S"){
        S <- do.call(survML, c(nuisance$surv, list(time = time[k != i], status = status[k != i], X = confounders[k != i,], A = treatment[k != i])))
        C <- do.call(survML, c(nuisance$cens, list(time = time[k != i], status = 1*(status[k != i] == 0), X = confounders[k != i,], A = treatment[k != i])))
        phi <- phiSurv(surv = S,
                       cens = C,
                       treat = pi,
                       newtimes = time,
                       newstatus = status,
                       newX = confounders,
                       newA = treatment,
                       predtimes = evaltimes)$phi
        tau_tmp <- predict(S, newX = confounders, newtimes = evaltimes, newA = 1)$surv -
          predict(S, newX = confounders, newtimes = evaltimes, newA = 0)$surv
        #taup <- mean(tau_tmp[k!=i])
        taup <- mean(phi[k!=i])
        tau <- tau_tmp[k == i]
        taus_tmp <- do.call(pseudoReg, c(nuisance$pseudoReg, list(pseudo = tau_tmp[k != i,], X = confounders[k != i, sminus])))
        taus <- predict(taus_tmp, newX = confounders[k == i, sminus])
      } else if(tauLearner == "T"){
        S1 <- do.call(survML, c(nuisance$surv, list(time = time[k != i & treatment == 1], status = status[k != i & treatment == 1], X = confounders[k != i & treatment == 1,])))
        S0 <- do.call(survML, c(nuisance$surv, list(time = time[k != i & treatment == 0], status = status[k != i & treatment == 0], X = confounders[k != i & treatment == 0,])))
        C1 <- do.call(survML, c(nuisance$cens, list(time = time[k != i & treatment == 1], status = 1*(status[k != i & treatment == 1] == 0), X = confounders[k != i & treatment == 1,])))
        C0 <- do.call(survML, c(nuisance$cens, list(time = time[k != i & treatment == 0], status = 1*(status[k != i & treatment == 0] == 0), X = confounders[k != i & treatment == 0,])))
        phi1 <- phiSurv(surv = S1,
                       cens = C1,
                       treat = pi,
                       newtimes = time,
                       newstatus = status,
                       newX = confounders,
                       newA = treatment,
                       predtimes = evaltimes)$phi1
        phi0 <- phiSurv(surv = S0,
                       cens = C0,
                       treat = pi,
                       newtimes = time,
                       newstatus = status,
                       newX = confounders,
                       newA = treatment,
                       predtimes = evaltimes)$phi0
        phi <- phi1 - phi0
        tau_tmp <- predict(S1, newX = confounders, newtimes = evaltimes)$surv -
          predict(S0, newX = confounders, newtimes = evaltimes)$surv
        #taup <- mean(tau_tmp[k!=i])
        taup <- mean(phi[k!=i])
        tau <- tau_tmp[k == i]
        taus_tmp <- do.call(pseudoReg, c(nuisance$pseudoReg, list(pseudo = tau_tmp[k != i,], X = confounders[k != i, sminus])))
        taus <- predict(taus_tmp, newX = confounders[k == i, sminus])
      } else if(tauLearner == "DR"){
        S <- do.call(survML, c(nuisance$surv, list(time = time[k != i], status = status[k != i], X = confounders[k != i,], A = treatment[k != i])))
        C <- do.call(survML, c(nuisance$cens, list(time = time[k != i], status = 1*(status[k != i] == 0), X = confounders[k != i,], A = treatment[k != i])))
        phi <- phiSurv(surv = S,
                       cens = C,
                       treat = pi,
                       newtimes = time,
                       newstatus = status,
                       newX = confounders,
                       newA = treatment,
                       predtimes = evaltimes)$phi$phi
        tau_tmps <- do.call(drSurv, list(time = time[k != i],
                                         status = status[k != i],
                                         treatment = treatment[k != i],
                                         confounders = confounders[k != i,],
                                         evaltimes = evaltimes,
                                         nuisance = nuisance,
                                         s = 0))
        tau_tmp <- predict(tau_tmps, newX = confounders)
        tau <- tau_tmp$tau[k == i]
        taup <- tau_tmps$taup
        taus_tmp <- do.call(pseudoReg, c(nuisance$pseudoReg, list(pseudo = tau_tmp$tau[k != i], X = confounders[k != i, sminus])))
        taus <- predict(taus_tmp, newX = confounders[k == i, sminus])
      }
      
      varTauPlug <- c(var(tau))
      varTausPlug <- c(var(taus))
      thetasPlug <- varTauPlug - varTausPlug
      psisPlug <- thetasPlug/varTauPlug
      
      #taup <- mean(phi[k != i, ])
      
      eifVarTaus_uncen <- (taup - phi[k == i, ])^2 - (taus - phi[k == i,])^2
      eifVarTau_uncen <- (taup - phi[k == i, ])^2 - (tau - phi[k == i, ])^2
      eifThetas_uncen <- (taus - phi[k == i, ])^2 - (tau - phi[k == i, ])^2
      
      varTau <- mean(eifVarTau_uncen)
      varTaus <- mean(eifVarTaus_uncen)
      thetas <- mean(eifThetas_uncen)
      psis <- thetas/varTau
      taup <- mean(phi[k == i])
      
      #vals <- list(eifVarTau_uncen = eifVarTau_uncen, eifVarTaus_uncen = eifVarTaus_uncen, eifThetas_uncen = eifThetas_uncen, phi = phi[k == i, ], tau = tau, taus = taus)
      vals <- list(varTau = varTau, varTaus = varTaus, thetas = thetas, psis = psis, taup = taup, phi = phi, thetasPlug = thetasPlug, psisPlug = psisPlug, varTauPlug = varTauPlug,
                   eifVarTau_uncen = eifVarTau_uncen, eifVarTaus_uncen = eifVarTaus_uncen, eifThetas_uncen = eifThetas_uncen,
                   nk = sum(k == i))
      vals
    })
    # eifVarTau_uncen <- unlist(lapply(1:folds, function(i){
    #   vals[[i]]$eifVarTau_uncen
    # }))
    # eifVarTaus_uncen <- unlist(lapply(1:folds, function(i){
    #   vals[[i]]$eifVarTaus_uncen
    # }))
    # eifThetas_uncen <- unlist(lapply(1:folds, function(i){
    #   vals[[i]]$eifThetas_uncen
    # }))
    # eifTaup_uncen <- unlist(lapply(1:folds, function(i){
    #   vals[[i]]$phi
    # }))
    # tau <- unlist(lapply(1:folds, function(i){
    #   vals[[i]]$tau
    # }))
    # taus <- unlist(lapply(1:folds, function(i){
    #   vals[[i]]$taus
    # }))
    varTau <- sum(sapply(vals, function(x) x$nk/N * x$varTau))
    varTaus <- sum(sapply(vals, function(x) x$nk/N * x$varTaus))
    thetas <- sum(sapply(vals, function(x) x$nk/N * x$thetas))
    taup <- sum(sapply(vals, function(x) x$nk/N * x$taup))
    psis <- thetas/varTau
    
    var_varTau <- sum(sapply(vals, function(x) x$nk/N * mean((x$eifVarTau_uncen - varTau)^2))) / N
    var_varTaus <- sum(sapply(vals, function(x) x$nk/N * mean((x$eifVarTaus_uncen - varTaus)^2))) / N
    var_thetas <- sum(sapply(vals, function(x) x$nk/N * mean((x$eifThetas_uncen - thetas)^2))) / N
    var_taup <- sum(sapply(vals, function(x) x$nk/N * mean((x$phi - taup)^2))) / N
    var_psis <- sum(sapply(vals, function(x) x$nk/N * mean(((x$eifThetas_uncen - thetas - psis*(x$eifVarTau_uncen - varTau))/varTau)^2))) / N
    
    thetasLog <- sum(sapply(vals, function(x) x$nk/N * (log(x$thetasPlug) + (x$thetas - x$thetasPlug))))
    varTauLog <- sum(sapply(vals, function(x) x$nk/N * (log(x$varTauPlug) + (x$varTau - x$varTauPlug))))
    psisLogit <- sum(sapply(vals, function(x) x$nk/N * (logit(x$psisPlug) + (x$varTau/x$varTauPlug) * ((x$psis - x$psisPlug)/(x$psisPlug - x$psisPlug^2)))))
    var_thetasLog <- var_thetas/exp(thetasLog)^2
    var_varTauLog <- var_varTau/exp(varTauLog)^2
    var_psisLogit <- var_psis/(expit(psisLogit) - expit(psisLogit)^2)^2
  }
  
  # varTau <- mean(eifVarTau_uncen)
  # varTaus <- mean(eifVarTaus_uncen)
  # thetas <- mean(eifThetas_uncen)
  # psis <- thetas/varTau
  # taup <- mean(eifTaup_uncen)
  # 
  # eifVarTau <- eifVarTau_uncen - varTau
  # eifVarTaus <- eifVarTaus_uncen - varTaus
  # eifThetas <- eifThetas_uncen - thetas
  # eifPsis <- (eifThetas - psis*eifVarTau)/varTau
  # eifTaup <- eifTaup_uncen - taup
  # 
  # var_varTau <- (1/N)*mean(eifVarTau^2)
  # var_varTaus <- (1/N)*mean(eifVarTaus^2)
  # var_thetas <- (1/N)*mean(eifThetas^2)
  # var_psis <- (1/N)*mean(eifPsis^2)
  # var_taup <- (1/N)*mean(eifTaup^2)
  
  #eifTaup <- phi - taup
  #seTaup <- sqrt((1/N)*mean(eifTaup^2))
  
  list(varTau = varTau,
       seVarTau = sqrt(var_varTau),
       varTaus = varTaus,
       seVarTaus = sqrt(var_varTaus),
       thetas = thetas,
       seThetas = sqrt(var_thetas),
       psis = psis,
       sePsis = sqrt(var_psis),
       ate = taup,
       seAte = sqrt(var_taup),
       thetasLog = thetasLog,
       seThetasLog = sqrt(var_thetasLog),
       varTauLog = varTauLog,
       seVarTauLog = sqrt(var_varTauLog),
       psisLogit = psisLogit,
       sePsisLogit = sqrt(var_psisLogit),
       treatment_diff = paste0(rev(treatLabs), collapse = " - "))
  
}




























teVimSurv_tmp2 <- function(time,
                          status,
                          treatment,
                          confounders,
                          s,
                          evaltimes,
                          nuisance = list(),
                          tauLearner = "S",
                          dr = FALSE,
                          taupFolds = 1,
                          folds = 1){
  #browser()
  confounders <- as.data.frame(confounders)
  N <- length(time)
  sminus <- setdiff(1:dim(confounders)[2], s)
  k <- sample(1:folds, length(status), replace = TRUE)
  if(folds == 1){
    pi <- do.call(propML, c(nuisance$prop, list(X = confounders, A = treatment)))
    if(tauLearner == "S"){
      S <- do.call(survML, c(nuisance$surv, list(time = time, status = status, X = confounders, A = treatment)))
      C <- do.call(survML, c(nuisance$cens, list(time = time, status = 1*(status==0), X = confounders, A = treatment)))
      phi <- phiSurv(surv = S,
                     cens = C,
                     treat = pi,
                     newtimes = time,
                     newstatus = status,
                     newX = confounders,
                     newA = treatment,
                     predtimes = evaltimes)$phi$phi
      tau <- predict(S, newX = confounders, newtimes = evaltimes, newA = 1)$surv -
        predict(S, newX = confounders, newtimes = evaltimes, newA = 0)$surv
      if(taupFolds == 1){
        taup <- mean(phi)
      } else{
        k_taup <- sample(1:taupFolds, length(status), replace = TRUE)
        valsTaup <- lapply(1:taupFolds, function(j){
          pi <- do.call(propML, c(nuisance$prop, list(X = confounders[k_taup != j, ], A = treatment[k_taup != j])))
          S <- do.call(survML, c(nuisance$surv, list(time = time[k_taup != j], status = status[k_taup != j], X = confounders[k_taup != j,], A = treatment[k_taup != j])))
          C <- do.call(survML, c(nuisance$cens, list(time = time[k_taup != j], status = 1*(status[k_taup != j] == 0), X = confounders[k_taup != j,], A = treatment[k_taup != j])))
          phi <- phiSurv(surv = S,
                         cens = C,
                         treat = pi,
                         newtimes = time[k_taup == j],
                         newstatus = status[k_taup == j],
                         newX = confounders[k_taup == j, ],
                         newA = treatment[k_taup == j],
                         predtimes = evaltimes)$phi$phi
          list(phi=phi)
        })
        taup <- mean(unlist(phi))
      }
      
      taus_tmp <- do.call(pseudoReg, c(nuisance$pseudoReg, list(pseudo = tau, X = confounders[, sminus, drop = FALSE])))
      taus <- predict(taus_tmp, newX = confounders[, sminus, drop = FALSE])
      varTauPlug <- var(tau)
      varTausPlug <- var(taus)
      thetasPlug <- varTauPlug - varTausPlug
      psisPlug <- thetasPlug/varTauPlug
    } else if(tauLearner == "T"){
      S1 <- do.call(survML, c(nuisance$surv, list(time = time[treatment == 1], status = status[treatment == 1], X = confounders[treatment == 1, ])))
      S0 <- do.call(survML, c(nuisance$surv, list(time = time[treatment == 0], status = status[treatment == 0], X = confounders[treatment == 0, ])))
      C1 <- do.call(survML, c(nuisance$cens, list(time = time[treatment == 1], status = 1*(status[treatment == 1]==0), X = confounders[treatment == 1, ])))
      C0 <- do.call(survML, c(nuisance$cens, list(time = time[treatment == 0], status = 1*(status[treatment == 0]==0), X = confounders[treatment == 0, ])))
      phi1 <- phiSurv(surv = S1,
                      cens = C1,
                      treat = pi,
                      newtimes = time,
                      newstatus = status,
                      newX = confounders,
                      newA = treatment,
                      predtimes = evaltimes)$phi$phi1
      phi0 <- phiSurv(surv = S0,
                      cens = C0,
                      treat = pi,
                      newtimes = time,
                      newstatus = status,
                      newX = confounders,
                      newA = treatment,
                      predtimes = evaltimes)$phi$phi0
      phi <- phi1 - phi0
      tau <- predict(S1, newX = confounders, newtimes = evaltimes)$surv -
        predict(S0, newX = confounders, newtimes = evaltimes)$surv
      #taup <- mean(tau)
      taup <- mean(phi)
      taus_tmp <- do.call(pseudoReg, c(nuisance$pseudoReg, list(pseudo = tau, X = confounders[, sminus, drop = FALSE])))
      taus <- predict(taus_tmp, newX = confounders[, sminus, drop = FALSE])
      varTauPlug <- c(var(tau))
      varTausPlug <- c(var(taus))
      thetasPlug <- varTauPlug - varTausPlug
      psisPlug <- thetasPlug/varTauPlug
    } else if(tauLearner == "DR"){
      S <- do.call(survML, c(nuisance$surv, list(time = time, status = status, X = confounders, A = treatment)))
      C <- do.call(survML, c(nuisance$cens, list(time = time, status = 1*(status==0), X = confounders, A = treatment)))
      phi <- phiSurv(surv = S,
                     cens = C,
                     treat = pi,
                     newtimes = time,
                     newstatus = status,
                     newX = confounders,
                     newA = treatment,
                     predtimes = evaltimes)$phi$phi
      tau_tmps <- do.call(drSurv, list(time = time,
                                       status = status,
                                       treatment = treatment,
                                       confounders = confounders,
                                       evaltimes = evaltimes,
                                       nuisance = nuisance,
                                       s = 0))
      tau_tmp <- predict(tau_tmps, newX = confounders)
      tau <- tau_tmp$tau
      taup <- tau_tmps$taup
      taus_tmp <- do.call(pseudoReg, c(nuisance$pseudoReg, list(pseudo = tau, X = confounders[, sminus, drop = FALSE])))
      taus <- predict(taus_tmp, newX = confounders[, sminus, drop = FALSE])
      varTauPlug <- var(tau)
      varTausPlug <- var(taus)
      thetasPlug <- varTauPlug - varTausPlug
      psisPlug <- thetasPlug/varTauPlug
    }
    
    #taup <- mean(phi)
    
    eifVarTaus_uncen <- (taup - phi)^2 - (taus - phi)^2
    eifVarTau_uncen <- (taup - phi)^2 - (tau - phi)^2
    eifThetas_uncen <- (taus - phi)^2 - (tau - phi)^2
    eifTaup_uncen <- phi
    
    varTau <- mean(eifVarTau_uncen)
    varTaus <- mean(eifVarTaus_uncen)
    thetas <- mean(eifThetas_uncen)
    psis <- thetas/varTau
    taup <- mean(eifTaup_uncen)
    
    eifVarTau <- eifVarTau_uncen - varTau
    eifVarTaus <- eifVarTaus_uncen - varTaus
    eifThetas <- eifThetas_uncen - thetas
    eifPsis <- (eifThetas - psis*eifVarTau)/varTau
    eifTaup <- eifTaup_uncen - taup
    
    var_varTau <- (1/N)*mean(eifVarTau^2)
    var_varTaus <- (1/N)*mean(eifVarTaus^2)
    var_thetas <- (1/N)*mean(eifThetas^2)
    var_psis <- (1/N)*mean(eifPsis^2)
    var_taup <- (1/N)*mean(eifTaup^2)
    
    thetasLog <- log(thetasPlug) + (thetas - thetasPlug)/thetasPlug
    varTauLog <- log(varTauPlug) + (varTau - varTauPlug)/varTauPlug
    psisLogit <- logit(psisPlug) + (varTau/varTauPlug) * ((psis - psisPlug)/(psisPlug - psisPlug^2))
    var_thetasLog <- var_thetas/exp(thetasLog)^2
    var_varTauLog <- var_varTau/exp(varTauLog)^2
    var_psisLogit <- var_psis/(expit(psisLogit) - expit(psisLogit)^2)^2
    
  } else{
    vals <- lapply(1:folds, function(i){
      pi <- do.call(propML, c(nuisance$prop, list(X = confounders[k != i,], A = treatment[k != i])))
      if(tauLearner == "S"){
        S <- do.call(survML, c(nuisance$surv, list(time = time[k != i], status = status[k != i], X = confounders[k != i,], A = treatment[k != i])))
        C <- do.call(survML, c(nuisance$cens, list(time = time[k != i], status = 1*(status[k != i] == 0), X = confounders[k != i,], A = treatment[k != i])))
        phi <- phiSurv(surv = S,
                       cens = C,
                       treat = pi,
                       newtimes = time,
                       newstatus = status,
                       newX = confounders,
                       newA = treatment,
                       predtimes = evaltimes)$phi$phi
        tau_tmp <- predict(S, newX = confounders, newtimes = evaltimes, newA = 1)$surv -
          predict(S, newX = confounders, newtimes = evaltimes, newA = 0)$surv
        #taup <- mean(tau_tmp[k!=i])
        taup <- mean(phi[k!=i])
        tau <- tau_tmp[k == i]
        taus_tmp <- do.call(pseudoReg, c(nuisance$pseudoReg, list(pseudo = tau_tmp[k != i,], X = confounders[k != i, sminus])))
        taus <- predict(taus_tmp, newX = confounders[k == i, sminus])
      } else if(tauLearner == "T"){
        S1 <- do.call(survML, c(nuisance$surv, list(time = time[k != i & treatment == 1], status = status[k != i & treatment == 1], X = confounders[k != i & treatment == 1,])))
        S0 <- do.call(survML, c(nuisance$surv, list(time = time[k != i & treatment == 0], status = status[k != i & treatment == 0], X = confounders[k != i & treatment == 0,])))
        C1 <- do.call(survML, c(nuisance$cens, list(time = time[k != i & treatment == 1], status = 1*(status[k != i & treatment == 1] == 0), X = confounders[k != i & treatment == 1,])))
        C0 <- do.call(survML, c(nuisance$cens, list(time = time[k != i & treatment == 0], status = 1*(status[k != i & treatment == 0] == 0), X = confounders[k != i & treatment == 0,])))
        phi1 <- phiSurv(surv = S1,
                        cens = C1,
                        treat = pi,
                        newtimes = time,
                        newstatus = status,
                        newX = confounders,
                        newA = treatment,
                        predtimes = evaltimes)$phi$phi1
        phi0 <- phiSurv(surv = S0,
                        cens = C0,
                        treat = pi,
                        newtimes = time,
                        newstatus = status,
                        newX = confounders,
                        newA = treatment,
                        predtimes = evaltimes)$phi$phi0
        phi <- phi1 - phi0
        tau_tmp <- predict(S1, newX = confounders, newtimes = evaltimes)$surv -
          predict(S0, newX = confounders, newtimes = evaltimes)$surv
        #taup <- mean(tau_tmp[k!=i])
        taup <- mean(phi[k!=i])
        tau <- tau_tmp[k == i]
        taus_tmp <- do.call(pseudoReg, c(nuisance$pseudoReg, list(pseudo = tau_tmp[k != i,], X = confounders[k != i, sminus])))
        taus <- predict(taus_tmp, newX = confounders[k == i, sminus])
      } else if(tauLearner == "DR"){
        S <- do.call(survML, c(nuisance$surv, list(time = time[k != i], status = status[k != i], X = confounders[k != i,], A = treatment[k != i])))
        C <- do.call(survML, c(nuisance$cens, list(time = time[k != i], status = 1*(status[k != i] == 0), X = confounders[k != i,], A = treatment[k != i])))
        phi <- phiSurv(surv = S,
                       cens = C,
                       treat = pi,
                       newtimes = time,
                       newstatus = status,
                       newX = confounders,
                       newA = treatment,
                       predtimes = evaltimes)$phi$phi
        tau_tmps <- do.call(drSurv, list(time = time[k != i],
                                         status = status[k != i],
                                         treatment = treatment[k != i],
                                         confounders = confounders[k != i,],
                                         evaltimes = evaltimes,
                                         nuisance = nuisance,
                                         s = 0))
        tau_tmp <- predict(tau_tmps, newX = confounders)
        tau <- tau_tmp$tau[k == i]
        taup <- tau_tmps$taup
        taus_tmp <- do.call(pseudoReg, c(nuisance$pseudoReg, list(pseudo = tau_tmp$tau[k != i], X = confounders[k != i, sminus])))
        taus <- predict(taus_tmp, newX = confounders[k == i, sminus])
      }
      
      varTauPlug <- c(var(tau))
      varTausPlug <- c(var(taus))
      thetasPlug <- varTauPlug - varTausPlug
      psisPlug <- thetasPlug/varTauPlug
      
      #taup <- mean(phi[k != i, ])
      
      eifVarTaus_uncen <- (taup - phi[k == i, ])^2 - (taus - phi[k == i,])^2
      eifVarTau_uncen <- (taup - phi[k == i, ])^2 - (tau - phi[k == i, ])^2
      eifThetas_uncen <- (taus - phi[k == i, ])^2 - (tau - phi[k == i, ])^2
      
      varTau <- mean(eifVarTau_uncen)
      varTaus <- mean(eifVarTaus_uncen)
      thetas <- mean(eifThetas_uncen)
      psis <- thetas/varTau
      taup <- mean(phi[k == i])
      
      #vals <- list(eifVarTau_uncen = eifVarTau_uncen, eifVarTaus_uncen = eifVarTaus_uncen, eifThetas_uncen = eifThetas_uncen, phi = phi[k == i, ], tau = tau, taus = taus)
      vals <- list(varTau = varTau, varTaus = varTaus, thetas = thetas, psis = psis, taup = taup, phi = phi, thetasPlug = thetasPlug, psisPlug = psisPlug, varTauPlug = varTauPlug,
                   eifVarTau_uncen = eifVarTau_uncen, eifVarTaus_uncen = eifVarTaus_uncen, eifThetas_uncen = eifThetas_uncen,
                   nk = sum(k == i))
      vals
    })
    # eifVarTau_uncen <- unlist(lapply(1:folds, function(i){
    #   vals[[i]]$eifVarTau_uncen
    # }))
    # eifVarTaus_uncen <- unlist(lapply(1:folds, function(i){
    #   vals[[i]]$eifVarTaus_uncen
    # }))
    # eifThetas_uncen <- unlist(lapply(1:folds, function(i){
    #   vals[[i]]$eifThetas_uncen
    # }))
    # eifTaup_uncen <- unlist(lapply(1:folds, function(i){
    #   vals[[i]]$phi
    # }))
    # tau <- unlist(lapply(1:folds, function(i){
    #   vals[[i]]$tau
    # }))
    # taus <- unlist(lapply(1:folds, function(i){
    #   vals[[i]]$taus
    # }))
    varTau <- sum(sapply(vals, function(x) x$nk/N * x$varTau))
    varTaus <- sum(sapply(vals, function(x) x$nk/N * x$varTaus))
    thetas <- sum(sapply(vals, function(x) x$nk/N * x$thetas))
    taup <- sum(sapply(vals, function(x) x$nk/N * x$taup))
    psis <- thetas/varTau
    
    var_varTau <- sum(sapply(vals, function(x) x$nk/N * mean((x$eifVarTau_uncen - varTau)^2))) / N
    var_varTaus <- sum(sapply(vals, function(x) x$nk/N * mean((x$eifVarTaus_uncen - varTaus)^2))) / N
    var_thetas <- sum(sapply(vals, function(x) x$nk/N * mean((x$eifThetas_uncen - thetas)^2))) / N
    var_taup <- sum(sapply(vals, function(x) x$nk/N * mean((x$phi - taup)^2))) / N
    var_psis <- sum(sapply(vals, function(x) x$nk/N * mean(((x$eifThetas_uncen - thetas - psis*(x$eifVarTau_uncen - varTau))/varTau)^2))) / N
    
    thetasLog <- sum(sapply(vals, function(x) x$nk/N * (log(x$thetasPlug) + (x$thetas - x$thetasPlug))))
    varTauLog <- sum(sapply(vals, function(x) x$nk/N * (log(x$varTauPlug) + (x$varTau - x$varTauPlug))))
    psisLogit <- sum(sapply(vals, function(x) x$nk/N * (logit(x$psisPlug) + (x$varTau/x$varTauPlug) * ((x$psis - x$psisPlug)/(x$psisPlug - x$psisPlug^2)))))
    var_thetasLog <- var_thetas/exp(thetasLog)^2
    var_varTauLog <- var_varTau/exp(varTauLog)^2
    var_psisLogit <- var_psis/(expit(psisLogit) - expit(psisLogit)^2)^2
  }
  
  # varTau <- mean(eifVarTau_uncen)
  # varTaus <- mean(eifVarTaus_uncen)
  # thetas <- mean(eifThetas_uncen)
  # psis <- thetas/varTau
  # taup <- mean(eifTaup_uncen)
  # 
  # eifVarTau <- eifVarTau_uncen - varTau
  # eifVarTaus <- eifVarTaus_uncen - varTaus
  # eifThetas <- eifThetas_uncen - thetas
  # eifPsis <- (eifThetas - psis*eifVarTau)/varTau
  # eifTaup <- eifTaup_uncen - taup
  # 
  # var_varTau <- (1/N)*mean(eifVarTau^2)
  # var_varTaus <- (1/N)*mean(eifVarTaus^2)
  # var_thetas <- (1/N)*mean(eifThetas^2)
  # var_psis <- (1/N)*mean(eifPsis^2)
  # var_taup <- (1/N)*mean(eifTaup^2)
  
  #eifTaup <- phi - taup
  #seTaup <- sqrt((1/N)*mean(eifTaup^2))
  
  list(varTau = varTau,
       seVarTau = sqrt(var_varTau),
       varTaus = varTaus,
       seVarTaus = sqrt(var_varTaus),
       thetas = thetas,
       seThetas = sqrt(var_thetas),
       psis = psis,
       sePsis = sqrt(var_psis),
       ate = taup,
       seAte = sqrt(var_taup),
       thetasLog = thetasLog,
       seThetasLog = sqrt(var_thetasLog),
       varTauLog = varTauLog,
       seVarTauLog = sqrt(var_varTauLog),
       psisLogit = psisLogit,
       sePsisLogit = sqrt(var_psisLogit))
  
}

logit <- function(x){
  log(x/(1 - x))
}

expit <- function(x){
  exp(x)/(1 + exp(x))
}


