covVimSurv <- function(time,
                          status,
                          treatment,
                          confounders,
                          l,
                          evaltimes,
                          nuisance = list(),
                          tauLearner = "S",
                          folds = 1){
  #browser()
  require(riskRegression)
  A <- treatment
  treatLabs <- levels(treatment)
  treatment <- 1*(treatment == treatLabs[2])
  confounders <- as.data.frame(confounders)
  xl_numeric <- list()
  for(i in seq_along(l)){
    if(class(confounders[, i]) %in% c("numeric", "integer")){
      xl_numeric[[i]] <- confounders[, i]
    } else if(class(confounders[, i]) == "factor"){
      xl_numeric_labs <- levels(confounders[, i])
      xl_numeric[[i]] <- 1*(confounders[, i] == xl_numeric_labs[2])
    }
  }
  
  N <- length(time)
  lminus <- lapply(l, function(x){setdiff(1:dim(confounders)[2], x)})
  k <- sample(1:folds, length(status), replace = TRUE)
  if(folds == 1){
    pi <- do.call(propML, c(nuisance$prop, list(X = confounders, A = A)))
    if(tauLearner == "S"){
      S <- do.call(survML, c(nuisance$surv, list(time = time, status = 1*(status == 1), X = confounders, A = treatment)))
      C <- do.call(survML, c(nuisance$cens, list(time = time, status = 1*(status == 0), X = confounders, A = treatment)))
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
      taup <- mean(phi)
      
      taul_tmp <- lapply(seq_along(l), function(x){
        do.call(pseudoReg, c(nuisance$pseudoRegTau, list(pseudo = tau[,1], X = confounders[, lminus[[x]], drop = FALSE])))
      })
      taul <- lapply(seq_along(l), function(x){
        predict(taul_tmp[[x]], newX = confounders[, lminus[[x]], drop = FALSE])
      })
      
      EXl_Xlminus_tmp <- lapply(seq_along(l), function(x){
        do.call(pseudoReg, c(nuisance$pseudoRegX, list(pseudo = confounders[, l[x]], X = confounders[, lminus[[x]], drop = FALSE])))
      })
      EXl_Xlminus <- lapply(seq_along(l), function(x){
        predict(EXl_Xlminus_tmp[[x]], newX = confounders[, lminus[[x]], drop = FALSE])
      })
      
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
      taup <- mean(phi)
      taul_tmp <- lapply(seq_along(l), function(x){
        do.call(pseudoReg, c(nuisance$pseudoRegTau, list(pseudo = tau[,1], X = confounders[, lminus[[x]], drop = FALSE])))
      })
      taul <- lapply(seq_along(l), function(x){
        predict(taul_tmp[[x]], newX = confounders[, lminus[[x]], drop = FALSE])
      })
      
      EXl_Xlminus_tmp <- lapply(seq_along(l), function(x){
        do.call(pseudoReg, c(nuisance$pseudoRegX, list(pseudo = confounders[, l[x]], X = confounders[, lminus[[x]], drop = FALSE])))
      })
      EXl_Xlminus <- lapply(seq_along(l), function(x){
        predict(EXl_Xlminus_tmp[[x]], newX = confounders[, lminus[[x]], drop = FALSE])
      })
    }
    
    res_tmp <- lapply(seq_along(l), function(x){
      eif_Gamma_j_uncen <- (phi - taul[[x]])*(xl_numeric[[x]] - EXl_Xlminus[[x]])
      Gamma_j <- mean(eif_Gamma_j_uncen)
      eif_Gamma_j <- eif_Gamma_j_uncen - Gamma_j
      var_Gamma_j <- (1/N)*mean(eif_Gamma_j^2)
      se_Gamma_j <- sqrt(var_Gamma_j)
      pval_Gamma_j <- 2*(1 - pnorm(abs(Gamma_j)/sqrt(var_Gamma_j)))
      
      eif_chi_j_uncen <- (xl_numeric[[x]] - EXl_Xlminus[[x]])^2
      chi_j <- mean(eif_chi_j_uncen)
      eif_chi_j <- eif_chi_j_uncen - chi_j
      var_chi_j <- (1/N)*mean(eif_chi_j^2)
      se_chi_j <- sqrt(var_chi_j)
      
      omega_t <- Gamma_j/chi_j
      eif_omega_t <- (eif_Gamma_j - omega_t*eif_chi_j)/chi_j
      var_omega_t <- (1/N)*mean(eif_omega_t^2)
      se_omega_t <- sqrt(var_omega_t)
      pval_omega_t <- 2*(1 - pnorm(abs(omega_t)/sqrt(var_omega_t)))
      
      data.frame(l = names(confounders)[l[x]],
                 Gamma_j = Gamma_j, se_Gamma_j = se_Gamma_j, pval_Gamma_j = pval_Gamma_j,
                 chi_j = chi_j, se_chi_j = se_chi_j,
                 omega_t = omega_t, se_omega_t = se_omega_t, pval_omega_t = pval_omega_t)
    })
    res <- do.call(rbind, res_tmp)
    
    eif_taup_uncen <- phi
    taup <- mean(phi)
    eif_taup <- eif_taup_uncen - taup
    var_taup <- (1/N)*mean(eif_taup^2)
    
    eif_est.1_uncen <- phi1
    est.1 <- mean(eif_est.1_uncen)
    eif_est.1 <- eif_est.1_uncen - est.1
    var_est.1 <- (1/N)*mean(eif_est.1^2)
    
    eif_est.0_uncen <- phi0
    est.0 <- mean(eif_est.0_uncen)
    eif_est.0 <- eif_est.0_uncen - est.0
    var_est.0 <- (1/N)*mean(eif_est.0^2)
  } else{
    vals <- lapply(1:folds, function(i){
      pi <- do.call(propML, c(nuisance$prop, list(X = confounders[k != i,], A = A[k != i])))
      if(tauLearner == "S"){
        S <- do.call(survML, c(nuisance$surv, list(time = time[k != i], status = status[k != i], X = confounders[k != i,], A = treatment[k != i])))
        C <- do.call(survML, c(nuisance$cens, list(time = time[k != i], status = 1*(status[k != i] == 0), X = confounders[k != i,], A = treatment[k != i])))
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
        tau_tmp <- predict(S, newX = confounders, newtimes = evaltimes, newA = 1)$surv -
          predict(S, newX = confounders, newtimes = evaltimes, newA = 0)$surv
        tau <- tau_tmp[k == i]
        #taup <- mean(tau_tmp[k != i])
        taup <- mean(phi[k != i])
        
        taul_tmp <- lapply(seq_along(l), function(x){
          do.call(pseudoReg, c(nuisance$pseudoRegTau, list(pseudo = tau_tmp[k != i, 1], X = confounders[k != i, lminus[[x]], drop = FALSE])))
        })
        taul <- lapply(seq_along(l), function(x){
          predict(taul_tmp[[x]], newX = confounders[k == i, lminus[[x]], drop = FALSE])
        })
        
        EXl_Xlminus_tmp <- lapply(seq_along(l), function(x){
          do.call(pseudoReg, c(nuisance$pseudoRegX, list(pseudo = confounders[k != i, l[x]], X = confounders[k != i, lminus[[x]], drop = FALSE])))
        })
        EXl_Xlminus <- lapply(seq_along(l), function(x){
          predict(EXl_Xlminus_tmp[[x]], newX = confounders[k == i, lminus[[x]], drop = FALSE])
        })
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
        tau <- tau_tmp[k == i]
        
        taul_tmp <- lapply(seq_along(l), function(x){
          do.call(pseudoReg, c(nuisance$pseudoRegTau, list(pseudo = tau_tmp[k != i, 1], X = confounders[k != i, lminus[[x]], drop = FALSE])))
        })
        taul <- lapply(seq_along(l), function(x){
          predict(taul_tmp[[x]], newX = confounders[k == i, lminus[[x]], drop = FALSE])
        })
        
        EXl_Xlminus_tmp <- lapply(seq_along(l), function(x){
          do.call(pseudoReg, c(nuisance$pseudoRegX, list(pseudo = confounders[k != i, l[x]], X = confounders[k != i, lminus[[x]], drop = FALSE])))
        })
        EXl_Xlminus <- lapply(seq_along(l), function(x){
          predict(EXl_Xlminus_tmp[[x]], newX = confounders[k == i, lminus[[x]], drop = FALSE])
        })
      }
      
      eif_Gamma_j_uncen <- lapply(seq_along(l), function(x){
        (phi[k == i, ] - taul[[x]])*(xl_numeric[[x]][k == i] - EXl_Xlminus[[x]])
      })
      Gamma_j <- lapply(eif_Gamma_j_uncen, mean)
      
      eif_chi_j_uncen <- lapply(seq_along(l), function(x){
        (xl_numeric[[x]][k == i] - EXl_Xlminus[[x]])^2
      })  
      chi_j <- lapply(eif_chi_j_uncen, mean)
      
      eif_taup_uncen <- phi[k == i, ]
      taup <- mean(eif_taup_uncen)
      
      eif_est.1_uncen <- phi1
      eif_est.0_uncen <- phi0
      
      est.1 <- mean(eif_est.1_uncen[k == i, ])
      est.0 <- mean(eif_est.0_uncen[k == i, ])
      
      vals <- list(nk = sum(k == i), eif_Gamma_j_uncen = eif_Gamma_j_uncen,
                   eif_taup_uncen = eif_taup_uncen, eif_chi_j_uncen = eif_chi_j_uncen,
                   Gamma_j = Gamma_j, chi_j = chi_j,
                   taup = taup, est.1 = est.1, est.0 = est.0,
                   eif_est.1_uncen  = eif_est.1_uncen[k == i, ], eif_est.0_uncen = eif_est.0_uncen[k == i, ])
      vals
    })
    
    res_tmp <- lapply(seq_along(l), function(j){
      Gamma_j <- sum(sapply(vals, function(x) x$nk/N * x$Gamma_j[[j]]))
      chi_j <- sum(sapply(vals, function(x) x$nk/N * x$chi_j[[j]]))
      omega_t <- Gamma_j/chi_j
      
      var_Gamma_j <- sum(sapply(vals, function(x) x$nk/N * mean((x$eif_Gamma_j_uncen[[j]] - Gamma_j)^2)))/N
      var_chi_j <- sum(sapply(vals, function(x) x$nk/N * mean((x$eif_chi_j_uncen[[j]] - chi_j)^2)))/N
      var_omega_t <- sum(sapply(vals, function(x) x$nk/N * mean(((x$eif_Gamma_j_uncen[[j]] - Gamma_j - omega_t*(x$eif_chi_j_uncen[[j]] - chi_j))/chi_j)^2)))/N
      
      pval_Gamma_j <- 2*(1 - pnorm(abs(Gamma_j)/sqrt(var_Gamma_j)))
      pval_omega_t <- 2*(1 - pnorm(abs(omega_t)/sqrt(var_omega_t)))
      data.frame(l = names(confounders)[l[j]],
                 Gamma_j = Gamma_j, se_Gamma_j = sqrt(var_Gamma_j), pval_Gamma_j = pval_Gamma_j,
                 chi_j = chi_j, se_chi_j = sqrt(var_chi_j),
                 omega_t, se_omega_t = sqrt(var_omega_t), pval_omega_t = pval_omega_t)
    })
    res <- do.call(rbind, res_tmp)
    
    taup <- sum(sapply(vals, function(x) x$nk/N * x$taup))
    var_taup <- sum(sapply(vals, function(x) x$nk/N * mean((x$eif_taup_uncen - taup)^2)))/N
    
    est.1 <- sum(sapply(vals, function(x) x$nk/N * x$est.1))
    var_est.1 <- sum(sapply(vals, function(x) x$nk/N * mean((x$eif_est.1_uncen - est.1)^2)))/N
    
    est.0 <- sum(sapply(vals, function(x) x$nk/N * x$est.0))
    var_est.0 <- sum(sapply(vals, function(x) x$nk/N * mean((x$eif_est.0_uncen - est.0)^2)))/N
    
    
  }
  
  list(modifiers = res,
       ate = data.frame(taup = taup, seAte = sqrt(var_taup),
                        est.1 = est.1, se_est.1 = sqrt(var_est.1),
                        est.0 = est.0, se_est.0 = sqrt(var_est.0)),
       treatment_diff = paste0(rev(treatLabs), collapse = " - "))
  
}





















covVimSurv_old <- function(time,
                       status,
                       treatment,
                       confounders,
                       s,
                       evaltimes,
                       nuisance = list(),
                       tauLearner = "S",
                       folds = 1){
  #browser()
  A <- treatment
  treatLabs <- levels(treatment)
  treatment <- 1*(treatment == treatLabs[2])
  confounders <- as.data.frame(confounders)
  if(class(confounders[, s]) == "numeric"){
    xs_numeric <- confounders[, s]
  } else if(class(confounders[, s]) == "factor"){
    xs_numeric_labs <- levels(confounders[, s])
    xs_numeric <- 1*(confounders[, s] == xs_numeric_labs[2])
  }
  N <- length(time)
  sminus <- setdiff(1:dim(confounders)[2], s)
  k <- sample(1:folds, length(status), replace = TRUE)
  if(folds == 1){
    pi <- do.call(propML, c(nuisance$prop, list(X = confounders, A = A)))
    if(tauLearner == "S"){
      S <- do.call(survML, c(nuisance$surv, list(time = time, status = status, X = confounders, A = treatment)))
      C <- do.call(survML, c(nuisance$cens, list(time = time, status = 1*(status == 0), X = confounders, A = treatment)))
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
      taup <- mean(phi)
      taus_tmp <- do.call(pseudoReg, c(nuisance$pseudoReg, list(pseudo = tau, X = confounders[, sminus, drop = FALSE])))
      taus <- predict(taus_tmp, newX = confounders[, sminus, drop = FALSE])
      EXs_Xsminus_tmp <- do.call(pseudoReg, c(nuisance$pseudoReg, list(pseudo = confounders[, s], X = confounders[, sminus, drop = FALSE])))
      EXs_Xsminus <- predict(EXs_Xsminus_tmp, newX = confounders[, sminus, drop = FALSE])
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
      taup <- mean(phi)
      taus_tmp <- do.call(pseudoReg, c(nuisance$pseudoReg, list(pseudo = tau, X = confounders[, sminus, drop = FALSE])))
      taus <- predict(taus_tmp, newX = confounders[, sminus, drop = FALSE])
      EXs_Xsminus_tmp <- do.call(pseudoReg, c(nuisance$pseudoReg, list(pseudo = confounders[, s], X = confounders[, sminus, drop = FALSE])))
      EXs_Xsminus <- predict(EXs_Xsminus_tmp, newX = confounders[, sminus, drop = FALSE])
    }
    
    eif_psi_t_uncen <- (phi - taus)*(xs_numeric - EXs_Xsminus)
    psi_t <- mean(eif_psi_t_uncen)
    eif_psi_t <- eif_psi_t_uncen - psi_t
    var_psi_t <- (1/N)*mean(eif_psi_t^2)
    
    eif_taup_uncen <- phi
    taup <- mean(phi)
    eif_taup <- eif_taup_uncen - taup
    var_taup <- (1/N)*mean(eif_taup^2)
    
    eif_phi_t_uncen <- (xs_numeric - EXs_Xsminus)^2
    phi_t <- mean(eif_phi_t_uncen)
    eif_phi_t <- eif_phi_t_uncen - phi_t
    var_phi_t <- (1/N)*mean(eif_phi_t^2)
    
    omega_t <- psi_t/phi_t
    eif_omega_t <- (eif_psi_t - omega_t*eif_phi_t)/phi_t
    var_omega_t <- (1/N)*mean(eif_omega_t^2)
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
                       predtimes = evaltimes)$phi$phi
        tau_tmp <- predict(S, newX = confounders, newtimes = evaltimes, newA = 1)$surv -
          predict(S, newX = confounders, newtimes = evaltimes, newA = 0)$surv
        tau <- tau_tmp[k == i]
        #taup <- mean(tau_tmp[k != i])
        taup <- mean(phi[k != i])
        taus_tmp <- do.call(pseudoReg, c(nuisance$pseudoReg, list(pseudo = tau_tmp[k != i], X = confounders[k != i, sminus, drop = FALSE])))
        taus <- predict(taus_tmp, newX = confounders[k == i, sminus, drop = FALSE])
        EXs_Xsminus_tmp <- do.call(pseudoReg, c(nuisance$pseudoReg, list(pseudo = confounders[k != i, s], X = confounders[k != i, sminus, drop = FALSE])))
        EXs_Xsminus <- predict(EXs_Xsminus_tmp, newX = confounders[k == i, sminus, drop = FALSE])
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
        EXs_Xsminus_tmp <- do.call(pseudoReg, c(nuisance$pseudoReg, list(pseudo = confounders[k != i, s], X = confounders[k != i, sminus, drop = FALSE])))
        EXs_Xsminus <- predict(EXs_Xsminus_tmp, newX = confounders[k == i, sminus, drop = FALSE])
      }
      
      eif_taup_uncen <- phi[k == i, ]
      taup <- mean(eif_taup_uncen)
      
      eif_psi_t_uncen <- (phi[k == i, ] - taus)*(xs_numeric[k == i] - EXs_Xsminus)
      psi_t <- mean(eif_psi_t_uncen)
      
      eif_phi_t_uncen <- (xs_numeric[k == i] - EXs_Xsminus)^2
      phi_t <- mean(eif_phi_t_uncen)
      
      vals <- list(nk = sum(k == i), eif_psi_t_uncen = eif_psi_t_uncen,
                   eif_taup_uncen = eif_taup_uncen, eif_phi_t_uncen = eif_phi_t_uncen,
                   psi_t = psi_t, phi_t = phi_t, taup = taup)
      vals
    })
    
    psi_t <- sum(sapply(vals, function(x) x$nk/N * x$psi_t))
    taup <- sum(sapply(vals, function(x) x$nk/N * x$taup))
    phi_t <- sum(sapply(vals, function(x) x$nk/N * x$phi_t))
    omega_t <- psi_t/phi_t
    
    var_psi_t <- sum(sapply(vals, function(x) x$nk/N * mean((x$eif_psi_t_uncen - psi_t)^2)))/N
    var_phi_t <- sum(sapply(vals, function(x) x$nk/N * mean((x$eif_phi_t_uncen - phi_t)^2)))/N
    var_taup <- sum(sapply(vals, function(x) x$nk/N * mean((x$eif_taup_uncen - taup)^2)))/N
    var_omega_t <- sum(sapply(vals, function(x) x$nk/N * mean(((x$eif_psi_t_uncen - psi_t - omega_t*(x$eif_phi_t - phi_t))/phi_t)^2)))/N
  }
  
  list(psi_t = psi_t,
       se_Psi_t = sqrt(var_psi_t),
       pval_psi_t = 2*(1 - pnorm(abs(psi_t)/sqrt(var_psi_t))),
       phi_t = phi_t,
       se_phi_t = sqrt(var_phi_t),
       omega_t = omega_t,
       se_omega_t = sqrt(var_omega_t),
       pval_omega_t = 2*(1 - pnorm(abs(omega_t)/sqrt(var_omega_t))),
       taup = taup,
       se_taup = sqrt(var_taup),
       treatment_diff = paste0(rev(treatLabs), collapse = " - "))
  
}


