# This file contains the functions to extract results from the simulations

# This function takes in unput the stan object
# and return the posterior sample (beta + pi)

ExtractPost <- function(stan.obj, 
                        comm_seas = FALSE,  comm_trend = FALSE, comm_AR = FALSE, AR = FALSE,
                        K = 2, last.regime = c("N", "G", "rW")){
  beta <- list()
  for(k in 1:K){
    label1L <- paste0("betaL", k, "_new")
    label2L <- paste0("betaL", k)
    beta[[k]] <- extract(stan.obj, pars = label1L, permuted = FALSE)
    
    # Get the first regime betas
    if(k==1){
      beta[[k]] <- do.call(rbind, lapply(1:4, function(x) beta[[k]][,x,]))
    } else {     # Get the other regimes commn parts
      
      if(AR){
        if(comm_AR & comm_seas & comm_trend){
          beta[[k]] <-  do.call(c, lapply(1:4, function(x) beta[[k]][,x,]))
        } else if(comm_seas & !comm_trend & !comm_AR){
          beta[[k]] <- do.call(rbind, lapply(1:4, function(x) beta[[k]][,x,]))
        } else if(!comm_seas & !comm_trend & !comm_AR){
          beta[[k]] <- do.call(rbind, lapply(1:4, function(x) beta[[k]][,x,]))
        }    
      } else {
        if(comm_seas & comm_trend){
          beta[[k]] <-  do.call(c, lapply(1:4, function(x) beta[[k]][,x,]))
        } else if(comm_seas & !comm_trend){
          beta[[k]] <- do.call(rbind, lapply(1:4, function(x) beta[[k]][,x,]))
        } else if(!comm_seas & !comm_trend){
          beta[[k]] <- do.call(rbind, lapply(1:4, function(x) beta[[k]][,x,]))
        } 
      }
    }  
    
    if(k > 1){ # Padd common part with seasonality
      if(AR){
        if(comm_AR & comm_seas & comm_trend){
          beta[[k]] <- cbind(beta[[k]],  beta[[1]][, 2:ncol(beta[[1]])])
        } else if(comm_seas & !comm_trend & !comm_AR){
          beta[[k]] <- cbind(beta[[k]],  beta[[1]][,(ncol(beta[[k]])+1):ncol(beta[[1]])])
        }  
      } else {
        if(comm_seas & comm_trend){
          beta[[k]] <- cbind(beta[[k]],  beta[[1]][, 2:ncol(beta[[1]])])
        } else if(comm_seas & !comm_trend){
          beta[[k]] <- cbind(beta[[k]],  beta[[1]][,(ncol(beta[[k]])+1):ncol(beta[[1]])])
        } 
      }
    }  
    names(beta)[k] <- label2L
  }
  
  for(k in 1:K){
    label1S <- paste0("betaS", k, "_new")
    label2S <- paste0("betaS", k)
    beta[[k+K]] <- extract(stan.obj, pars = label1S, permuted = FALSE)
    beta[[k+K]] <- do.call(c, lapply(1:4, function(x) beta[[k+K]][,x,]))
    names(beta)[k + K] <- label2S
  }
  
  if(last.regime == "rW"){
    label1Sh <- paste0("betaSh", K)
    beta[[2*K + 1]] <- extract(stan.obj, pars = label1Sh, permuted = FALSE)
    beta[[2*K + 1]] <- do.call(c, lapply(1:4, function(x) beta[[2*K + 1]][,x,]))
    names(beta)[2*K + 1] <- label1Sh
  }
  
  aux_subd <- matrix(0, nrow = K, ncol = K)
  aux_subd[row(aux_subd) + 1 == col(aux_subd)] <- seq(1, 2 * (K-1),2)
  aux_subd[row(aux_subd)  == col(aux_subd) + 1] <- seq(2, 2 * (K-1),2)
  
  
  count <- length(beta) + 1  
  for(j in 1:(2*(K-1))){
    idx <- which(aux_subd == j, arr.ind = TRUE)
    label1Tp <- paste0("betaTp", idx[,1] - 1, idx[,2] - 1, "_new")
    beta[[count]] <- extract(stan.obj, pars = label1Tp, permuted = FALSE)
    beta[[count]] <- do.call(c, lapply(1:4, function(x) beta[[count]][,x,]))
    names(beta)[count] <- label1Tp
    count <- count + 1
  }
  
  pit <- extract(stan.obj, pars = c("pit"), permuted = FALSE)
  pit <- do.call(rbind, lapply(1:4, function(x) pit[,x,]))
  n <- ncol(pit)/K
  pi <- array(NA, dim = c(4000, n, K))
  for(j in 1:K){
    pi[,,j] <- pit[,(n*(j-1) + 1) :(n*j)] 
  }
  
  return(list(beta = beta, pi = pi)) 
}  


# This function gets the posterior and the data and provide the 
# posterior model paramters (mu, sigma, xi)

extract_etaDist <- function(Posterior, data, AR = FALSE, 
                            K = 2, last.regime = c("N", "G", "rW")){
  etaL <- list()
  etaS <- list()
  
  if(AR){
    fitF <- gam(y ~ MaxT_Lag1c + GST_c + s(sc_doy, bs = "bs", k=5), data = data, fit = FALSE)
  } else {
    fitF <- gam(y ~ GST_c + s(sc_doy, bs = "bs", k=5), data = data, fit = FALSE)
  }
  
  
  for(k in 1:K){
    labelL <- paste0("betaL", k)
    etaL[[k]] <- fitF$X %*% t(eval(parse(text=paste0("Posterior$beta$", labelL))))
    labelS <- paste0("betaS", k)
    etaS[[k]] <- fitF$X[,1] %*% t(eval(parse(text=paste0("Posterior$beta$", labelS))))
  }
  
  if(last.regime == "rW"){
    label1Sh <- paste0("betaSh", K)
    etaSh <- fitF$X[,1] %*% t(eval(parse(text=paste0("Posterior$beta$", label1Sh))))
  }
  
  if(last.regime == "rW"){
    return(list(etaL = etaL, 
                etaS = etaS,
                etaSh = etaSh))
  } else {
    return(list(etaL = etaL, 
                etaS = etaS))
  } 
}


# This function gets the posterior and the data and provide 
# the posterior estimates for the filtered state probabilities 
# according to algorithm 1 (Supplementary Material)

extract_ProbSummary <-  function(Posterior, dataTrain, dataTest, AR = FALSE, 
                                 K = 2, last.regime = c("N", "G", "rW")){
  
  
  TPM <- ExtractTPM(beta = Posterior$beta, dataTrain, K = K, model.foo.Tp = NULL, vtp = FALSE, last.regime = last.regime)$delta
  idx <- seq(1, nrow(dataTrain), 153)
  
  Prob_May1 <-  apply(exp(Posterior$pi),c(1,3), mean)
  
  n_p <- 2 + 2 + (K-2) * 3
  eta <- extract_etaDist(Posterior, data = dataTest, AR = AR, 
                         K = K, last.regime = last.regime)
  etaL <- eta$etaL
  etaS <- eta$etaS
  
  f_ys_pred <- list()
  for(j in 1:(K-1)){
    f_ys_pred[[j]] <- dnorm(dataTest$y, etaL[[j]], exp(etaS[[j]]), log = T)
  }
  
  if(last.regime == "N"){
    f_ys_pred[[K]] <- dnorm(dataTest$y, etaL[[K]], exp(etaS[[K]]), log = T)
  } else if (last.regime == "G"){
    f_ys_pred[[K]] <- devd(dataTest$y, loc = etaL[[K]], scale = exp(etaS[[K]]), shape = 0, log = T)
  } else if (last.regime == "rW" ){
    etaSh <- eta$etaSh
    f_ys_pred[[K]] <- devd(dataTest$y, loc = etaL[[K]], scale = exp(etaS[[K]]), shape = -exp(etaSh), log = T) #Add -exp()
  } 
  
  
  pit <- array(0, c(4000,nrow(dataTest),K))
  p_pred <- matrix(0, nrow(dataTest),K)
  pit[,1,] <-  log(Prob_May1)
  
  
  p <- matrix(NA, 4000, 2*(K-1)+K)
  f <- c()
  
  for(t in 2:nrow(dataTest)){
    p[,1] = log(TPM[,1,1]) + pit[, t-1, 1] + f_ys_pred[[1]][t,]
    p[,2] = log(TPM[,2,1]) + pit[, t-1, 2]  + f_ys_pred[[1]][t,]
    
    if(K>2){
      p[, 3] = log(TPM[,1,2]) + pit[, t - 1, 1] + f_ys_pred[[2]][t,]
      p[, 4] = log(TPM[,2,2]) + pit[, t - 1, 2] + f_ys_pred[[2]][t,]
      p[, 5] = log(TPM[,3,2]) + pit[, t - 1, 3] + f_ys_pred[[2]][t,]
    }
    if(K>3){
      count <- 0
      for(k in 3:(K-1)){
        p[, k * 2 + count] = log(TPM[,k-1,k]) + pit[, t - 1, k-1] + f_ys_pred[[k]][t,]
        p[, k * 2 + 1 + count] = log(TPM[,k,k]) + pit[, t - 1, k] + f_ys_pred[[k]][t,]
        p[, k * 2 + 2 + count] = log(TPM[, k+1,k]) + pit[, t - 1, k+1] + f_ys_pred[[k]][t,]
        count <- count + 1
      }
    }   
    
    p[, n_p - 1] = log(TPM[,K-1, K]) + pit[, t - 1, K - 1] + f_ys_pred[[K]][t,]
    p[, n_p] = log(TPM[,K, K]) + pit[, t - 1, K] + f_ys_pred[[K]][t,]
    
    f <- apply(p, 1, function(x) log(sum(exp(x))))  
    
    pit[,t,1]  = apply(p[, 1:2], 1, function(x) log(sum(exp(x)))) - f
    
    if(K>2){
      for(j in 2:(K-1)){ 
        pit[,t,j] =  apply(p[, ((j-1)*3): ((j-1)*3+2)], 1, function(x) log(sum(exp(x)))) - f
      }
    }
    pit[,t,K] =  apply(p[, (n_p-1): (n_p)], 1, function(x) log(sum(exp(x)))) - f 
  } 
  
  p_pred <- apply(exp(pit), c(2,3), mean)
  return(list(pit = pit, p_pred = p_pred))
}


# This function gets the posterior and the data and return 
# the regime mean estimates and credible intervals of distribution

extract_PredSummary <-  function(Posterior, data, AR = FALSE,
                                 K = 2, last.regime = c("N", "G", "rW")){
  
  eta <- extract_etaDist(Posterior, data = data, AR = AR,
                         K = K, last.regime = last.regime)
  etaL <- eta$etaL
  etaS <- eta$etaS
  if(last.regime == "rW"){
    etaSh  <- eta$etaSh
  }
  
  mean_reg <- array(NA, dim = c( nrow(data), 4000, K))
  mean_reg_quantiles <- array(NA, dim=c(5, nrow(data), K))
  
  for(j in 1:(K-1)){
    mean_reg[,,j] <- etaL[[j]]
    mean_reg_quantiles[,,j] <- apply(mean_reg[,,j],1, quantile, p = c(0.05, 0.25, 0.5, 0.75, 0.95))
  }
  
  if(last.regime == "N"){
    mean_reg[,,K] <- etaL[[K]]
  }
  
  if(last.regime == "G"){
    mean_reg[,,K] <-  etaL[[K]] + exp(etaS[[K]]) * 0.5772
  }
  
  if(last.regime == "rW"){
    mean_reg[,,K] <-  etaL[[K]] + exp(etaS[[K]]) * (gamma(1 + exp(etaSh)) - 1)/ (-exp(etaSh)) 
  }  
  
  mean_reg_quantiles[,,K] <- apply(mean_reg[,,K],1, quantile, p = c(0.05, 0.25, 0.5, 0.75, 0.95))
  
  return(list(mean_reg = mean_reg,
              mean_reg_quantiles = mean_reg_quantiles))  
  
}


# This function gets the posterior and the data and return 
# the regime quantile estimates of distribution

extract_Quantiles <-  function(Posterior, data, AR = FALSE,
                               K = 2, last.regime = c("N", "G", "rW"), quantiles = NULL){
  
  eta <- extract_etaDist(Posterior, data = data, AR = AR,
                         K = K, last.regime = last.regime)
  etaL <- eta$etaL
  etaS <- eta$etaS
  quantile_summary <- list()
  quantile_reg_summary <- matrix(0,  nrow(data), K)
  for(i in 1:length(quantiles)){
    quantile_reg <- array(NA, dim = c( nrow(data), 4000, K))
    for(j in 1:K){
      quantile_reg[,,j] <- etaL[[j]] + exp(etaS[[j]]) * qnorm(quantiles[i])#fitF$X %*% t(eval(parse(text=paste0("Posterior$beta$", labelL)))) 
      quantile_reg_summary[,j] <- apply(quantile_reg[,,j],1, median)
    }
    quantile_summary[[i]] <- quantile_reg_summary
  }
  
  return(quantile_summary)  
  
}

# This function gets the beta posterior sample and the data and provide 
# the multi-logit transition probabilities (nu)

extract_etaTp <- function(beta, dataTrain, 
                          K = 2, last.regime = c("N", "G", "rW")){
  etaTp <- list()
  aux_subd <- matrix(0, nrow = K, ncol = K)
  aux_subd[row(aux_subd) + 1 == col(aux_subd)] <- seq(1, 2 * (K-1),2)
  aux_subd[row(aux_subd)  == col(aux_subd) + 1] <- seq(2, 2 * (K-1),2)
  
  if(last.regime == "N" | last.regime == "G"){
    count <- K * 2 + 1
  } else {
    count <- (K * 2 + 1) + 1
  }
  
  for(j in 1:(2*(K-1))){
    etaTp[[j]] <- t(beta[[count]])
    count <- count + 1
  }
  return(etaTp)
}


# This function gets the beta posterior sample 
# and return the transition probabilities (delta)

ExtractTPM <- function(beta, data, K = 2, model.foo.Tp = NULL, vtp = FALSE, 
                       last.regime = c("N", "G", "rW")){
  
  aux_subd <- matrix(0, nrow = K, ncol = K)
  aux_subd[row(aux_subd) + 1 == col(aux_subd)] <- seq(1, 2 * (K-1),2)
  aux_subd[row(aux_subd)  == col(aux_subd) + 1] <- seq(2, 2 * (K-1),2)
  
  etaTp <- extract_etaTp(beta, data, 
                         K = K, last.regime = last.regime)
  
  delta <- array(0, c(4000, K,K))
  delta_sum <- array(0, c(K,K))
  
  idx <- which(aux_subd == 1, arr.ind = TRUE)
  delta[,idx[,1],idx[,2]] = invlogit(etaTp[[1]])
  delta_sum[idx[,1],idx[,2]] = invlogit(apply((etaTp[[1]]), 1, median))
  
  if(K>2){
    for(j in 2: (2 + 2*(K-2)-1)){
      idx <- which(aux_subd == j, arr.ind = TRUE)
      if((j %% 2) == 0){
        delta[,idx[,1],idx[,2]] = exp(etaTp[[j]])/(1+exp(etaTp[[j]]) + exp(etaTp[[j+1]]))
        delta_sum[idx[,1],idx[,2]] = apply(exp(etaTp[[j]])/(1+exp(etaTp[[j]]) + exp(etaTp[[j+1]])), 1, median)
      } else {
        delta[,idx[,1],idx[,2]] = exp(etaTp[[j]])/(1+exp(etaTp[[j]]) + exp(etaTp[[j-1]]))
        delta_sum[idx[,1],idx[,2]] = apply(exp(etaTp[[j]])/(1+exp(etaTp[[j]]) + exp(etaTp[[j-1]])), 1, median)
      }
    }
  }
  
  idx <- which(aux_subd == 2 + 2*(K-2), arr.ind = TRUE)
  delta[,idx[,1],idx[,2]] = invlogit(etaTp[[2 + 2*(K-2)]])
  delta_sum[idx[,1],idx[,2]] = invlogit(apply((etaTp[[2 + 2*(K-2)]]), 1, median))
  
  for(j in 1:K){
    delta[,j,j] <- 1 - apply(delta[,j,], 1,sum)
    delta_sum[j,j] <- 1 - sum(delta_sum[j,])
  }
  
  return(list(delta = delta,
              delta_sum = delta_sum))
}


# This function gets the posterior and the data and return the forecasted 
# in-sample temperatures and the forecasted probabilities and states

predict_temp <-  function(Posterior, dataTrain, 
                          comm_seas = comm_seas, comm_trend = comm_trend, comm_AR = comm_AR, AR = AR,
                          Nsim = 100, 
                          K = 2, last.regime = c("N", "G", "rW"), refilter = FALSE, 
                          proj = FALSE){
  
  Nsim <- 4000
  n_p <- 2 + 2 + (K-2) * 3
  
  idx_year <- ((unique(dataTrain$Year) - 1995) * 153 + 1) : ((unique(dataTrain$Year) - 1995 + 1) * 153)   
  
  hat_prob_in <- exp(Posterior$pi[,idx_year,, drop = FALSE])
  hat_regime <- matrix(0, nrow(dataTrain), Nsim) 
  y_pred <- matrix(NA, nrow = nrow(dataTrain), ncol = Nsim) 
  
  # Refiltering probability
  ProbFor <- array(NA, dim = c(nrow(dataTrain), K, Nsim))
  
  p <- rep(NA, 2*(K-1)+K)
  f <- c()
  f_ys_pred <- list()
  
  eta <- extract_etaDist(Posterior, data = dataTrain, 
                         K = K, last.regime = last.regime, AR)
  etaL <- eta$etaL
  etaS <- eta$etaS
  
  if(last.regime == "rW"){
    etaSh <- eta$etaSh
  }
  
  for(s in 1: Nsim){
    for(t in 1:nrow(dataTrain)){
      ProbFor[t,,s] <- hat_prob_in[s,t,] #t(TPM)%*%hat_prob_in
      hat_regime[t,s] <- sample(1:K,1, prob=ProbFor[t,,s]) 
      if(hat_regime[t,s] != K ){
        y_pred[t,s] <- rnorm(1,etaL[[hat_regime[t,s]]][t,s],exp(etaS[[hat_regime[t,s]]][t,s]))
      }  
      if(hat_regime[t,s] == K){
        if(last.regime == "N"){
          y_pred[t,s] <- rnorm(1,etaL[[K]][t,s],exp(etaS[[K]][t,s]))
        } else if (last.regime == "G"){
          y_pred[t,s] <- revd(1, loc = etaL[[K]][t,s],
                              scale = exp(etaS[[K]][t,s]), shape = 0, type = "GEV")  
        } else if (last.regime == "rW"){
          y_pred[t,s] <- revd(1, loc = etaL[[K]][t,s], 
                              scale = exp(etaS[[K]][t,s]), shape = -exp(etaSh[t,s]), type = "GEV")  
        }
      }
    }
    print(s)
  }  
  
  return(list(y_pred = y_pred, 
              ProbFor = ProbFor, 
              hat_regime = hat_regime))
} 


# This function gets the posterior and the data and return the forecasted 
# out-of-sample temperatures and the forecasted probabilities and states

forecast_temp <-  function(Posterior, dataTrain, dataTest, 
                           comm_seas = comm_seas, comm_trend = comm_trend, 
                           comm_AR = comm_AR, AR = AR, Nsim = 100, K = 2, 
                           last.regime = c("N", "G", "rW"), refilter = FALSE, 
                           proj = FALSE){
  
  
  Nsim <- 4000
  n_p <- 2 + 2 + (K-2) * 3
  
  if(AR){
    fitF <- gam(y ~ MaxT_Lag1c + GST_c + s(sc_doy, bs = "bs", k=5), data = dataTest, fit = FALSE)
    Xl <- fitF$X
  } else {
    fitF <- gam(y ~ GST_c + s(sc_doy, bs = "bs", k=5), data = dataTest, fit = FALSE)
    Xl <- fitF$X
  }
  
  idx <- seq(1, nrow(dataTrain), 153)
  hat_prob_in <- colMeans(apply(exp(Posterior$pi),c(1,3), mean))
  
  hat_regime <- matrix(0, 153, Nsim) 
  hat_regime[1,] <- do.call(c, lapply(1:Nsim, function(x) sample(1:K,1, prob=hat_prob_in))) 
  
  y_pred <- matrix(NA, nrow = 153, ncol = Nsim) 
  y_pred[1,] <- dataTest$y[1]#[start]
  
  # Refiltering probability
  ProbFor <- array(NA, dim = c(153, K, Nsim))
  ProbFor[1,,] <- hat_prob_in
  
  p <- rep(NA, 2*(K-1)+K)
  f <- c()
  f_ys_pred <- list()
  
  TPM <- ExtractTPM(beta = Posterior$beta, dataTrain, K = K, model.foo.Tp = NULL, vtp = FALSE, last.regime = last.regime)$delta#_sum
  
  eta <- extract_etaDist(Posterior, data = dataTest, 
                         K = K, last.regime = last.regime, AR)
  etaL <- eta$etaL
  etaS <- eta$etaS
  
  if(last.regime == "rW"){
    etaSh <- eta$etaSh
  }
  
  for(s in 1: Nsim){
    for(t in 2:153){
      if(t == 2){
        ProbFor[t,,s] <- t(TPM[s,,]) %*% hat_prob_in
      } else {
        ProbFor[t,,s] <- t(TPM[s,,]) %*% ProbFor[t-1,,s]
      }
      
      hat_regime[t,s] <- sample(1:K,1, prob=ProbFor[t,,s]) 
      
      if(!AR){
        if(hat_regime[t,s] != K ){
          y_pred[t,s] <- rnorm(1,etaL[[hat_regime[t,s]]][t,s],exp(etaS[[hat_regime[t,s]]][t,s]))
        }  
      } else if(AR){
        if(hat_regime[t,s] !=K){
          y_pred[t,s] <- rnorm(1, median(Posterior$beta[[hat_regime[t,s]]][,1]) +
                                 median(Posterior$beta[[hat_regime[t,s]]][,2]) *  (y_pred[t-1, s]  - mean(dataTrain$MaxT_Lag1)) +
                                 median(Posterior$beta[[hat_regime[t,s]]][,3]) *  Xl[t,3] +
                                 sum(Xl[t,4:ncol(Xl)] * t(apply(Posterior$beta[[1]][,4:ncol(Xl)],2,median))),
                               exp(etaS_median[[hat_regime[t,s]]][t]))
        }
      }
      
      
      
      if(!AR){
        if(hat_regime[t,s] == K){
          if(last.regime == "N"){
            y_pred[t,s] <- rnorm(1,etaL[[K]][t,s],exp(etaS[[K]][t,s]))
          } else if (last.regime == "G"){
            y_pred[t,s] <- revd(1, loc = etaL[[K]][t],
                                scale = exp(etaS[[K]][t]), shape = 0, type = "GEV")  
          } else if (last.regime == "rW"){
            y_pred[t,s] <- revd(1, loc = etaL[[K]][t], scale = exp(etaS[[K]][t]), shape = -exp(etaSh[t]), type = "GEV")  
          }
        }
        
        for(j in 1:(K-1)){
          f_ys_pred[[j]] <- dnorm(y_pred[t,s], etaL[[j]][t,s],
                                  exp(etaS[[j]][t,s]), log = T)  
        }
        
        if(last.regime == "N"){
          f_ys_pred[[K]] <- dnorm(y_pred[t,s], etaL[[K]][t,s],
                                  exp(etaS[[K]][t,s]), log = T)  
        } else if (last.regime == "G"){
          f_ys_pred[[K]] <- devd(y_pred[t,s], etaL[[K]][t,s], 
                                 scale = exp(etaS[[K]][t,s]), shape = 0, type = "GEV", log = T)
        } else if (last.regime == "rW"){
          f_ys_pred[[K]] <- devd(y_pred[t,s], loc = etaL[[K]][t,s],
                                 scale = exp(etaS[[K]][t,s]), shape = -exp(etaSh[t,s]), type = "GEV", log = T)
        }        
      } else {
        if(hat_regime[t,s] == K){
          if(last.regime == "N"){
            y_pred[t,s] <- rnorm(1, median(Posterior$beta[[K]][,1]) + 
                                   median(Posterior$beta[[K]][,2]) *  (y_pred[t-1, s]  - mean(dataTrain$MaxT_Lag1)) +  
                                   median(Posterior$beta[[K]][,3]) *  Xl[t,3] + 
                                   sum(Xl[t,4:ncol(Posterior$beta[[1]])] * t(apply(Posterior$beta[[1]][,4:ncol(Posterior$beta[[1]])],2,median))),
                                 exp(etaS_median[[K]][t]))
          } else if (last.regime == "G"){
            y_pred[t,s] <- revd(1, median(Posterior$beta[[K]][,1]) + 
                                  median(Posterior$beta[[K]][,2]) *  (y_pred[t-1, s]  - mean(dataTrain$MaxT_Lag1)) +  
                                  median(Posterior$beta[[K]][,3]) *  Xl[t,3] + 
                                  sum(Xl[t,4:ncol(Posterior$beta[[1]])] * t(apply(Posterior$beta[[1]][,4:ncol(Posterior$beta[[1]])],2,median))),
                                scale = exp(etaS_median[[K]][t]), shape = 0, type = "GEV")  
          } else if (last.regime == "rW"){
            y_pred[t,s] <- revd(1, median(Posterior$beta[[K]][,1]) + 
                                  median(Posterior$beta[[K]][,2]) *  (y_pred[t-1, s]  - mean(dataTrain$MaxT_Lag1)) +  
                                  median(Posterior$beta[[K]][,3]) *  Xl[t,3] + 
                                  sum(Xl[t,4:ncol(Posterior$beta[[1]])] * t(apply(Posterior$beta[[1]][,4:ncol(Posterior$beta[[1]])],2,median))),
                                scale = exp(etaS_median[[K]][t]), shape = -exp(etaSh_median[t]), type = "GEV")  
          }
        }
        
        for(j in 1:(K-1)){
          f_ys_pred[[j]] <- dnorm(y_pred[t,s], median(Posterior$beta[[j]][,1]) + 
                                    median(Posterior$beta[[j]][,2]) *  (y_pred[t-1, s]  - mean(dataTrain$MaxT_Lag1)) +  
                                    median(Posterior$beta[[j]][,3]) *  Xl[t,3] + 
                                    sum(Xl[t,4:ncol(Posterior$beta[[1]])] * t(apply(Posterior$beta[[1]][,4:ncol(Posterior$beta[[1]])],2,median))),
                                  exp(etaS_median[[j]][t]), log = T)  
        }
        
        if(last.regime == "N"){
          f_ys_pred[[K]] <- dnorm(y_pred[t,s], median(Posterior$beta[[K]][,1]) + 
                                    median(Posterior$beta[[K]][,2]) *  (y_pred[t-1, s]  - mean(dataTrain$MaxT_Lag1)) +  
                                    median(Posterior$beta[[K]][,3]) *  Xl[t,3] + 
                                    sum(Xl[t,4:ncol(Posterior$beta[[1]])] * t(apply(Posterior$beta[[1]][,4:ncol(Posterior$beta[[1]])],2,median))),
                                  exp(etaS_median[[K]][t]), log = T)  
        } else if (last.regime == "G"){
          f_ys_pred[[K]] <- devd(y_pred[t,s], median(Posterior$beta[[K]][,1]) + 
                                   median(Posterior$beta[[K]][,2]) *  (y_pred[t-1, s]  - mean(dataTrain$MaxT_Lag1)) +  
                                   median(Posterior$beta[[K]][,3]) *  Xl[t,3] + 
                                   sum(Xl[t,4:ncol(Posterior$beta[[1]])] * t(apply(Posterior$beta[[1]][,4:ncol(Posterior$beta[[1]])],2,median))),
                                 scale = exp(etaS_median[[K]][t]), shape = 0, type = "GEV", log = T)
        } else if (last.regime == "rW"){
          f_ys_pred[[K]] <- devd(y_pred[t,s], median(Posterior$beta[[K]][,1]) + 
                                   median(Posterior$beta[[K]][,2]) *  (y_pred[t-1, s]  - mean(dataTrain$MaxT_Lag1)) +  
                                   median(Posterior$beta[[K]][,3]) *  Xl[t,3] + 
                                   sum(Xl[t,4:ncol(Posterior$beta[[1]])] * t(apply(Posterior$beta[[1]][,4:ncol(Posterior$beta[[1]])],2,median))),
                                 scale = exp(etaS_median[[K]][t]), shape = -exp(etaSh_median[t]), type = "GEV", log = T)
        }
      }
      p[1] = log(TPM[s,1,1]) + log(ProbFor[t-1,1,s])  + f_ys_pred[[1]]
      p[2] = log(TPM[s,2,1]) + log(ProbFor[t-1,2,s])  + f_ys_pred[[1]]
      if(K>2){
        p[3] = log(TPM[s,1,2]) + log(ProbFor[t - 1, 1, s]) + f_ys_pred[[2]]
        p[4] = log(TPM[s,2,2]) + log(ProbFor[t - 1, 2, s]) + f_ys_pred[[2]]
        p[5] = log(TPM[s,3,2]) + log(ProbFor[t - 1, 3, s]) + f_ys_pred[[2]]
      }
      if(K>3){
        count <- 0
        for(k in 3:(K-1)){
          p[k * 2 + count] = log(TPM[s,k-1,k]) + log(ProbFor[t - 1, k-1, s]) + f_ys_pred[[k]]
          p[k * 2 + 1 + count] = log(TPM[s,k,k]) + log(ProbFor[t - 1, k, s]) + f_ys_pred[[k]]
          p[k * 2 + 2 + count] = log(TPM[s,k+1,k]) + log(ProbFor[t - 1, k+1, s]) + f_ys_pred[[k]]
          count <- count + 1
        }
      }   
      
      p[n_p - 1] = log(TPM[s,K-1, K]) + log(ProbFor[t-1,K-1,s]) + f_ys_pred[[K]]
      p[n_p] = log(TPM[s,K, K]) + log(ProbFor[t-1,K,s]) + f_ys_pred[[K]]
      
      f <-  log(sum(exp(p)))
      
      ProbFor[t,1,s] = exp(log(sum(exp(p[1:2]))) - f)
      
      if(K>2){
        for(j in 2:(K-1)){ 
          ProbFor[t,j,s] =  exp(log(sum(exp(p[((j-1)*3): ((j-1)*3+2)]))) - f)
        }
      }
      
      ProbFor[t,K,s] = exp(log(sum(exp(p[(n_p-1): (n_p)]))) - f)
    }
    print(s)
  }  
  
  return(list(y_pred = y_pred, 
              ProbFor = ProbFor, 
              hat_regime = hat_regime))
}


# This function gets the posterior and the data and return the 
# in-sample temperature trajectories (using predict_temp)

TempIn <- function(Posterior, dataTrain, AR = AR, K = 2, 
                   last.regime = c("N", "G", "rW"), Nsim = 100,
                   comm_seas, comm_trend, comm_AR, proj, refilter = FALSE){
  yearTrain <- unique(dataTrain$Year)
  
  idx_yearTrain <- lapply(yearTrain, function(x) which(dataTrain$Year == x))
  
  compute_insample <- function(Posterior, dataTrain ,  comm_seas, comm_trend, 
                                    comm_AR, AR,  Nsim, K, last.regime,
                                    refilter, proj){
    TempTrajectories <- predict_temp(Posterior, dataTrain,  
                                     comm_seas = comm_seas, comm_trend = comm_trend,
                                     comm_AR = comm_AR, AR = AR,  Nsim = Nsim,
                                     K = K, last.regime = last.regime, refilter = FALSE, 
                                     proj = proj)
    return(list(TempTrajectories = TempTrajectories))
  }
  
  ncores <- 2 
  cl <- makePSOCKcluster(ncores)
  setDefaultCluster(cl)
  clusterExport(NULL, c("Posterior", "dataTrain", "comm_seas", "comm_trend",
                        "comm_AR", "AR", 
                        "Nsim", "K", "last.regime", "refilter", "proj", "compute_insample",
                        "predict_temp", "ExtractTPM", "extract_etaDist", "extract_etaTp",
                        "idx_yearTrain","invlogit","gam","extract_ProbSummary"), envir = environment())
  clusterEvalQ(NULL, {
    library(mgcv)  
    library(loo)
    library(extRemes)
  })
  
  
  
  out <- function(.x){
    out2 <- compute_insample(Posterior, dataTrain = dataTrain[idx_yearTrain[[.x]],], comm_seas, comm_trend, 
                                  comm_AR, AR, Nsim, K, last.regime,
                                  refilter, proj)
    return(list(gen = out2))
  }
  
  environment(out) <- .GlobalEnv
  
  res <- list()
  res <- parLapply(NULL, 1 : length(idx_yearTrain), out)
  
  stopCluster(cl)
  rm(cl)
  
  return(res)
}


# This function gets the posterior and the data and return the 
# out-of-sample temperature trajectories (using forecast_temp)

TempOut <- function(Posterior, dataTrain, dataTest, AR = AR, K = 2, 
                       last.regime = c("N", "G", "rW"), Nsim = 100,
                       comm_seas, comm_trend, comm_AR, proj, refilter = FALSE){

  yearTest <- unique(dataTest$Year)
  
  idx_yearTest <- lapply(yearTest, function(x) which(dataTest$Year == x))
  
  compute_TempTrajectories <- function(Posterior, dataTrain, dataTest, 
                                       comm_seas, comm_trend, comm_AR, AR,  
                                       Nsim, K, last.regime, refilter, proj){
    TempTrajectories <- forecast_temp(Posterior, dataTrain, dataTest, 
                                                        comm_seas = comm_seas, comm_trend = comm_trend,
                                                        comm_AR = comm_AR, AR = AR,  Nsim = Nsim,
                                                        K = K, last.regime = last.regime, refilter = FALSE, 
                                                        proj = proj)
    return(list(TempTrajectories = TempTrajectories))
  }
  
  ncores <- 2
  cl <- makePSOCKcluster(ncores)
  setDefaultCluster(cl)
  clusterExport(NULL, c("Posterior", "dataTrain", "dataTest", "comm_seas", "comm_trend",
                        "comm_AR", "AR", 
                        "Nsim", "K", "last.regime", "refilter", "proj", "compute_TempTrajectories",
                        "forecast_temp", 
                        "ExtractTPM", "extract_etaDist", "extract_etaTp",
                        "idx_yearTest","invlogit","gam","extract_ProbSummary"), envir = environment())
  clusterEvalQ(NULL, {
    library(mgcv)  
    library(loo)
    library(extRemes)
  })
  
  
  out <- function(.x){
    out2 <- compute_TempTrajectories(Posterior, dataTrain,  dataTest[idx_yearTest[[.x]],], comm_seas, comm_trend, 
                         comm_AR, AR, Nsim, K, last.regime,
                         refilter, proj)
    return(list(gen = out2))
  }
  
  environment(out) <- .GlobalEnv
  
  
  res <- list()
  res <- parLapply(NULL, 1 : length(idx_yearTest), out)
  
  stopCluster(cl)
  rm(cl)
  
  return(res)
}


# This function returns the heatwaves detected by the Markov-switching model 
# (max-prob classification for three consecutive days)

get_excedance <- function(K, pi, year_start, year_end, last.regime = c("N", "G", "rW")){
  
  year <- rep(year_start:year_end, each = 153)
  data <- data.frame(Year = year, doy = 1:153, t = 1:(153 * (year_end - year_start + 1)), 
                     month = c(rep(5,31), rep(6,30), rep(7,31), rep(8,31), rep(9,30)),
                     dom = c(1:31,1:30,1:31,1:31,1:30))
  data$date <- with(data,as.Date(paste0(year,"-", month,"-", dom)))
  
  dataHWs <-  data.frame(year = year_start:year_end)
  dataHWdays <-  data.frame(year = year_start:year_end)
  get_prob <- apply(exp(pi), c(2,3), mean)
  get_reg <- apply(get_prob, 1, which.max)
  data$Reg <- get_reg 
  data$MaxP <- get_reg == K
  data$flag_HW <- rep(0, length(year))
  HW_MaxP <-list() 
  
  count <- 1
  for(j in year_start:year_end){
    idx_year <- which(year == j)
    tmp <- rle(data$MaxP[idx_year])
    data$flag_HW[idx_year] <- with(tmp, rep(lengths * values, lengths))
    HW_MaxP[[count]] <- data[idx_year,][which(data$flag_HW[idx_year] >= 3),] 
    count <- count + 1
  }
  
  HWprob <- do.call("rbind", HW_MaxP)
  HWprob <- as.data.frame(HWprob)
  
  ExcessYear <- tapply(data$flag_HW, data$Year, table)
  HWsYear <- lapply(1: length(ExcessYear), 
                    function(x) {idx_table <- which(as.numeric(names(ExcessYear[[x]]))>=3); 
                    as.numeric(ExcessYear[[x]])[idx_table] / (as.numeric(names(ExcessYear[[x]])))[idx_table]} )
  dataHWs$Freq <- sapply(HWsYear, sum)
  dataHWdays$Freq <- sapply(1: length(ExcessYear), function(x) {idx_table <- which(as.numeric(names(ExcessYear[[x]]))>=3); sum(as.numeric(ExcessYear[[x]])[idx_table]) })
  colnames(dataHWdays)[2] <- colnames(dataHWs)[2] <- "Freq" 
  
  
  HWs <- data[,7:9]
  colnames(HWs) <- c(paste0("Reg_", paste0(rep("N", K -1), collapse = ""), last.regime),
                     "excedance", "exc") 
  return(list(HWs =  HWs, NHWs = as.data.frame(dataHWs), DHWs = as.data.frame(dataHWdays))) 
}


# This function gets the simulated trajectories and the empirical threshold
# and returns the simulated number of heat waves and days in heat wave 
# (obtained as exceedence of the empirical based threshold)

getForeHWdist <- function(trajectories, dataHWr, 
                          year_start, year_end, 
                          AR, model_names,
                          K, last.regime,
                          site, Nsim, year_thresh){ 
  
  year_seq <- year_start:year_end
  dataSimulation <- list()
  for(k in 1:length(trajectories)){
    dataSimulation[[k]] <- list()
    for(i in 1:length(year_seq)){
      heatwaveRdata <- dataHWr[dataHWr$years ==  year_thresh, ]
      dataSimulation[[k]][[i]] <- data.frame(T = 1:(2*Nsim), 
                                             year = year_seq[i], 
                                             type = model_names[k],
                                             AR = AR,
                                             K = label_K[k],
                                             last.regime = last.regime[k], 
                                             site = site,
                                             measure = rep(c("HWd", "HWn"), each = Nsim),
                                             Freq = 0)
      for(j in 1:Nsim){
        aux1 <- data.frame("doy" = 1:153)
        aux2 <- trajectories[[k]][[i]][,j] >    heatwaveRdata$tresh  
        tmp <- rle(aux2)
        aux3 <- with(tmp, rep(lengths * values, lengths))
        aux1$getdays <- aux3
        aux1$flag_excedance <- aux1$getdays>=3
        dataSimulation[[k]][[i]]$Freq[j] <- sum(aux1$flag_excedance)
        dataSimulation[[k]][[i]]$Freq[Nsim + j] <- length(tmp$lengths[tmp$lengths>=3 & tmp$values == TRUE])  
      }
    }
  }
  return(dataSimulation)
}


# This function gets the simulated trajectories and regime state membership 
# and returns the simulated number of heat waves and days in heat wave 
# (obtained when temperatures are in the last state for at least three days)

getForeHWdist2 <- function(trajectories, states, year_start, year_end, 
                           AR, model_names, K, last.regime,
                           site, Nsim, year_thresh){ 
  
  year_seq <- year_start:year_end
  dataSimulation <- list()
  for(k in 1:length(trajectories)){
    dataSimulation[[k]] <- list()
    for(i in 1:length(year_seq)){
      dataSimulation[[k]][[i]] <- data.frame(T = 1:(2*Nsim), 
                                             year = year_seq[i], 
                                             type = model_names[k],
                                             AR = AR,
                                             K = label_K[k],
                                             last.regime = last.regime[k], 
                                             site = site,
                                             measure = rep(c("HWd", "HWn"), each = Nsim),
                                             Freq = 0)
      for(j in 1:Nsim){
        aux1 <- data.frame("doy" = 1:153)
        aux2 <- ( states[[k]][[i]][,j] == K[k])  
        tmp <- rle(aux2)
        aux3 <- with(tmp, rep(lengths * values, lengths))
        aux1$getdays <- aux3
        aux1$flag_excedance <- aux1$getdays>=3
        dataSimulation[[k]][[i]]$Freq[j] <- sum(aux1$flag_excedance)
        dataSimulation[[k]][[i]]$Freq[Nsim + j] <- length(tmp$lengths[tmp$lengths>=3 & tmp$values == TRUE])  
      }
    }
  }
  return(dataSimulation)
}


# This function assigns a progressive number to yearly heat waves and
# it is used to obtain some plots

numberingHW <- function(obj){
  if(nrow(obj) !=0 ){
    idx_diff <- diff(obj$doy)
    progNumbHW <- rep(0, length(idx_diff))
    progNumbHW[1] <- 1
    count <- 2
    for(j in 1:length(idx_diff)){
      if(idx_diff[j]==1){
        progNumbHW[j + 1] <- count - 1
      } else {
        progNumbHW[j + 1] <- count
        count <- count  + 1
      }
    }
    obj$progNumbHW <-   progNumbHW
    return(obj)
  } else {
    
    obj$progNumbHW <- numeric(0)
    return(obj)    
  }
}

