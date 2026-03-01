#This code generates the stan file according to the specified model formula

# Call the fit_stan function to run the model
model_fitting <- function(data, modfoo, trainStart, trainEnd, K = 2, 
                          regimes = "N-N", iter = 2000, initialise_param = NULL,
                          piStart = c(0.5, 0.5), AR = FALSE, comm_trend = FALSE, comm_AR = FALSE, comm_seas = FALSE) {
  
  datiTrain <- data[data$Year <= trainEnd,]
  datiTest <- data[data$Year > trainEnd,]

  fit <- fit_stan(data = datiTrain, data_pred = datiTest, regimes = regimes, 
                  comm_trend = comm_trend, comm_AR = comm_AR, comm_seas = comm_seas, AR = AR, 
                  mod_foo = modfoo, pi00init = piStart, K = K,
                  iter = iter, init_par = initialise_param, ncores=4,
                  adapt_delta = 0.99, max_treedepth = 10)


  return(fit)
}

# Creates a list with model formulas
get_foo <- function(formula, K = 2){
  out <- list()
  
  aaa <- expand.grid(c(0:(K-1)), c(0:(K-1)))
  aaa <- aaa[abs(apply(aaa, 1, diff))==1,]
  aaa <- apply(format(aaa), 1, paste, collapse = "")
  
  # Working with distributions of regimes 
  Tp_length <- length(aaa)
  
  dist <- rep(NA, K)

  dist <- sapply(formula[1:K], function(x) rhs.vars(as.formula(gsub('"', "'", deparse(x, width.cutoff = 100)))))
  nf <- c(rep(2,K-1),ifelse(tail(dist, 1)=="rWeibull", 3, 2))
  out <- lapply(1:K, function(x) as.list(strsplit(gsub('"', "'",formula[[x]]), split = "[|]")[[2]][1:nf[x]], ))
  out <- lapply(out, function(x) lapply(x, function(y) as.formula(y)))
  out <- lapply(out, function(x) lapply(x, function(y) as.formula(paste("y", paste(as.character(rhs(y)), collapse = "")))))
  
  out2 <- c(out, list(formula("y~1"))[rep(1,Tp_length)])

  return(list(out2 = out2, regimes = dist)) 
}

# Processes the model formulas and returns the related objects needed for estimate
process_modfoo <- function(obj, data, K = 2){
  Dim <- length(obj)
  mat <- vector("list", Dim)   
  S_jkl <- vector("list", Dim)           # Penalty matrices for the j-th regime, the k -th linear predictor and the l-th smooth term
  X_jk <- vector("list", Dim)            # Model matrix for the j-th regime and the k-th linear predictor 
  lin_eff <- vector("list", Dim)         # linear effects for the j-th regime and the k-th linear predictor
  smooth_eff <- vector("list", Dim)      # smooth effects for the j-th regime and the k-th linear predictor
  vs_eff <- vector("list", Dim)          # variables for smooth effects  for the j-th regime and the k-th linear predictor
  idx_lin_eff <- vector("list", Dim)     #
  idx_smooth_eff <- vector("list", Dim)
  n_eff <- vector("list", Dim)
  n_eff_lin <- vector("list", Dim)
  n_eff_smooth <- vector("list", Dim)
  
  for(i in 1 : Dim){
    Dim2 <- length(list(obj[[i]]))
    mat[[i]] <- vector("list", Dim2)
    lin_eff[[i]] <- vector("list", Dim2)
    smooth_eff[[i]] <- vector("list", Dim2) 
    vs_eff[[i]] <- vector("list", Dim2)
    idx_lin_eff[[i]] <- vector("list", Dim2)
    idx_smooth_eff[[i]] <- vector("list", Dim2)
    X_jk[[i]] <- vector("list", Dim2)
    S_jkl[[i]] <- vector("list", Dim2)
    n_eff[[i]] <- vector("list", Dim2)
    n_eff_lin[[i]] <- vector("list", Dim2)
    n_eff_smooth[[i]] <- vector("list", Dim2)
    
    if(i <= K){
      
      for(j in 1:length(obj[[i]])){
        
        if(obj[[i]][[j]] == formula(y ~ 1)){
          X_jk[[i]][[j]] <- matrix(1, nrow(data), 1)
          S_jkl[[i]][[j]] <- matrix(1,1,1)
          n_eff[[i]][[j]] <- 1
          n_eff_lin[[i]][[j]] <- 0
          n_eff_smooth[[i]][[j]] <- 0
          mat[[i]][[j]] <- matrix(0,1,3) 
          lin_eff[[i]][[j]] <-list(NULL)
          smooth_eff[[i]][[j]] <- list(NULL)
          idx_smooth_eff[[i]][[j]] <- list(NULL)
          idx_lin_eff[[i]][[j]] <- list(NULL)
          vs_eff[[i]][[j]] <- list(NULL)
        } else {
          b <- gam(obj[[i]][[j]], data = data, fit =FALSE)
          flag_lin <- identical(attr(b$pterms, "term.labels"), character(0))      # FALSE means there are linear effect
          flag_smooth <- is.null(unlist(lapply(b$smooth, function(x) x$label)))   # FALSE means there are smooth effects
          if(flag_lin & (!flag_smooth)) {  # Only smooth effects in linear predictor
            lin_eff[[i]][[j]] <-list(NULL)
            idx_lin_eff[[i]][[j]] <- list(NULL)
            smooth_eff[[i]][[j]] <- unlist(lapply(b$smooth, function(x) x$label))
            vs_eff[[i]][[j]] <- unlist(lapply(b$smooth, function(x) x$vn))
            
            idx2 <- rep(NA, length(vs_eff[[i]][[j]]))
            for(s in 1:length(vs_eff[[i]][[j]])){
              idx2[s] <- grepl(vs_eff[[i]][[j]][s], smooth_eff[[i]][[j]][s])
            }
            
            idx3 <- rep(FALSE, length( smooth_eff[[i]][[j]])) 
            idx3[!(1:length(idx3)) %in% idx_smooth_eff[[i]][[j]] ] <- idx2
            
            idx_smooth_eff[[i]][[j]] <- which(idx3)
            
            
            X_jk[[i]][[j]] <- b$X
            mat_S <- matrix(0, ncol(b$X)-1, ncol(b$X)-1)
            mat[[i]][[j]] <- matrix(0,length(b$smooth), 3)
            mat[[i]][[j]][,1] <- c(unlist(lapply(b$smooth, function(x) x$first.para)))
            mat[[i]][[j]][,2] <- c(unlist(lapply(b$smooth, function(x) x$last.para)))
            mat[[i]][[j]][,3] <- 1:length(unlist(lapply(b$smooth, function(x) x$last.para)))
            mat_S <- as.matrix(do.call("bdiag", lapply(1:nrow(mat[[i]][[j]]), function(x) b$S[[x]])))
            S_jkl[[i]][[j]] <-  mat_S
            n_eff[[i]][[j]] <- 1 + length(smooth_eff[[i]][[j]])
            n_eff_lin[[i]][[j]] <- 0
            n_eff_smooth[[i]][[j]] <- length(smooth_eff[[i]][[j]])
            mat[[i]][[j]] <- rbind(matrix(0,1,3), mat[[i]][[j]])
            
          } else if ((!flag_lin) & flag_smooth){                                  # Only linear effects in linear predictor 
            lin_eff[[i]][[j]] <- attr(b$pterms, "term.labels")
            smooth_eff[[i]][[j]] <- list(NULL)
            idx_smooth_eff[[i]][[j]] <- list(NULL)
            idx_lin_eff[[i]][[j]] <- which(rhs.vars(obj[[i]][[j]]) %in% lin_eff[[i]][[j]])
            X_jk[[i]][[j]] <- b$X
            S_jkl[[i]][[j]] <-  diag(1,ncol(b$X)-1)
            n_eff[[i]][[j]] <- 1 + length(lin_eff[[i]][[j]])
            n_eff_lin[[i]][[j]] <- length(lin_eff[[i]][[j]])
            n_eff_smooth[[i]][[j]] <- 0
            vs_eff[[i]][[j]] <- list(NULL)
            mat[[i]][[j]] <- matrix(0,length(attr(b$pterms, "term.labels")), 3)
            mat[[i]][[j]][,1] <- c(2:(length(attr(b$pterms, "term.labels") )+1))
            mat[[i]][[j]][,2] <- c(2:(length(attr(b$pterms, "term.labels") )+1))
            mat[[i]][[j]][,3] <- rep(0,length(attr(b$pterms, "term.labels") )) 
            mat[[i]][[j]] <- rbind(matrix(0,1,3), mat[[i]][[j]])
          } else {                                                                # Both linear and smooth effects in linear predictor
            lin_eff[[i]][[j]] <- attr(b$pterms, "term.labels")  
            idx_lin_eff[[i]][[j]] <- which(rhs.vars(obj[[i]][[j]]) %in% lin_eff[[i]][[j]])
            smooth_eff[[i]][[j]] <- unlist(lapply(b$smooth, function(x) x$label))
            vs_eff[[i]][[j]] <- unlist(lapply(b$smooth, function(x) x$vn))   
            idx2 <- rep(NA, length(vs_eff[[i]][[j]]))
            for(s in 1:length(vs_eff[[i]][[j]])){
              idx2[s] <- grepl(vs_eff[[i]][[j]][s], smooth_eff[[i]][[j]][s])
            }
            idx3 <- rep(FALSE, length(lin_eff[[i]][[j]]) + length( smooth_eff[[i]][[j]]))
            idx3[!(1:length(idx3)) %in% idx_lin_eff[[i]][[j]] ] <- idx2
            
            idx_smooth_eff[[i]][[j]] <- which(rhs.vars(obj[[i]][[j]]) %in% smooth_eff[[i]][[j]] | idx3)    
            
            mat[[i]][[j]] <- matrix(0, length(attr(b$pterms, "term.labels") ) + length( b$smooth),3)
            mat[[i]][[j]][,1] <- c(2:(length(attr(b$pterms, "term.labels") )+1), unlist(lapply(b$smooth, function(x) x$first.para)))               #indices for accessing 
            mat[[i]][[j]][,2] <- c(2:(length(attr(b$pterms, "term.labels") )+1), unlist(lapply(b$smooth, function(x) x$last.para)))
            mat[[i]][[j]][,3] <- c(rep(0,length(attr(b$pterms, "term.labels") ) ),  1:length(unlist(lapply(b$smooth, function(x) x$last.para))))   #
            X_jk[[i]][[j]] <- cbind(b$X[,"(Intercept)"], do.call("cbind",lapply(1:nrow(mat[[i]][[j]]), function(x) b$X[, mat[[i]][[j]][x,1]: mat[[i]][[j]][x,2], drop = FALSE])))
            mat_S <- matrix(0, ncol(b$X)-1, ncol(b$X)-1)
            mat_S <- as.matrix(do.call("bdiag", lapply(1:nrow(mat[[i]][[j]]), function(x) if( mat[[i]][[j]][x,1]== mat[[i]][[j]][x,2]){matrix(1)} else {b$S[[ mat[[i]][[j]][x,3]]]})))
            S_jkl[[i]][[j]] <-  mat_S
            n_eff[[i]][[j]] <- 1 + length(lin_eff[[i]][[j]]) + length(smooth_eff[[i]][[j]])
            n_eff_lin[[i]][[j]] <- length(lin_eff[[i]][[j]])
            n_eff_smooth[[i]][[j]] <- length(smooth_eff[[i]][[j]])
            mat[[i]][[j]] <- rbind(matrix(0,1,3), mat[[i]][[j]])
          }
        }
      }
    } else {
        if(obj[[i]] == formula(y ~ 1)){
        X_jk[[i]] <- matrix(1, nrow(data), 1)
        S_jkl[[i]] <- matrix(1,1,1)
        n_eff[[i]] <- 1
        n_eff_lin[[i]] <- 0
        n_eff_smooth[[i]] <- 0
        mat[[i]] <- matrix(0,1,3) 
        lin_eff[[i]] <-list(NULL)
        smooth_eff[[i]] <- list(NULL)
        idx_smooth_eff[[i]] <- list(NULL)
        idx_lin_eff[[i]] <- list(NULL)
        vs_eff[[i]] <- list(NULL)
      } else {
        b <- gam(obj[[i]], data = data, fit =FALSE)
        flag_lin <- identical(attr(b$pterms, "term.labels"), character(0))      # FALSE means there are linear effect
        flag_smooth <- is.null(unlist(lapply(b$smooth, function(x) x$label)))   # FALSE means there are smooth effects
        if(flag_lin & (!flag_smooth)) {  # Only smooth effects in linear predictor
          lin_eff[[i]] <-list(NULL)
          idx_lin_eff[[i]] <- list(NULL)
          smooth_eff[[i]] <- unlist(lapply(b$smooth, function(x) x$label))
          vs_eff[[i]] <- unlist(lapply(b$smooth, function(x) x$vn))
          
          idx2 <- rep(NA, length(vs_eff[[i]]))
          for(s in 1:length(vs_eff[[i]])){
            idx2[s] <- grepl(vs_eff[[i]][s], smooth_eff[[i]][s])
          }
          
          idx3 <- rep(FALSE, length( smooth_eff[[i]])) # length(lin_eff[[i]][[j]]) +
          idx3[!(1:length(idx3)) %in% idx_smooth_eff[[i]] ] <- idx2
          
          idx_smooth_eff[[i]] <- which(idx3)
          
          
          X_jk[[i]] <- b$X
          mat_S <- matrix(0, ncol(b$X)-1, ncol(b$X)-1)
          mat[[i]] <- matrix(0,length(b$smooth), 3)
          mat[[i]][,1] <- c(unlist(lapply(b$smooth, function(x) x$first.para)))
          mat[[i]][,2] <- c(unlist(lapply(b$smooth, function(x) x$last.para)))
          mat[[i]][,3] <- 1:length(unlist(lapply(b$smooth, function(x) x$last.para)))
          mat_S <- as.matrix(do.call("bdiag", lapply(1:nrow(mat[[i]]), function(x) b$S[[x]])))
          S_jkl[[i]] <-  mat_S
          n_eff[[i]] <- 1 + length(smooth_eff[[i]])
          n_eff_lin[[i]] <- 0
          n_eff_smooth[[i]] <- length(smooth_eff[[i]])
          mat[[i]] <- rbind(matrix(0,1,3), mat[[i]])
        } else if ((!flag_lin) & flag_smooth){      # Only linear effects in linear predictor 
          lin_eff[[i]] <- attr(b$pterms, "term.labels")
          smooth_eff[[i]] <- list(NULL)
          idx_smooth_eff[[i]] <- list(NULL)
          idx_lin_eff[[i]] <- which(rhs.vars(obj[[i]]) %in% lin_eff[[i]])
          X_jk[[i]] <- b$X
          S_jkl[[i]] <-  diag(1,ncol(b$X)-1)
          n_eff[[i]] <- 1 + length(lin_eff[[i]])
          n_eff_lin[[i]] <- length(lin_eff[[i]])
          n_eff_smooth[[i]] <- 0
          vs_eff[[i]] <- list(NULL)
          mat[[i]] <- matrix(0,length(attr(b$pterms, "term.labels")), 3)
          mat[[i]][,1] <- c(2:(length(attr(b$pterms, "term.labels") )+1))
          mat[[i]][,2] <- c(2:(length(attr(b$pterms, "term.labels") )+1))
          mat[[i]][,3] <- rep(0,length(attr(b$pterms, "term.labels") )) 
          mat[[i]] <- rbind(matrix(0,1,3), mat[[i]])
        } else {                                    # Both linear and smooth effects in linear predictor
          lin_eff[[i]] <- attr(b$pterms, "term.labels")  
          idx_lin_eff[[i]] <- which(rhs.vars(obj[[i]]) %in% lin_eff[[i]])
          smooth_eff[[i]] <- unlist(lapply(b$smooth, function(x) x$label))
          vs_eff[[i]] <- unlist(lapply(b$smooth, function(x) x$vn))   
          idx2 <- rep(NA, length(vs_eff[[i]]))
          for(s in 1:length(vs_eff[[i]])){
            idx2[s] <- grepl(vs_eff[[i]][s], smooth_eff[[i]][s])
          }
          idx3 <- rep(FALSE, length(lin_eff[[i]]) + length( smooth_eff[[i]]))
          idx3[!(1:length(idx3)) %in% idx_lin_eff[[i]] ] <- idx2
          
          idx_smooth_eff[[i]] <- which(rhs.vars(obj[[i]]) %in% smooth_eff[[i]] | idx3)    
          
          mat[[i]] <- matrix(0, length(attr(b$pterms, "term.labels") ) + length( b$smooth),3)
          mat[[i]][,1] <- c(2:(length(attr(b$pterms, "term.labels") )+1), unlist(lapply(b$smooth, function(x) x$first.para)))               #indices for accessing 
          mat[[i]][,2] <- c(2:(length(attr(b$pterms, "term.labels") )+1), unlist(lapply(b$smooth, function(x) x$last.para)))
          mat[[i]][,3] <- c(rep(0,length(attr(b$pterms, "term.labels") ) ),  1:length(unlist(lapply(b$smooth, function(x) x$last.para))))   #
          X_jk[[i]] <- cbind(b$X[,"(Intercept)"], do.call("cbind",lapply(1:nrow(mat[[i]]), function(x) b$X[, mat[[i]][x,1]: mat[[i]][x,2], drop = FALSE])))
          mat_S <- matrix(0, ncol(b$X)-1, ncol(b$X)-1)
          mat_S <- as.matrix(do.call("bdiag", lapply(1:nrow(mat[[i]]), function(x) if( mat[[i]][x,1]== mat[[i]][x,2]){matrix(1)} else {b$S[[ mat[[i]][x,3]]]})))
          S_jkl[[i]] <-  mat_S
          n_eff[[i]] <- 1 + length(lin_eff[[i]]) + length(smooth_eff[[i]])
          n_eff_lin[[i]] <- length(lin_eff[[i]])
          n_eff_smooth[[i]] <- length(smooth_eff[[i]])
          mat[[i]] <- rbind(matrix(0,1,3), mat[[i]])
        }
        
      }
    }
    
  }  
  
  return(list(lin_eff = lin_eff,
              idx_lin_eff = idx_lin_eff,  
              smooth_eff = smooth_eff,
              idx_smooth_eff = idx_smooth_eff,   
              vs_eff = vs_eff, 
              X_jk = X_jk, 
              S_jkl = S_jkl, 
              n_eff = n_eff, 
              n_eff_lin = n_eff_lin,
              n_eff_smooth = n_eff_smooth,
              mat = mat) )
}

# Initialises the model parameters
intialise_function <- function(X, data, regimes = "N-N", nchains = 4, K = NULL, neff_l1_lin = 0, 
                               comm_seas = FALSE, comm_trend =  FALSE, comm_AR = FALSE, AR = FALSE ){
  
  pl <- list() # number of columns 
  ps <- list()
  ptpm <- list() 
  
  pl <- lapply(X[1:K], function(x) ncol(x[[1]]))
  ps <- lapply(X[1:K], function(x) ncol(x[[2]]))
  
  aaa <- expand.grid(c(0:(K-1)), c(0:(K-1)))
  aaa <- aaa[abs(apply(aaa, 1, diff))==1,]
  aaa <- apply(format(aaa), 1, paste, collapse = "")
  
  ptpm <- lapply((K+1):(length(aaa)+K), function(x) ncol(X[[x]]))
  
  initialise_param <- vector("list", 4)
  set.seed(23) #Seed used to obtain initial values 
  
  qy <- quantile(data$y, prob = c(1:K/(K+1)))
  ord_quant <- sort(c(0,(1:(K-1))/(K), (1:(K))/(K + 1),1))
  quant2 <- quantile(data$y, prob = ord_quant)
  quant2[1] <- quant2[1] - 0.0001
  sd <- lapply(c(1,(2*(1:(K-1))+1)), function(x) sd(data$y[data$y > quant2[x] &  data$y <= quant2[x + 2]]))
  
  for(j in 1:nchains){
    initialise_param_temp <- list()
    initialise_param_tempS <- list()
    initialise_param_tempTp <- list()
    
    if ( regimes == "N-rW" ){
      initialise_param_tempSh <- list(array((runif(1, -3, -2))))  
      names(initialise_param_tempSh) <- paste0("betaSh", K)
    }
    
    # for(k in 1:K){
    #   if(pl[[k]] > 1){
    #     if(AR){
    #       if(comm_AR & comm_trend & comm_seas){
    #         if(k == 1){
    #           initialise_param_temp[[k]] <- c(rnorm(1, qy[k], .1), rep(0, pl[[k]] - 1))
    #         } else {
    #           initialise_param_temp[[k]] <- array(rnorm(1, qy[k], .1))
    #         } 
    #       } else if(comm_seas & !comm_trend & !comm_AR){
    #         if(k == 1){
    #           initialise_param_temp[[k]] <- c(rnorm(1, qy[k], .1), rep(0, pl[[k]] - 1))
    #         } else {
    #           initialise_param_temp[[k]] <- c(rnorm(1, qy[k], .1), rep(0, neff_l1_lin))
    #         }  
    #       } else if(!comm_AR & !comm_trend & !comm_seas){
    #         initialise_param_temp[[k]] <- c(rnorm(1, qy[k], .1), rep(0, pl[[k]] - 1))
    #       }
    #     } else {
    #       if(comm_seas & comm_trend){
    #         if(k == 1){
    #            initialise_param_temp[[k]] <- c(rnorm(1, qy[k], .1), rep(0, pl[[k]] - 1))
    #         } else {
    #           initialise_param_temp[[k]] <- array(rnorm(1, qy[k], .1))
    #         }  
    #       } else if(comm_seas & !comm_trend){
    #         if(k == 1){
    #           initialise_param_temp[[k]] <- c(rnorm(1, qy[k], .1), rep(0, pl[[k]] - 1))
    #         } else {
    #           initialise_param_temp[[k]] <- c(rnorm(1, qy[k], .1), rep(0, neff_l1_lin))
    #         }  
    #       } else if(!comm_seas & !comm_trend){
    #         initialise_param_temp[[k]] <- c(rnorm(1, qy[k], .1), rep(0, pl[[k]] - 1))
    #       }
    #     } 
    #   } else {
    #     initialise_param_temp[[k]] <- array(rnorm(1, qy[k], .1)) 
    #   }
    #   if(ps[[k]] > 1){
    #     initialise_param_tempS[[k]] <-  c(rnorm(1, log(sd[[k]]), .1), rep(0, ps[[k]] - 1))
    #   } else {
    #     initialise_param_tempS[[k]] <-  array(rnorm(1, log(sd[[k]]), .1))
    #   }
    # }

    for(k in 1:K){
      if(pl[[k]] > 1){
        if(AR){
          if(comm_AR & comm_trend & comm_seas){
            if(k == 1){
              initialise_param_temp[[k]] <- c(runif(1, qy[k]-1, qy[k]+1), rep(0, pl[[k]] - 1))
            } else {
              initialise_param_temp[[k]] <- array(runif(1, qy[k]-1, qy[k]+1))
            } 
          } else if(comm_seas & !comm_trend & !comm_AR){
            if(k == 1){
              initialise_param_temp[[k]] <- c(runif(1, qy[k]-1, qy[k]+1), rep(0, pl[[k]] - 1))
            } else {
              initialise_param_temp[[k]] <- c(runif(1, qy[k]-1, qy[k]+1), rep(0, neff_l1_lin))
            }  
          } else if(!comm_AR & !comm_trend & !comm_seas){
            initialise_param_temp[[k]] <- c(runif(1, qy[k]-1, qy[k]+1), rep(0, pl[[k]] - 1))
          }
        } else {
          if(comm_seas & comm_trend){
            if(k == 1){
              initialise_param_temp[[k]] <- c(runif(1, qy[k]-1, qy[k]+1), rep(0, pl[[k]] - 1))
            } else {
              initialise_param_temp[[k]] <- array(runif(1, qy[k]-1, qy[k]+1))
            }  
          } else if(comm_seas & !comm_trend){
            if(k == 1){
              initialise_param_temp[[k]] <- c(runif(1, qy[k]-1, qy[k]+1), rep(0, pl[[k]] - 1))
            } else {
              initialise_param_temp[[k]] <- c(runif(1, qy[k]-1, qy[k]+1), rep(0, neff_l1_lin))
            }  
          } else if(!comm_seas & !comm_trend){
            initialise_param_temp[[k]] <- c(runif(1, qy[k]-1, qy[k]+1), rep(0, pl[[k]] - 1))
          }
        } 
      } else {
        initialise_param_temp[[k]] <- array(runif(1, qy[k]-1, qy[k]+1)) 
      }
      if(ps[[k]] > 1){
        initialise_param_tempS[[k]] <-  c(runif(1, 0.5, 1.5), rep(0, ps[[k]] - 1))
      } else {
        initialise_param_tempS[[k]] <-  array(runif(1, 0.5, 1.5))
      }
    }
    
    names(initialise_param_temp) <- paste0("betaL", 1:K)
    names(initialise_param_tempS) <- paste0("betaS", 1:K)
    
    for(k in 1:length(ptpm)){
      if(ptpm[[k]] > 1){
        initialise_param_tempTp[[k]] <- c(runif(1, -3, -2), rep(0, ptpm[[k]]-1))
      } else {
        initialise_param_tempTp[[k]] <- array(runif(1, -3, -2))
      }   
    }
    names(initialise_param_tempTp) <-   paste0("betaTp", aaa)
    
    if ( regimes == "N-rW" ){
      initialise_param[[j]] <- c(initialise_param_temp, initialise_param_tempS, initialise_param_tempTp, initialise_param_tempSh)
    } else {
      initialise_param[[j]] <- c(initialise_param_temp, initialise_param_tempS, initialise_param_tempTp)
    }  
  }  
  
  return(initialise_param)
} 

# Distinguing the last regime distribution
assign_type <- function(reg){
  if(reg == "Gaussian") {
    typeD <- 1
  } else if (reg == "Gumbel"){
    typeD <- 2
  } else if (reg == "Frechet"){
    typeD <- 3
  } else if (reg == "rWeibull"){
    typeD <- 4
  } else if (reg == "GEV"){
    typeD <- 5
  }
  
  return(typeD)
}

# Crates all the objects needed for fitting the model, fits the model and 
# returns the stan object

fit_stan <- function(data, data_pred, mod_foo, regimes,  pi00init = c(0.5, 0.5), K = 2, iter, init_par = "random", 
                     ncores = 4, adapt_delta = 0.95, max_treedepth = 10, stepsize = 0.1, 
                     comm_seas = FALSE, comm_trend = FALSE, comm_AR = FALSE, AR = FALSE){
  
  res <- get_foo(mod_foo, K = K)
  idx_1stMay <- which(data$sc_doy==0)

  # TP combinations
  aaa <- expand.grid(c(0:(K-1)), c(0:(K-1)))
  aaa <- aaa[abs(apply(aaa, 1, diff))==1,]
  aaa <- apply(format(aaa), 1, paste, collapse = "")
  
  # Setting hyperparameters according the previous estimate
  
  out <- process_modfoo(res$out2, data = data, K = K)

  dtype <- lapply(1:K, function(x) assign_type(res$regimes[x]))
  names(dtype) <- paste0("typeD", 1:K)
  
  pl <- lapply(out$X_jk[1:K][, drop = FALSE], function(x) ncol(x[[1]]))
  names(pl) <- paste0("pl",1:K)
  ps <- lapply(out$X_jk[1:K][, drop = FALSE], function(x) ncol(x[[2]]))
  names(ps) <- paste0("ps",1:K)
  ptpm <- lapply((K+1):(length(aaa)+K), function(x) ncol(out$X_jk[[x]]))
  names(ptpm) <- paste0("ptpm", aaa)
  
  Xl <- lapply(out$X_jk[1:K][, drop = FALSE],"[[",1) 
  names(Xl) <- paste0("Xl", 1:K)
  Xs <- lapply(out$X_jk[1:K][, drop = FALSE],"[[",2) 
  names(Xs) <- paste0("Xs", 1:K)
  Xtpm <- out$X_jk[(K+1):(length(aaa)+K)][, drop = FALSE]
  names(Xtpm) <- paste0("Xtpm", aaa)
  
  neff_l <- lapply(out$n_eff[1:K],"[[",1) # includes intercept
  names(neff_l) <- paste0("neff_l", 1:K)
  neff_s <- lapply(out$n_eff[1:K],"[[",2) # includes intercept
  names(neff_s) <- paste0("neff_s", 1:K)
  neff_ptpm <- out$n_eff[(K+1):(length(aaa)+K)]
  names(neff_ptpm) <- paste0("neff_ptpm", aaa)
  
   quant <- 1:K/(K+1)
   bl <- lapply(1:K, function(x) array(c(quantile(data$y, prob = quant[x]), rep(0,pl[[x]]-1))))
   names(bl) <- paste0("hp_mu_bl", 1:K)
   bs <- lapply(1:K, function(x) array(c(log(sd(data$y)), rep(0,ps[[x]]-1))))
   names(bs) <- paste0("hp_mu_bs", 1:K)
   btp <- lapply(1:length(aaa),function(x) array(c(-1,rep(0,ptpm[[x]]-1))))
   names(btp) <- paste0("hp_mu_btpm", aaa)
   
   Sbl <- lapply(1:K, function(x) t(as.matrix(chol(nearPD(out$S_jkl[[x]][[1]])$mat))))
   names(Sbl) <- paste0("S_bl", 1:K)
   Sbs <- lapply(1:K, function(x) t(as.matrix(chol(nearPD(out$S_jkl[[x]][[2]])$mat))))
   names(Sbs) <- paste0("S_bs", 1:K)
   Sbtpm <- lapply((K+1):(length(aaa)+K), function(x) t(as.matrix(chol(nearPD(out$S_jkl[[x]])$mat))))
   names(Sbtpm) <- paste0("S_btpm", aaa)
   
   ql <- lapply(1:K, function(x) ncol(out$S_jk[[x]][[1]]))
   names(ql) <- paste0("ql", 1:K)
   qs <- lapply(1:K, function(x) ncol(out$S_jk[[x]][[2]]))
   names(qs) <- paste0("qs", 1:K)
   qtpm <- lapply((K+1):(length(aaa)+K), function(x) ncol(out$S_jk[[x]]))
   names(qtpm) <- paste0("qtpm", aaa)

   matL <- lapply(1:K, function(x) out$mat[[x]][[1]])
   names(matL) <- paste0("matL", 1:K)
   matS <- lapply(1:K, function(x) out$mat[[x]][[2]])
   names(matS) <- paste0("matS", 1:K)
   matTpm <- lapply((K+1):(length(aaa)+K), function(x) out$mat[[x]])
   names(matTpm) <- paste0("matTpm", aaa)
   
   neff_l_lin <- lapply(out$n_eff_lin[1:K],"[[",1) 
   names(neff_l_lin) <- paste0("neff_l", 1:K, "_lin")
   neff_s_lin <- lapply(out$n_eff_lin[1:K],"[[",2) 
   names(neff_s_lin) <- paste0("neff_s", 1:K, "_lin")
   neff_ptpm_lin <- out$n_eff_lin[(K+1):(length(aaa)+K)]
   names(neff_ptpm_lin) <- paste0("neff_ptpm", aaa, "_lin")

   neff_l_smooth <- lapply(out$n_eff_smooth[1:K],"[[",1) 
   names(neff_l_smooth) <- paste0("neff_l", 1:K, "_smooth")
   neff_s_smooth <- lapply(out$n_eff_smooth[1:K],"[[",2) 
   names(neff_s_smooth) <- paste0("neff_s", 1:K, "_smooth")
   neff_ptpm_smooth <- out$n_eff_smooth[(K+1):(length(aaa)+K)]
   names(neff_ptpm_smooth) <- paste0("neff_ptpm", aaa, "_smooth")
   
   stan_dat_rW <- NULL
   
   if(dtype[[K]] == 4) {
     psh <- list(ncol(out$X_jk[[K]][[3]]))
     names(psh) <- paste0("psh",K)
     
     Xsh <- list(out$X_jk[[K]][[3]])
     names(Xsh) <- paste0("Xsh", K)
     
     neff_sh <- lapply(out$n_eff[K],"[[",3) # includes intercept
     names(neff_sh) <- paste0("neff_sh", K)
     
     bsh <- list(array(c(-2, rep(0,ncol(out$X_jk[[K]][[3]])-1))))
     names(bsh) <- paste0("hp_mu_bsh", K)
     
     Sbsh <- list(t(chol(out$S_jkl[[K]][[3]])))
     names(Sbsh) <- paste0("S_bsh", K)
     
     qsh <- list(ncol(out$S_jk[[K]][[3]]))
     names(qsh) <- paste0("qsh", K)
     
     matSh <- list(out$mat[[K]][[3]])
     names(matSh) <- paste0("matSh", K)
     
     neff_sh_lin <- list(out$n_eff_lin[[K]][[3]])
     names(neff_sh_lin) <- paste0("neff_sh", K, "_lin")
     
     neff_sh_smooth <- list(out$n_eff_smooth[[K]][[3]])
     names(neff_sh_smooth) <- paste0("neff_sh", K, "_smooth")
     
     stan_dat_rW <- c(psh, Xsh, neff_sh, bsh, Sbsh, qsh, matSh, neff_sh_lin, neff_sh_smooth)
   }
   
  stan_dat <- c(dtype, 
     pl, ps, ptpm, 
     Xl, Xs, Xtpm, 
     neff_l, neff_s, neff_ptpm, 
     bl, bs, btp,
     Sbl, Sbs, Sbtpm,
     ql, qs, qtpm, 
     matL, matS, matTpm, 
     neff_l_lin, neff_s_lin, neff_ptpm_lin, 
     neff_l_smooth, neff_s_smooth, neff_ptpm_smooth, 
     stan_dat_rW, 
     with(data, list(T = length(y), 
                     y = y, 
                     pi0 = pi00init,
                     n_may1 = length(idx_1stMay),
                     may1 = idx_1stMay))
   )
  
  if(dtype[[K]] == 1) model <- "Normal"
  if(dtype[[K]] == 4) model <- "rWeibull"
  if(dtype[[K]] == 2) model <- "Gumbel"
  
  init_par_type <- init_par
  if(init_par == "custom" | init_par == "MAP"){
    if(AR){
      if(comm_seas & comm_trend & comm_AR){
        init_par <- intialise_function(X = out$X_jk, data, neff_l1_lin =  0, K = K, regimes = regimes,
                                       AR = AR, comm_seas = comm_seas,  comm_trend = comm_trend, comm_AR = comm_AR)
      } else if(comm_seas & !comm_trend & !comm_AR){
        init_par <- intialise_function(X = out$X_jk, data, neff_l1_lin =  out$n_eff_lin[[1]][[1]], K = K, regimes = regimes,
                                       AR = AR, comm_seas = comm_seas,  comm_trend = comm_trend, comm_AR = comm_AR)
      } else if (!comm_seas & !comm_trend & !comm_AR){
        init_par <- intialise_function(X = out$X_jk, data, neff_l1_lin = pl[[1]][1], K = K, regimes = regimes,
                                       AR = AR, comm_seas = comm_seas,  comm_trend = comm_trend, comm_AR = comm_AR)
      }
    } else {
      if(comm_seas & comm_trend){
        init_par <- intialise_function(X = out$X_jk, data, neff_l1_lin =  0, K = K, regimes = regimes,
                                       AR = AR, comm_seas = comm_seas,  comm_trend = comm_trend, comm_AR = comm_AR)
      } else if(comm_seas & !comm_trend){
        init_par <- intialise_function(X = out$X_jk, data, neff_l1_lin =  out$n_eff_lin[[1]][[1]], K = K, regimes = regimes,
                                       AR = AR, comm_seas = comm_seas,  comm_trend = comm_trend, comm_AR = comm_AR)
      } else if (!comm_seas & !comm_trend){
        init_par <- intialise_function(X = out$X_jk, data, neff_l1_lin = pl[[1]][1], K = K, regimes = regimes,
                                       AR = AR, comm_seas = comm_seas,  comm_trend = comm_trend, comm_AR = comm_AR)
      }   
    }
  } else {
    init_par <- "random"  
  }
  

  build_stan_model(K, regimes, comm_seas, comm_trend, comm_AR, AR, data)
  
  
  if(init_par_type == "MAP"){
    if(AR){
      comp_model <- cmdstan_model(paste0("functions/MS_", regimes, "_K", K, "_commseas", comm_seas, "_commtrend", comm_trend, "commAR", comm_AR, ".stan"))
    } else {
      comp_model <- cmdstan_model(paste0("functions/MS_", regimes, "_K", K, "_commseas", comm_seas, "_commtrend", comm_trend, ".stan"))
    }
    out <- comp_model$optimize(data = stan_dat, 
                               init = list(init_par[[1]]), 
                               jacobian=TRUE)
    init_mle <- out$mle()
    idxInterc <- order(sapply(1:K, function(x) init_mle[grepl(paste0("betaL", x, ".*new"), names(init_mle))][1]))
    
    for(i in 1:4){
    init <- list()
    if(AR){
      if(comm_seas & comm_trend & comm_AR){
        idx <- c(ql[[1]] -  (neff_l_lin[[1]] - 1) - 1, lapply(2:K, function(x) 0))
        betaL <- c(lapply(idxInterc[1], function(x) c(init_mle[grepl(paste0("betaL",x,".*new"), 
                                                          names(init_mle))][1:(neff_l[[x]]-1)]+runif(1,-.1,.1), 
                                           rep(0, idx[[x]]))),
                   lapply(idxInterc[-1], function(x) c(init_mle[grepl(paste0("betaL",x,".*new"), 
                                                            names(init_mle))][1]+runif(1,-.1,.1), 
                                             rep(0, idx[[x]]))))
      } else if(comm_seas & !comm_trend & !comm_AR){
        idx <- c(ql[[1]] -  (neff_l_lin[[1]] - 1) - 1, lapply(2:K, function(x) (1-as.numeric(comm_seas)) * (ql[[x]] - (neff_l_lin[[x]] - 1) - 1)))
        betaL <- lapply(1:K, function(x) c(init_mle[grepl(paste0("betaL",idxInterc[x],".*new"), 
                                                          names(init_mle))][1:(neff_l[[x]]-1)]+runif(1,-.1,.1), 
                                           rep(0, idx[[x]])))
      } else if (!comm_seas & !comm_trend & !comm_AR){
        idx <- c(ql[[1]] -  (neff_l_lin[[1]] - 1) - 1, lapply(2:K, function(x) (1-as.numeric(comm_seas)) * (ql[[x]] - (neff_l_lin[[x]] - 1) - 1)))
        betaL <- lapply(idxInterc, function(x) c(init_mle[grepl(paste0("betaL",x,".*new"), 
                                                          names(init_mle))][1:(neff_l[[x]]-1)]+runif(1,-.1,.1), 
                                           rep(0, idx[[x]])))
      } 
    } else {
      if(comm_seas & comm_trend){
        idx <- c(ql[[1]] -  (neff_l_lin[[1]] - 1) - 1, lapply(2:K, function(x) 0))
        betaL <- c(lapply(idxInterc[1], function(x) c(init_mle[grepl(paste0("betaL",x,".*new"), 
                                                          names(init_mle))][1:(neff_l[[x]]-1)]+runif(1,-.1,.1), 
                                           rep(0, idx[[x]]))),
                 lapply(idxInterc[-1], function(x) c(init_mle[grepl(paste0("betaL",x,".*new"), 
                                               names(init_mle))][1]+runif(1,-.1,.1), 
                                rep(0, idx[[x]]))))  
      } else if(comm_seas & !comm_trend){
        idx <- c(ql[[1]] -  (neff_l_lin[[1]] - 1) - 1, lapply(2:K, function(x) (1-as.numeric(comm_seas)) * (ql[[x]] - (neff_l_lin[[x]] - 1) - 1)))
        betaL <- lapply(idxInterc, function(x) c(init_mle[grepl(paste0("betaL",x,".*new"), 
                                                          names(init_mle))][1:(neff_l[[x]]-1)]+runif(1,-.1,.1), 
                                           rep(0, idx[[x]])))
      } else if(!comm_seas & !comm_trend){
        idx <- c(ql[[1]] -  (neff_l_lin[[1]] - 1) - 1, lapply(2:K, function(x) (1-as.numeric(comm_seas)) * (ql[[x]] - (neff_l_lin[[x]] - 1) - 1)))
        betaL <- lapply(idxInterc, function(x) c(init_mle[grepl(paste0("betaL",x,".*new"), 
                                                          names(init_mle))][1:(neff_l[[x]]-1)]+runif(1,-.1,.1), 
                                           rep(0, idx[[x]])))
      }
    }
  
    names(betaL) <- paste0("betaL",1:K)
    betaS <- lapply(idxInterc, function(x) c(init_mle[grepl(paste0("betaS",x,".*new"), 
                                                      names(init_mle))][1:(neff_s[[x]]-1)]+runif(1,-.1,.1), 
                                       rep(0, qs[[x]]-1)))
    names(betaS) <- paste0("betaS",1:K)
   
     if(any(idxInterc!=c(1:K))){
      print(paste("MAP switch", idxInterc))
    }
     
    betaTp <- lapply(1:length(aaa), function(x) array(runif(1, -3, -2)))
    
    names(betaTp) <- paste0("betaTp",aaa)
    
    init <- c(betaL, betaS, betaTp)
    
    if (dtype[[K]] == 4 ){
     betaSh  <-  list(init_mle[grepl(paste0("betaSh", K) , names(init_mle))]+runif(1,-.1,.1))
    names(betaSh) <- paste0("betaSh", K)
    init <- c(init, betaSh) 
    }
    init[sapply(init, length)==1] <- lapply(init[sapply(init, length)==1], as.array)
    

      init_par[[i]] <- init
    }
  }
  
  pars <- c(paste0("betaL",1:K, "_new"),
            paste0("betaS",1:K, "_new"),
            if(regimes == "N-rW")   paste0("betaSh",K),
            paste0("betaTp", aaa, "_new"),
            "ltauL1", "pit", "log_lik", "lp__")

  if(AR){
    out2 <- stan(file = paste0("functions/MS_", regimes, "_K", K, "_commseas", comm_seas, "_commtrend", comm_trend, "commAR", comm_AR, ".stan"),
                 data = stan_dat, iter = iter,
                 cores = ncores, init = init_par, include = TRUE,
                 pars = pars)  
  } else {  
    out2 <- stan(file = paste0("functions/MS_", regimes, "_K", K, "_commseas", comm_seas, "_commtrend", comm_trend, ".stan"),
                 data = stan_dat, iter = iter,
                 cores = ncores, 
                 init = init_par, include = TRUE,
                 pars = pars)
  }

  return(out2)
}

# Creates the stan file to fit the model
build_stan_model <- function(K, regimes, comm_seas, comm_trend, comm_AR, AR, data){
  aaa <- expand.grid(c(0:(K-1)), c(0:(K-1)))
  aaa <- aaa[abs(apply(aaa, 1, diff))==1,]
  aaa <- apply(format(aaa), 1, paste, collapse = "")
  n_p <- 2 + 2 + (K-2) * 3
  
  if(regimes == "N-N" | regimes == "N-G"){
    function_block <- c(
      "functions{",
      "real gev_gumbel_lpdf(real y, real mu, real sigma) {",
      "    real z;",
      "    real lp;",
      "    z = -(y - mu) / sigma;",
      "    lp = z - exp(z);",
      "    return lp -  log(sigma);",
      "  }",
      " ",
      "  real gev_gumbel_rng(real mu, real sigma){",
      "    return mu - sigma * log(-log(uniform_rng(0,1)));",
      "  }",
      "}")
  } else if (regimes == "N-rW"){
    function_block <- c(
      "functions{",
      "  real gev_rWei_lpdf(real y, real mu, real sigma, real zeta) {",
      "    real z;",
      "    real lp;",
      "    z = 1 + (y - mu) * zeta / sigma;",
      "    if(z > 0) lp = -log(sigma) - (1 + 1 / zeta) * log(z) - z^( 1 / abs(zeta)); // pow(z,);",
      "    if(z < 0) lp = -1e6;",
      "    return lp;",
      "  }", 
      "",
      " real gev_rWei_rng(real mu, real sigma, real zeta){",
      "   return mu - (sigma / zeta) * (1 - (-log(uniform_rng(0,1))) ^ abs(zeta) );",
      " }",
      "}",
      ""
    )
  }
  
  type <- paste0("  int typeD", 1:K, ";")
  pl <- paste0("  int<lower=0> pl", 1:K, ";")
  ps <- paste0("  int<lower=0> ps", 1:K, ";")
  if(regimes == "N-rW") psh <- paste0("  int<lower=0> psh", K, ";")
  ptpm <- paste0("  int<lower=1> ptpm", aaa, ";")
  
  neff_l <- paste0("  int<lower=0> neff_l",1:K,";")
  neff_s <- paste0("  int<lower=0> neff_s",1:K,";")
  if(regimes == "N-rW") neff_sh <- paste0("  int<lower=0> neff_sh", K,";")
  neff_ptpm <- paste0("  int<lower=0> neff_ptpm", aaa,";")
  
  neff_l_lin <- paste0("  int<lower=0> neff_l", 1:K, "_lin;")
  neff_s_lin <- paste0("  int<lower=0> neff_s", 1:K, "_lin;")
  if(regimes == "N-rW") neff_sh_lin <- paste0("  int<lower=0> neff_sh", K,"_lin;")
  neff_ptpm_lin <- paste0("  int<lower=0> neff_ptpm", aaa, "_lin;")
  
  neff_l_smooth <- paste0("  int<lower=0> neff_l", 1:K, "_smooth;")
  neff_s_smooth <- paste0("  int<lower=0> neff_s", 1:K, "_smooth;")
  if(regimes == "N-rW") neff_sh_smooth  <- paste0("  int<lower=0> neff_sh",K,"_smooth;")
  neff_ptpm_smooth <- paste0(" int<lower=0> neff_ptpm", aaa, "_smooth;")
  
  pi0 <- paste0("  simplex[", K, "] pi0;");
  
  Xl <- paste0("  matrix[T, pl", 1:K,"] Xl", 1:K, ";")
  Xs <- paste0("  matrix[T, ps", 1:K,"] Xs", 1:K, ";")
  if(regimes == "N-rW") Xsh <- paste0("  matrix[T, psh", K,"] Xsh", K, ";")
  Xtpm <-  paste0("  matrix[T, ptpm", aaa,"] Xtpm", aaa, ";")
  
  matL <- paste0("  array[neff_l", 1:K, ", 3] int matL", 1:K, ";")
  matS <- paste0("  array[neff_s", 1:K, ", 3] int matS", 1:K, ";")
  if(regimes == "N-rW") matSh <- paste0("  array[neff_sh", K, ", 3] int matSh", K, ";")
  matTpm <- paste0("  array[neff_ptpm", aaa, ", 3] int matTpm", aaa, ";")
  
  hp_mu_bl <- paste0("  vector[pl", 1:K,"] hp_mu_bl", 1:K, ";")
  hp_mu_bs <- paste0("  vector[ps", 1:K,"] hp_mu_bs", 1:K, ";")
  if(regimes == "N-rW") hp_mu_bsh <- paste0("  vector[psh", K,"] hp_mu_bsh", K, ";")
  hp_mu_btpm <- paste0("  vector[ptpm", aaa,"] hp_mu_btpm", aaa, ";")
  
  ql <- paste0("  int<lower=0> ql", 1:K, ";")
  qs <- paste0("  int<lower=0> qs", 1:K, ";")
  if(regimes == "N-rW") qsh <- paste0("  int<lower=0> qsh", K, ";")
  qtpm <- paste0("  int<lower=0> qtpm", aaa, ";")
  
  S_bl <- paste0("  matrix[ql", 1:K, ", ql",1:K,"] S_bl", 1:K, ";")
  S_bs <- paste0("  matrix[qs", 1:K, ", qs",1:K,"] S_bs", 1:K, ";")
  if(regimes == "N-rW") S_bsh <- paste0("  matrix[qsh", K, ", qsh", K,"] S_bsh", K, ";")
  S_btpm <- paste0("  matrix[qtpm", aaa, ", qtpm", aaa,"] S_btpm", aaa, ";")
  
  
  data_block <- c(
    "data {", 
    "  int<lower=1> T;",  
    "  vector[T] y;", "",
    c(type, "", pl, "", ps, "", if(regimes == "N-rW") psh, "", ptpm, "", 
      neff_l, "", neff_s, "", if(regimes == "N-rW") neff_sh, "", neff_ptpm, "", 
      neff_l_lin, "", neff_s_lin, "", if(regimes == "N-rW") neff_sh_lin, "", neff_ptpm_lin, "",
      neff_l_smooth, "", neff_s_smooth, "", if(regimes == "N-rW") neff_sh_smooth, "", neff_ptpm_smooth, "",
      pi0, "", Xl, "", Xs, "", if(regimes == "N-rW") Xsh, "", Xtpm, "", 
      matL, "", matS, "", if(regimes == "N-rW") matSh, "", matTpm, "", 
      hp_mu_bl, "", hp_mu_bs, "", if(regimes == "N-rW") hp_mu_bsh, "", hp_mu_btpm, "", 
      ql, "", qs, "", if(regimes == "N-rW") qsh, "", qtpm, "", 
      S_bl, "", S_bs, "", if(regimes == "N-rW") S_bsh, "",S_btpm), "",
    "  int<lower=1> n_may1;", 
    "  array[n_may1] int may1; ",
    "}",
    ""
  )
  
  
  
  transformed_data_block <- c( 
    "transformed data {",
    "  //Nothing",  
    "}",
    ""
  )
  
  if(AR){ # Up to now, this flag is unuseful but if we want implent other cases reveals to be useful
    if(comm_seas & comm_trend & comm_AR){ # comm_seas = TRUE, comm_AR = TRUE, comm_trend = TRUE 
      betaL <-  c(paste0("  vector[pl", 1, "] betaL",1, ";"), 
                  paste0("  vector[1] betaL", 2:K, ";"))
    } else if (comm_seas & !comm_trend & comm_AR){
      stop("Not implemented") 
    } else if (comm_seas & comm_trend & !comm_AR){
      stop("Not implemented") 
    } else if(comm_seas & !comm_trend & !comm_AR){ # comm_seas = TRUE, comm_AR = FALSE, comm_trend = FALSE 
      betaL <-  c(paste0("  vector[pl", 1, "] betaL",1, ";"), 
                  paste0("  vector[1 + neff_l", 2:K, "_lin] betaL", 2:K, ";"))
    } else if (!comm_seas & comm_trend & !comm_AR){
      stop("Not implemented") 
    } else if (!comm_seas & !comm_trend & comm_AR){
      stop("Not implemented") 
    } else if (!comm_seas & comm_trend & comm_AR){
      stop("Not implemented") 
    } else { # comm_seas = FALSE, comm_AR = FALSE, comm_trend = FALSE 
      betaL <- paste0("  vector[pl", 1:K, "] betaL", 1:K, ";")
    }
  } else {
    if(comm_seas & comm_trend){ #comm_seas = TRUE, comm_trend = TRUE
      betaL <-  c(paste0("  vector[pl", 1, "] betaL",1, ";"), 
                  paste0("  vector[1] betaL", 2:K, ";"))
    } else if(comm_seas & !comm_trend){  #comm_seas = TRUE, comm_trend = FALSE
      betaL <-  c(paste0("  vector[pl", 1, "] betaL",1, ";"), 
                  paste0("  vector[1 + neff_l", 2:K, "_lin] betaL", 2:K, ";"))
    } else if(!comm_seas & comm_trend){  #comm_seas = FALSE, comm_trend = TRUE
      stop("Not implemented") 
    } else { #comm_seas = FALSE, comm_trend = FALSE 
      betaL <- paste0("  vector[pl", 1:K, "] betaL", 1:K, ";")  
    }
  }
    
  betaS <- paste0("  vector[ps", 1:K, "] betaS",1:K, ";")
  
  if(regimes == "N-rW") betaSh <- paste0("  vector[psh", K, "] betaSh", K, ";")
  
  betaTp <- paste0("  vector[ptpm",aaa,"] betaTp", aaa, ";")
  if(comm_seas){ 
    ltauL <- paste0("  real<lower=0> ltauL", 1, ";")
  } else {
    ltauL <- paste0("  real<lower=0> ltauL", 1:K, ";")  
  }  
  ltauS <- paste0("  real<lower=0> ltauS", 1:K, ";")
  ltauSh <- paste0("  real<lower=0> ltauSh", K, ";")
  ltauTp <- paste0("  real<lower=0> ltauTp", aaa, ";")
  
  parameters_block <- c(
    "parameters{",
    betaL, "", betaS, "", if(regimes == "N-rW") betaSh, "", betaTp, "", ltauL, "", ltauS, "", if(regimes == "N-rW") ltauSh, "", ltauTp,  
    "}",
    ""
  )
  
  eta_loc <- paste0("  real eta_loc", 1:K, ";")
  eta_scale <- paste0("  real eta_scale", 1:K, ";")
  if(regimes == "N-rW")  eta_shape <- paste0("  real eta_shape", K, ";")
  eta_tpm <- paste0("  real eta_tpm", aaa, ";")
  pit <- paste0("  matrix[T,", K, "]", " pit;")
  
  
  delta1 <- "  vector<lower = 0, upper = 1>[1] delta1;"
  if(K > 2){
    deltak <- unlist(lapply((2):(K-1), function(x) paste0("  simplex[3] delta", x, ";")))
  }
  deltaK <- paste0("  vector<lower = 0, upper = 1>[1] delta", K, ";");
  
  f_ys <- paste0("  vector[", K, "] f_ys;")
  f <- "  vector[T] f;"
  
  p <- paste0("  vector[", n_p, "] p;")
  
  if(AR){ # Up to now, this flag is unuseful but if we want implent other cases reveals to be useful
    if(comm_seas & comm_trend & comm_AR){ # comm_seas = TRUE, comm_AR = TRUE, comm_trend = TRUE 
      betaL_new <- c("  vector[pl1] betaL1_new = betaL1;",
                     paste0("  vector[1] betaL", 2:K, "_new = betaL", 2:K, ";" ))
    } else if(comm_seas & !comm_trend & !comm_AR){ # comm_seas = TRUE, comm_AR = FALSE, comm_trend = FALSE 
      betaL <-        betaL_new <- c("  vector[pl1] betaL1_new = betaL1;",
                                     paste0("  vector[1 + neff_l", 2:K, "_lin] betaL", 2:K, "_new = betaL", 2:K, ";" ))
    } else if(!comm_seas & !comm_trend & !comm_AR) { # comm_seas = FALSE, comm_AR = FALSE, comm_trend = FALSE 
      betaL_new <- paste0("  vector[pl", 1:K, "] betaL", 1:K, "_new = betaL", 1:K, ";" )
    }
  } else {  
    if(comm_seas & comm_trend){ #comm_seas = TRUE, comm_trend = TRUE
      betaL_new <- c("  vector[pl1] betaL1_new = betaL1;",
                     paste0("  vector[1] betaL", 2:K, "_new = betaL", 2:K, ";" ))
    } else  if(comm_seas & !comm_trend){ # comm_seas = TRUE & comm_trend = FALSE
      betaL_new <- c("  vector[pl1] betaL1_new = betaL1;",
                     paste0("  vector[1 + neff_l", 2:K, "_lin] betaL", 2:K, "_new = betaL", 2:K, ";" ))
    } else if(!comm_seas & !comm_trend){ # comm_seas = FALSE & comm_trend = FALSE
      betaL_new <- paste0("  vector[pl", 1:K, "] betaL", 1:K, "_new = betaL", 1:K, ";" )
    }  
  }
  
  betaS_new <- paste0("  vector[ps", 1:K, "] betaS", 1:K, "_new = betaS", 1:K, ";")
  betaTp_new <- paste0("  vector[ptpm", aaa,"] betaTp", aaa, "_new = betaTp", aaa,";") 
  
  
  bbb <- list()
  
  if(AR){# Up to now, this flag is unuseful but if we want implent other cases reveals to be useful
    if(comm_seas & comm_trend & comm_AR){ # comm_seas = TRUE, comm_AR = TRUE, comm_trend = TRUE 
      bbb[[1]] <-  c("    if(neff_l1_smooth >= 1){", 
                     "      for(k in (1 + neff_l1_lin + 1) : neff_l1){",
                     "        int idxrow1 = matL1[k,1];",
                     "        int idxrow2 = matL1[k,2];",
                     "        betaL1_new[idxrow1:idxrow2] =  betaL1[idxrow1 : idxrow2] * ltauL1;",
                     "      }",
                     "    }",
                     "",
                     paste0("    eta_loc1 = Xl1[t,] * betaL1_new;"),
                     "") 
      for(j in 2 : K) {
        bbb[[j]] <- c(
          paste0("eta_loc", j, " = Xl", j, "[t,1] * betaL", j, "_new[1] + Xl1[t,2:pl1] * betaL1_new[2:pl1];") 
        )
      } 
    } else if(comm_seas & !comm_trend & !comm_AR){ # comm_seas = TRUE, comm_AR = FALSE, comm_trend = FALSE 
      bbb[[1]] <-  c("    if(neff_l1_smooth >= 1){", 
                     "      for(k in (1 + neff_l1_lin + 1) : neff_l1){",
                     "        int idxrow1 = matL1[k,1];",
                     "        int idxrow2 = matL1[k,2];",
                     "        betaL1_new[idxrow1:idxrow2] =  betaL1[idxrow1 : idxrow2] * ltauL1;",
                     "      }",
                     "    }",
                     "",
                     paste0("    eta_loc1 = Xl1[t,] * betaL1_new;"),
                     "") 
      for(j in 2 : K) {
        bbb[[j]] <- c(
          paste0("eta_loc", j, " = Xl", j, "[t,1:(1+neff_l", j, "_lin)] * betaL", j, "_new[1:(1+neff_l", j, "_lin)] + Xl1[t,(neff_l1_lin+2):pl1] * betaL1_new[(neff_l1_lin+2):pl1];") 
        )
      }
    } else if(!comm_seas & !comm_trend & !comm_AR){
      for(j in 1:K){
        bbb[[j]] <- c(paste0("    if(neff_l", j, "_smooth >= 1){"), 
                      paste0("      for(k in (1 + neff_l", j, "_lin + 1) : neff_l", j, "){"),
                      paste0("        int idxrow1 = matL", j, "[k,1];"),
                      paste0("        int idxrow2 = matL", j, "[k,2];"),
                      paste0("        betaL", j,"_new[idxrow1:idxrow2] =  betaL", j, "[idxrow1 : idxrow2] * ltauL", j, ";"),
                      "      }",
                      "    }",
                      "",
                      paste0("    eta_loc", j,  " = Xl", j, "[t,] * betaL", j, "_new;"), 
                      "")
      }  
    }
  } else { 
    if(comm_seas & comm_trend){ #comm_seas = TRUE, comm_trend = TRUE
      bbb[[1]] <-  c("    if(neff_l1_smooth >= 1){", 
                     "      for(k in (1 + neff_l1_lin + 1) : neff_l1){",
                     "        int idxrow1 = matL1[k,1];",
                     "        int idxrow2 = matL1[k,2];",
                     "        betaL1_new[idxrow1:idxrow2] =  betaL1[idxrow1 : idxrow2] * ltauL1;",
                     "      }",
                     "    }",
                     "",
                     paste0("    eta_loc1 = Xl1[t,] * betaL1_new;"),
                     "") 
      for(j in 2 : K) {
        bbb[[j]] <- c(
          paste0("eta_loc", j, " = Xl", j, "[t,1] * betaL", j, "_new[1] + Xl1[t,2:pl1] * betaL1_new[2:pl1];") 
        )
      } 
    } else  if(comm_seas & !comm_trend){ # comm_seas = TRUE & comm_trend = FALSE
      bbb[[1]] <-  c("    if(neff_l1_smooth >= 1){", 
                     "      for(k in (1 + neff_l1_lin + 1) : neff_l1){",
                     "        int idxrow1 = matL1[k,1];",
                     "        int idxrow2 = matL1[k,2];",
                     "        betaL1_new[idxrow1:idxrow2] =  betaL1[idxrow1 : idxrow2] * ltauL1;",
                     "      }",
                     "    }",
                     "",
                     paste0("    eta_loc1 = Xl1[t,] * betaL1_new;"),
                     "") 
      for(j in 2 : K) {
        bbb[[j]] <- c(
          paste0("eta_loc", j, " = Xl", j, "[t,1:(1+neff_l", j, "_lin)] * betaL", j, "_new[1:(1+neff_l", j, "_lin)] + Xl1[t,(neff_l1_lin+2):pl1] * betaL1_new[(neff_l1_lin+2):pl1];") 
        )
      }
    } else if(!comm_seas & !comm_trend){
      for(j in 1:K){
        bbb[[j]] <- c(paste0("    if(neff_l", j, "_smooth >= 1){"), 
                      paste0("      for(k in (1 + neff_l", j, "_lin + 1) : neff_l", j, "){"),
                      paste0("        int idxrow1 = matL", j, "[k,1];"),
                      paste0("        int idxrow2 = matL", j, "[k,2];"),
                      paste0("        betaL", j,"_new[idxrow1:idxrow2] =  betaL", j, "[idxrow1 : idxrow2] * ltauL", j, ";"),
                      "      }",
                      "    }",
                      "",
                      paste0("    eta_loc", j,  " = Xl", j, "[t,] * betaL", j, "_new;"), 
                      "")
      }  
    }
  }
  
  build_eta_loc <- do.call(c, bbb)
  
  ccc <- list()
  
  for(j in 1:K){
    ccc[[j]] <- c(paste0("    if(neff_s", j, "_smooth >= 1){"), 
                  paste0("      for(k in (1 + neff_s", j, "_lin + 1) : neff_s", j, "){"),
                  paste0("        int idxrow1 = matS", j, "[k,1];"),
                  paste0("        int idxrow2 = matS", j, "[k,2];"),
                  paste0("        betaS", j,"_new[idxrow1:idxrow2] =  betaS", j, "[idxrow1 : idxrow2] * ltauS", j, ";"),
                  "      }",
                  "    }",
                  paste0("    eta_scale", j,  " = Xs", j, "[t,] * betaS", j, "_new;"), 
                  "")
  }  
  
  build_eta_scale <- do.call(c, ccc)
  
  if(regimes == "N-rW") build_eta_shape <- paste0("    eta_shape", K,  " = Xsh", K, "[t,] * betaSh", K, ";")
  
  ddd <- list()
  for(j in 1: (K - 1)){ 
    ddd[[j]] <- paste0("    f_ys[", j,"] = normal_lpdf(y[t]| eta_loc", j, ", exp(eta_scale", j,"));")
  }
  
  if(regimes == "N-N"){ 
    ddd[[K]] <- paste0("      f_ys[", K, "] = normal_lpdf(y[t]| eta_loc", K, ", exp(eta_scale", K, "));")
    
  } 
  if(regimes == "N-G"){
    ddd[[K]] <- paste0("      f_ys[", K, "] = gev_gumbel_lpdf(y[t]| eta_loc", K, ", exp(eta_scale", K, "));")
  }
  
  if (regimes == "N-rW"){
    ddd[[K]] <-  paste0("      f_ys[", K, "] = gev_rWei_lpdf(y[t]| eta_loc", K, ", exp(eta_scale", K, "), -exp(eta_shape", K, "));")
  } 
  
  build_f_ys <- do.call(c, ddd)
  
  
  fff <- list()
  
  for(j in 1:length(aaa)){
    fff[[j]] <- c(paste0("    if(neff_ptpm", aaa[j], "_smooth >= 1){"), 
                  paste0("      for(k in (1 + neff_ptpm", aaa[j], "_lin + 1) : neff_ptpm", aaa[j], "){"),
                  paste0("        int idxrow1 = matTpm", aaa[j], "[k,1];"),
                  paste0("        int idxrow2 = matTpm", aaa[j], "[k,2];"),
                  paste0("        betaTp", aaa[j],"_new[idxrow1:idxrow2] =  betaTp", aaa[j], "[idxrow1 : idxrow2] * ltauTp", aaa[j], ";"),
                  "      }",
                  "    }",
                  paste0("    eta_tpm", aaa[j],  " = Xtpm", aaa[j], "[t,] * betaTp", aaa[j], "_new;"), 
                  "")
  }  
  
  build_eta_Tp <- do.call(c, fff)
  
  delta_res <- list()
  delta_res[[1]] <- "    delta1[1] = inv_logit(eta_tpm01);"
  
  
  if(K > 2){ 
    for(k in 2:(K-1)){
      delta_res[[k]] <- rep("", 3)
      delta_res[[k]][1]  <- paste0("    delta", k, 
                                   "[1] = exp(eta_tpm", paste0(k-1, k-2), 
                                   ")/(1+exp(eta_tpm", paste0(k-1, k-2), ") + exp(eta_tpm", paste0(k-1, k), "));")
      delta_res[[k]][2]  <-  paste0("    delta", k, 
                                    "[3] = exp(eta_tpm", paste0(k-1, k), 
                                    ")/(1+exp(eta_tpm", paste0(k-1, k-2), ") + exp(eta_tpm", paste0(k-1, k), "));")
      delta_res[[k]][3] = paste0("    delta", k, "[2] = 1 - delta", k, "[1] - delta", k, "[3];")    
    }
  }
  
  delta_res[[K]] <- paste0("    delta", K, "[1] = inv_logit(eta_tpm", paste0(K-1, K-2), ");")
  
  
  build_delta <- do.call(c, delta_res)
  
  
  p_res <- list()
  p_res[[1]] <- rep("", 2)
  p_res[[1]][1] = "      p[1] = log1m(delta1[1]) + log(pi0[1]) + f_ys[1];"
  p_res[[1]][2] = "      p[2] = log(delta2[1]) + log(pi0[2]) + f_ys[1];"
  
  if(K >2){
    p_res[[2]] <- rep("",3)
    p_res[[2]][1] <- paste0("      p[3] = log(delta", 1, "[1]) + log(pi0[", 1, "]) + f_ys[", 2, "];")
    p_res[[2]][2] <-  paste0("      p[4] = log(delta", 2, "[2]) + log(pi0[", 2, "]) + f_ys[", 2, "];")
    p_res[[2]][3] <-  paste0("      p[5] = log(delta", 3, "[1]) + log(pi0[", 3, "]) + f_ys[", 2, "];")
  }
  if(K > 3){
    count <- 0
    for(k in 3:(K-1)){
        p_res[[k]] <- rep("",3)
        p_res[[k]][1] <- paste0("      p[", k*2 + count, "] = log(delta", k-1, "[3]) + log(pi0[", k-1, "]) + f_ys[", k, "];")
        p_res[[k]][2] <-  paste0("      p[", k*2 + 1 + count, "] = log(delta", k, "[2]) + log(pi0[", k, "]) + f_ys[", k, "];")
        p_res[[k]][3] <-  paste0("      p[", k*2 + 2 + count, "] = log(delta", k+1, "[1]) + log(pi0[", k+1, "]) + f_ys[", k, "];")
      count <- count + 1
        }
    }  
  
  if(K > 2){
    p_res[[K]] <- rep("", 2)
    p_res[[K]][1] <- paste0("      p[", n_p - 1, "] = log(delta", K-1, "[3]) + log(pi0[", K - 1, "]) + f_ys[", K, "];")
    p_res[[K]][2] = paste0("      p[", n_p, "] = log1m(delta", K, "[1]) + log(pi0[", K, "]) + f_ys[", K, "];")
  } else {
    p_res[[K]] <- rep("", 2)
    p_res[[K]][1] <- paste0("      p[", n_p - 1, "] = log(delta", K-1, "[1]) + log(pi0[", K - 1, "]) + f_ys[", K, "];")
    p_res[[K]][2] = paste0("      p[", n_p, "] = log1m(delta", K, "[1]) + log(pi0[", K, "]) + f_ys[", K, "];") 
  }
  p_t1 <- do.call(c, p_res)
  
  ft_res_t1 <- "      f[t] = log_sum_exp(p);"
  
  pit_res_t1 <- list()
  pit_res_t1[[1]] <- "      pit[1,1]  = log_sum_exp(p[1:2])-f[t];"
  


  if(K>2){
    for(j in 2:(K-1)){ 
      pit_res_t1[[j]] <- paste0("      pit[1,", j,"] = log_sum_exp(p[", paste0(((j-1)*3), ":", ((j-1)*3+2)),"])-f[t];")
    }
  }
  
  pit_res_t1[[K]] <- paste0("      pit[1,", K,"] = log_sum_exp(p[", paste0((n_p-1),":", n_p), "])-f[t];")
  
  pit_t1 <- do.call(c, pit_res_t1)
  
  
  
  
  
  pt_res <- list()
  pt_res[[1]] <- rep("", 2)
  pt_res[[1]][1] = "      p[1] = log1m(delta1[1]) + pit[t-1, 1] + f_ys[1];"
  pt_res[[1]][2] = "      p[2] = log(delta2[1]) + pit[t-1, 2] + f_ys[1];"
  
  if(K >2){
    pt_res[[2]] <- rep("",3)
    pt_res[[2]][1] <- paste0("      p[3] = log(delta", 1, "[1]) + pit[t - 1,", 1, "] + f_ys[", 2, "];")
    pt_res[[2]][2] <-  paste0("      p[4] = log(delta", 2, "[2]) + pit[t - 1,", 2, "] + f_ys[", 2, "];")
    pt_res[[2]][3] <-  paste0("      p[5] = log(delta", 3, "[1]) + pit[t - 1,", 3, "] + f_ys[", 2, "];")
  }
  if(K>3){
    count <- 0
    for(k in 3:(K-1)){
        pt_res[[k]] <- rep("",3)
        pt_res[[k]][1] <- paste0("      p[", k*2 + count, "] = log(delta", k-1, "[3]) + pit[t - 1,", k-1, "] + f_ys[", k, "];")
        pt_res[[k]][2] <-  paste0("      p[", k*2 + 1 + count, "] = log(delta", k, "[2]) + pit[t - 1,", k, "] + f_ys[", k, "];")
        pt_res[[k]][3] <-  paste0("      p[", k*2 + 2 + count, "] = log(delta", k+1, "[1]) + pit[t - 1,", k+1, "] + f_ys[", k, "];")
      count <- count + 1
        } 
    }  
  
  if(K>2){
    pt_res[[K]] <- rep("", 2)
    pt_res[[K]][1] <- paste0("      p[", n_p - 1, "] = log(delta", K-1, "[3]) + pit[t - 1,", K - 1, "] + f_ys[", K, "];")
    pt_res[[K]][2] = paste0("      p[", n_p, "] = log1m(delta", K, "[1]) + pit[t - 1,", K, "] + f_ys[", K, "];")
  } else {
    pt_res[[K]] <- rep("", 2)
    pt_res[[K]][1] <- paste0("      p[", n_p - 1, "] = log(delta", K-1, "[1]) + pit[t - 1,", K - 1, "] + f_ys[", K, "];")
    pt_res[[K]][2] = paste0("      p[", n_p, "] = log1m(delta", K, "[1]) + pit[t - 1,", K, "] + f_ys[", K, "];")
  }
  
  p_t <- do.call(c, pt_res)
  
  ft_res_t <- "      f[t] = log_sum_exp(p);"
  
  pit_res_t <- list()
  pit_res_t[[1]] <- "      pit[t,1]  = log_sum_exp(p[1:2])-f[t];"
  
  if(K>2){
    for(j in 2:(K-1)){ 
      pit_res_t[[j]] <- paste0("      pit[t,", j,"] = log_sum_exp(p[", paste0(((j-1)*3), ":", ((j-1)*3+2)),"])-f[t];")
    }
  }
  pit_res_t[[K]] <- paste0("      pit[t,", K,"] = log_sum_exp(p[", paste0((n_p-1),":", n_p), "])-f[t];")
  
  pit_t <- do.call(c, pit_res_t)
  
  
  transformed_parameters_block <- c(
    "transformed parameters{", 
    eta_loc, "", eta_scale, "", if(regimes == "N-rW")  eta_shape, eta_tpm, "", 
    pit, "", delta1, "", if(K>2) deltak, "", deltaK, "", f_ys, "", f, "", p, "",
    betaL_new, "", betaS_new, "",  betaTp_new, "", 
    "  for(t in 1:T) {",
    build_eta_loc, "", build_eta_scale, "", if(regimes == "N-rW") build_eta_shape, "",
    build_f_ys, "", build_eta_Tp, "", build_delta, "",
    "    if (t == 1) {", 
    p_t1, "", ft_res_t1, "", pit_t1,  
    "    } else {",
    p_t, "", ft_res_t, "", pit_t, 
    "    }",
    "  }",
    "}",
    ""
  )
  
  eee <- list()
  
  ord_quant <- sort(c(0,(1:(K-1))/(K), (1:(K))/(K + 1),1))
  quant2 <- quantile(data$y, prob = ord_quant)
  quant2[1] <- quant2[1] - 0.0001
  sd <- lapply(c(1,(2*(1:(K-1))+1)), function(x) sd(data$y[data$y > quant2[x] &  data$y <= quant2[x + 2]]))
  
  
  if(AR){
    if(comm_seas & comm_trend & comm_AR ){
      eee[[1]] <-  c(paste0("  if(pl", 1, "> 1){"),
                     paste0("    betaL", 1, "[1] ~ normal(hp_mu_bl", 1, "[1], sd(y));"),
                     paste0("    for(j in 2:neff_l", 1, "){"),
                     paste0("      ltauL", 1, " ~  student_t(7,0,.5);"),
                     paste0("      if(j <= (1  + neff_l", 1, "_lin)){"),
                     paste0("        betaL", 1, "[j] ~ normal(hp_mu_bl", 1, "[j], 2.5);"),
                     "      } else {",
                     paste0("        int idxrow1 = matL", 1, "[j, 1];"),
                     paste0("        int idxrow2 = matL", 1, "[j, 2];"),
                     paste0("        betaL", 1, "[idxrow1 : idxrow2] ~ multi_normal_prec(hp_mu_bl", 1, "[idxrow1 : idxrow2],  S_bl", 1, "[(idxrow1-1):(idxrow2-1),(idxrow1-1):(idxrow2-1)]  * S_bl", 1, "[(idxrow1-1):(idxrow2-1),(idxrow1-1):(idxrow2-1)]');"),
                     "      }",
                     "    }",
                     "  } else {",
                     paste0("    betaL", 1, " ~ normal(hp_mu_bl", 1, "[1], sd(y));"),
                     paste0("    ltauL", 1, " ~  student_t(7,0,.5);"),  
                     "  }", "")
      for(j in 2:K){
        eee[[j]] <- c(paste0("    betaL", j, "[1] ~ normal(hp_mu_bl", j, "[1], sd(y));"))
      }  
    } else if(comm_seas & !comm_trend & !comm_AR){
      eee[[1]] <-  c(paste0("  if(pl", 1, "> 1){"),
                     paste0("    betaL", 1, "[1] ~ normal(hp_mu_bl", 1,"[1], sd(y));"), 
                     paste0("    for(j in 2:neff_l", 1, "){"),
                     paste0("      ltauL", 1, " ~  student_t(7,0,.5);"),
                     paste0("      if(j <= (1  + neff_l", 1, "_lin)){"),
                     paste0("        betaL", 1, "[j] ~ normal(hp_mu_bl", 1, "[j], 2.5);"), 
                     "      } else {",
                     paste0("        int idxrow1 = matL", 1, "[j, 1];"),
                     paste0("        int idxrow2 = matL", 1, "[j, 2];"),
                     paste0("        betaL", 1, "[idxrow1 : idxrow2] ~ multi_normal_prec(hp_mu_bl", 1, "[idxrow1 : idxrow2],  S_bl", 1, "[(idxrow1-1):(idxrow2-1),(idxrow1-1):(idxrow2-1)]  * S_bl", 1, "[(idxrow1-1):(idxrow2-1),(idxrow1-1):(idxrow2-1)]');"),
                     "      }",
                     "    }",
                     "  } else {",
                     paste0("    betaL", 1, " ~ normal(hp_mu_bl", 1, "[1], sd(y));"),
                     paste0("    ltauL", 1, " ~  student_t(7,0,.5);"),  
                     "  }", "")
      for(j in 2:K){
        eee[[j]] <- c(paste0("  if(pl", j, "> 1){"),
                      paste0("    betaL", j, "[1] ~ normal(hp_mu_bl", j, "[1],.25);"), 
                      paste0("    for(j in 2:neff_l", j, "){"),
                      paste0("      if(j <= (1  + neff_l", j, "_lin)){"),
                      paste0("        betaL", j, "[j] ~ normal(hp_mu_bl", j, "[j], 2.5);"), 
                      "      } ",
                      "    }",
                      "  } else {",
                      paste0("    betaL", j, " ~ normal(hp_mu_bl", j, "[1], sd(y));"),
                      "  }", "")
      }  
    } else if(!comm_seas & !comm_trend & !comm_AR){
      for(k in 1:K){ 
        eee[[k]]  <- c(paste0("  if(pl", k, "> 1){"),
                       paste0("    betaL", k, "[1] ~ normal(hp_mu_bl", k,"[1],sd(y));"),
                       paste0("    for(j in 2:neff_l", k, "){"),
                       paste0("      ltauL", k, " ~  student_t(7,0,.5);"),
                       paste0("      if(j <= (1  + neff_l", k, "_lin)){"),
                       paste0("        betaL", k, "[j] ~ normal(hp_mu_bl", k, "[j], 2.5);"),
                       "      } else {",
                       paste0("        int idxrow1 = matL", k, "[j, 1];"),
                       paste0("        int idxrow2 = matL", k, "[j, 2];"),
                       paste0("        betaL", k, "[idxrow1 : idxrow2] ~ multi_normal_prec(hp_mu_bl", k, "[idxrow1 : idxrow2],  S_bl", k, "[(idxrow1-1):(idxrow2-1),(idxrow1-1):(idxrow2-1)]  * S_bl", k, "[(idxrow1-1):(idxrow2-1),(idxrow1-1):(idxrow2-1)]');"),
                       "      }",
                       "    }",
                       "  } else {",
                       paste0("    betaL", k, " ~ normal(hp_mu_bl", k, "[1], sd(y));"),
                       paste0("    ltauL", k, " ~  student_t(7,0,.5);"),  
                       "  }", ""
        )
      }
    }
  } else {  
    if(comm_seas & comm_trend){
      eee[[1]] <-  c(paste0("  if(pl", 1, "> 1){"),
                     paste0("    betaL", 1, "[1] ~ normal(hp_mu_bl", 1, "[1], sd(y));"),
                     paste0("    for(j in 2:neff_l", 1, "){"),
                     paste0("      ltauL", 1, " ~  student_t(7,0,.5);"),
                     paste0("      if(j <= (1  + neff_l", 1, "_lin)){"),
                     paste0("        betaL", 1, "[j] ~ normal(hp_mu_bl", 1, "[j], 2.5);"),
                     "      } else {",
                     paste0("        int idxrow1 = matL", 1, "[j, 1];"),
                     paste0("        int idxrow2 = matL", 1, "[j, 2];"),
                     paste0("        betaL", 1, "[idxrow1 : idxrow2] ~ multi_normal_prec(hp_mu_bl", 1, "[idxrow1 : idxrow2],  S_bl", 1, "[(idxrow1-1):(idxrow2-1),(idxrow1-1):(idxrow2-1)]  * S_bl", 1, "[(idxrow1-1):(idxrow2-1),(idxrow1-1):(idxrow2-1)]');"),
                     "      }",
                     "    }",
                     "  } else {",
                     paste0("    betaL", 1, " ~ normal(hp_mu_bl", 1, "[1], sd(y));"),
                     paste0("    ltauL", 1, " ~  student_t(7,0,.5);"),  
                     "  }", "")
      for(j in 2:K){
        eee[[j]] <- c(paste0("    betaL", j, "[1] ~ normal(hp_mu_bl", j, "[1], sd(y));"))
      }  
    } else if(comm_seas & !comm_trend){
      eee[[1]] <-  c(paste0("  if(pl", 1, "> 1){"),
                     paste0("    betaL", 1, "[1] ~ normal(hp_mu_bl", 1, "[1], sd(y));"),  
                     paste0("    for(j in 2:neff_l", 1, "){"),
                     paste0("      ltauL", 1, " ~  student_t(7,0,.5);"),
                     paste0("      if(j <= (1  + neff_l", 1, "_lin)){"),
                     paste0("        betaL", 1, "[j] ~ normal(hp_mu_bl", 1, "[j], 2.5);"),
                     "      } else {",
                     paste0("        int idxrow1 = matL", 1, "[j, 1];"),
                     paste0("        int idxrow2 = matL", 1, "[j, 2];"),
                     paste0("        betaL", 1, "[idxrow1 : idxrow2] ~ multi_normal_prec(hp_mu_bl", 1, "[idxrow1 : idxrow2],  S_bl", 1, "[(idxrow1-1):(idxrow2-1),(idxrow1-1):(idxrow2-1)]  * S_bl", 1, "[(idxrow1-1):(idxrow2-1),(idxrow1-1):(idxrow2-1)]');"),
                     "      }",
                     "    }",
                     "  } else {",
                     paste0("    betaL", 1, " ~ normal(hp_mu_bl", 1, "[1], sd(y));"),
                     paste0("    ltauL", 1, " ~  student_t(7,0,.5);"),  
                     "  }", "")
      for(j in 2:K){
        eee[[j]] <- c(paste0("  if(pl", j, "> 1){"),
                      paste0("    betaL", j, "[1] ~ normal(hp_mu_bl", j, "[1], sd(y));"),  
                      paste0("    for(j in 2:neff_l", j, "){"),
                      paste0("      if(j <= (1  + neff_l", j, "_lin)){"),
                      paste0("        betaL", j, "[j] ~ normal(hp_mu_bl", j, "[j], 2.5);"),
                      "      } ",
                      "    }",
                      "  } else {",
                      paste0("    betaL", j, " ~ normal(hp_mu_bl", j, "[1], sd(y));"),
                      "  }", "")
      }  
    } else if(!comm_seas & !comm_trend){
      for(k in 1:K){ 
        eee[[k]]  <- c(paste0("  if(pl", k, "> 1){"),
                       paste0("    betaL", k, "[1] ~ normal(hp_mu_bl", k, "[1], sd(y));"),
                       paste0("    for(j in 2:neff_l", k, "){"),
                       paste0("      ltauL", k, " ~  student_t(7,0,.5);"),
                       paste0("      if(j <= (1  + neff_l", k, "_lin)){"),
                       paste0("        betaL", k, "[j] ~ normal(hp_mu_bl", k, "[j], 2.5);"),
                       "      } else {",
                       paste0("        int idxrow1 = matL", k, "[j, 1];"),
                       paste0("        int idxrow2 = matL", k, "[j, 2];"),
                       paste0("        betaL", k, "[idxrow1 : idxrow2] ~ multi_normal_prec(hp_mu_bl", k, "[idxrow1 : idxrow2],  S_bl", k, "[(idxrow1-1):(idxrow2-1),(idxrow1-1):(idxrow2-1)]  * S_bl", k, "[(idxrow1-1):(idxrow2-1),(idxrow1-1):(idxrow2-1)]');"),
                       "      }",
                       "    }",
                       "  } else {",
                       paste0("    betaL", k, " ~ normal(hp_mu_bl", k, "[1], sd(y));"),
                       paste0("    ltauL", k, " ~  student_t(7,0,.5);"),  
                       "  }", ""
        )
      }
    }
  }
  prior_betaL <- do.call(c, eee)
  
  
  kkk <- list()
  
  for(k in 1:K){ 
    kkk[[k]]  <- c(paste0("  if(ps", k, "> 1){"),
                   paste0("    betaS", k, "[1] ~ normal(hp_mu_bs", k, "[1], 0.5);"),
                   paste0("    for(j in 2:neff_s", k, "){"),
                   paste0("      ltauS", k, " ~  student_t(7,0,.5);"),
                   paste0("      if(j <= (1  + neff_s", k, "_lin)){"),
                   paste0("        betaS", k, "[j] ~ normal(hp_mu_bs", k, "[j], 1);"),
                   "      } else {",
                   paste0("        int idxrow1 = matS", k, "[j, 1];"),
                   paste0("        int idxrow2 = matS", k, "[j, 2];"),
                   paste0("        betaS", k, "[idxrow1 : idxrow2] ~ multi_normal_prec(hp_mu_bs", k, "[idxrow1 : idxrow2],  S_bs", k, "[(idxrow1-1):(idxrow2-1),(idxrow1-1):(idxrow2-1)]  * S_bs", k, "[(idxrow1-1):(idxrow2-1),(idxrow1-1):(idxrow2-1)]');"),
                   "      }",
                   "    }",
                   "  } else {",
                   paste0("    betaS", k, " ~ normal(hp_mu_bs",k,"[1],1);"), 
                   paste0("    ltauS", k, " ~  student_t(7,0,.5);"),  
                   "  }", ""
    )
  }
  
  prior_betaS <- do.call(c, kkk)
  
  if(regimes == "N-rW"){
    prior_betaSh <-c(paste0("   if(psh", K,  "> 1){  "),
                     paste0("     betaSh", K, "[1] ~ normal(hp_mu_bsh", K, "[1], 0.5);"),
                     paste0("     ltauSh", K, " ~  student_t(7,0,0.5);"),   
                     paste0("     betaSh", K, "[2 : psh", K, "] ~ multi_normal_prec(hp_mu_bsh", K, "[2 : psh", K, "],  S_bsh", K, ");"),
                     "    } else {",  
                     paste0("     betaSh", K, " ~ normal(hp_mu_bsh", K, "[1], 0.5);"), 
                     paste0("     ltauSh", K, " ~  student_t(7,0,0.5);"),
                     "  }",
                     ""
    )
  }
  
  
  jjj <- list()
  
  for(k in 1:length(aaa)){ 
    jjj[[k]]  <- c(paste0("  if(ptpm", aaa[k], "> 1){"),
                   paste0("    betaTp", aaa[k], "[1] ~ normal(hp_mu_btpm", aaa[k], "[1], 0.5);"),
                   paste0("    for(j in 2:neff_ptpm", aaa[k], "){"),
                   paste0("      ltauTp", aaa[k], " ~  student_t(7,0,.5);"),
                   paste0("      if(j <= (1  + neff_ptpm", aaa[k], "_lin)){"),
                   paste0("        betaTp", aaa[k], "[j] ~ normal(hp_mu_btpm", aaa[k], "[j], 1);"),
                   "      } else {",
                   paste0("        int idxrow1 = matTpm", aaa[k], "[j, 1];"),
                   paste0("        int idxrow2 = matTpm", aaa[k], "[j, 2];"),
                   paste0("        betaTp", aaa[k], "[idxrow1 : idxrow2] ~ multi_normal_prec(hp_mu_btpm", aaa[k], "[idxrow1 : idxrow2],  S_btpm", aaa[k], "[(idxrow1-1):(idxrow2-1),(idxrow1-1):(idxrow2-1)]  * S_btpm", aaa[k], "[(idxrow1-1):(idxrow2-1),(idxrow1-1):(idxrow2-1)]');"),
                   "      }",
                   "    }",
                   "  } else {",
                   paste0("    betaTp", aaa[k], " ~ normal(hp_mu_btpm", aaa[k], "[1], 0.5);"),
                   paste0("    ltauTp", aaa[k], " ~  student_t(7,0,.5);"),  
                   "  }", ""
    )
  }
  
  
  prior_betaTp <- do.call(c, jjj)
  
  model_block <- c(
    "model {",
    prior_betaL, prior_betaS, if(regimes == "N-rW") prior_betaSh, prior_betaTp, 
    " target += f;",
    "}",
    ""
  )
  
  generated_quantities_block <- c(
    "generated quantities {",
    "  vector[T] log_lik;",
    "  for (t in 1:T) {",
    "    log_lik[t] = f[t];",
    "  }",
    "}"
  )  
  
  if(AR){
    fileConn<- file(paste0("functions/MS_", regimes, "_K", K, "_commseas", comm_seas, "_commtrend", comm_trend, "commAR", comm_AR, ".stan"))  
  } else {
    fileConn<- file(paste0("functions/MS_", regimes, "_K", K, "_commseas", comm_seas, "_commtrend", comm_trend, ".stan"))
  }
  
  writeLines(c(function_block, data_block, transformed_data_block, parameters_block, 
               transformed_parameters_block, model_block, generated_quantities_block), fileConn)
  close(fileConn)
}


