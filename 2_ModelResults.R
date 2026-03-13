# This code extracts the results analyzed in the paper for the fitted models:
# - loo 
# - predicted temperatures (in and out of sample)
# - predicted states (in and out of sample)
# - predicted probabilities (in sample)

source("functions/ExtractPosterior.R")
source("functions/load_Packages.R")

minK <- 2
maxK <- 5
year_start <- 1995
year_end <- 2018
year_end2 <- 2024

seq_month <- rep(0:153)
label_month <- seq_month[c(1, 32, 62, 93, 124)] 
label_LASTreg <- c("N", "G", "rW")
label_LASTreg2 <- rep(label_LASTreg, length(minK:maxK))
label_reg <- c("N", "NN", "NNN", "NNNN")
label_K <- rep(minK:maxK, each = length(label_LASTreg))

#set.seed(23)
idx <- sample(1:4000, 200, replace = F) 

idx_model <- c(1,4,7,10) # To get the results for full model set 1:12

year_seq <- 1995:2024
yearSel <- 2018 
idxYearSel <- ((yearSel - 1995) * 153 + 1) : ((yearSel - 1995 + 1) * 153)  
idx2024 <- ((2024 - 2019) * 153 + 1) : ((2024 - 2019 + 1) * 153)   

## Capriva----
location <- "Capriva"
stations <- "_C"

comm_seas <- TRUE
comm_trend <- FALSE
AR <- FALSE
comm_AR <- FALSE

load("data/dataFVG/datiCapriva.Rdata")

model_names <- paste0("fit_", as.vector(t(outer(label_reg, label_LASTreg, FUN = paste0))), stations, "_s")
model_names2 <- paste0("fit_", as.vector(t(outer(label_reg, label_LASTreg, FUN = paste0))), stations)
file_select <- paste0(model_names,  ".RData") 
fileName <- list.files(path = "ResultsGST", pattern = stations)
get_model <- which(file_select%in%fileName)

# Train and Validation 
dataTrain <- dati[dati$Year <= year_end &  dati$Month >= 5 & dati$Month <= 9,]
dataVal <- dati[dati$Year > year_end & dati$Year <= year_end2 & dati$Month >= 5 & dati$Month <= 9,]

log_lik <- loo <- list()
TempPred <- list()
StatesIn <- list()
TempTraj <- list()
StatesOut <- list()
HatProbIn <- list()
Exc <- list()
Quantiles <- Quantiles5 <- list()

for(i in get_model){
  # Load STAN object 
  load(paste0("ResultsGST/",file_select[i]))
  fit <- mget(model_names2[i], envir = globalenv())[[1]] 
  rm(list=grep(model_names2[i],ls(),value=TRUE))
  
  # Get log-lik and loo
  log_lik[[i]] <- extract_log_lik(fit)
  loo[[i]] <- loo(log_lik[[i]])

  # Get Posterior
  Posterior <- ExtractPost(stan.obj = fit, K = label_K[i], last.regime = label_LASTreg2[i], 
                           comm_seas = comm_seas, comm_trend = comm_trend)
  
  # Fig 2 results 
  HatProbIn[[i]] <- as.data.frame(apply(exp(Posterior$pi),c(2,3), mean)[idxYearSel,])
  colnames(HatProbIn[[i]]) <- 1:label_K[i]
  HatProbIn[[i]]$doy <- 1:153

  Exc[[i]] <- get_excedance(K = label_K[i], Posterior$pi, year_start, year_end, last.regime = label_LASTreg2[i])
   
  # Fig 3 results + MCMC samples
  if(i == 10){ # Results of Model NNNNN
    QuantilesNEW <- extract_Quantiles(Posterior, data = dataVal, AR = AR, K = label_K[i],
                                      last.regime = "N", quantiles = c(0.025,0.25,0.5,0.75,0.975))

    QuantilesNEW <- sapply(1:label_K[i], function(x) QuantilesNEW[[x]][,label_K[i]])
    Quantiles5 <- QuantilesNEW[idx2024,]
    colnames(Quantiles5) <- paste0("qNEW5",c("0.025","0.25","0.5","0.75", "0.975"))

    predictedValues <- extract_PredSummary(Posterior, dataVal, AR = AR,
                                           K = label_K[i], last.regime = label_LASTreg2[i])


    Quantiles1 <- predictedValues$mean_reg_quantiles[,idx2024,]
    Quantiles <- t(apply(Quantiles1, 2, c))
    colnames(Quantiles) <- paste0("q",rep(1:label_K[i],each=5), c("0.025","0.25","0.5","0.75", "0.975"))

    extracted_chains <- rstan::extract(fit,
                                       pars = c("betaL1_new", "betaL2_new", "betaL3_new", "betaL4_new", "betaL5_new",
                                                "betaS1_new", "betaS2_new", "betaS3_new", "betaS4_new", "betaS5_new",
                                                "betaTp01_new", "betaTp10_new",
                                                "betaTp12_new", "betaTp21_new",
                                                "betaTp23_new", "betaTp32_new",
                                                "betaTp34_new", "betaTp43_new",
                                                "ltauL1"),
                                       permuted = FALSE)
  }

  # Temperature trajectories and generated regimes (In Sample)
  getTempTraj_insample <- TempIn(Posterior, dataTrain, AR = AR, comm_AR = comm_AR,
                                 comm_seas = comm_seas, comm_trend = comm_trend,
                                 K = label_K[i], last.regime = label_LASTreg2[i], Nsim = 200, proj = FALSE)

  TempPred[[i]] <- lapply(1: (year_end - year_start + 1), function (x){
     getTempTraj_insample[[x]]$gen$TempTrajectories$y_pred
  } )

  StatesIn[[i]] <- lapply(1: (year_end - year_start + 1), function (x){
     getTempTraj_insample[[x]]$gen$TempTrajectories$hat_regime
  } )
  
  # Temperature trajectories and generated regimes (Out of Sample)
  getTempTraj_OOS <- TempOut(Posterior, dataTrain, dataTest = dataVal, AR = AR, comm_AR = comm_AR,
                            comm_seas = comm_seas, comm_trend = comm_trend,
                            K = label_K[i], last.regime = label_LASTreg2[i], Nsim = 200, proj = TRUE)
  TempTraj[[i]] <- lapply(1: (year_end2 - year_end), function (x){
    getTempTraj_OOS[[x]]$gen$TempTrajectories$y_pred
  } )
  
  StatesOut[[i]]<-lapply(1: (year_end2 - year_end), function (x){
     getTempTraj_OOS[[x]]$gen$TempTrajectories$hat_regime
  } )

  print(i)
}

names(loo) <- names(log_lik) <- model_names
names(TempTraj) <-  model_names
names(TempPred) <-  model_names
names(StatesIn) <-  model_names
names(StatesOut) <-  model_names
names(HatProbIn) <- model_names
names(Exc) <- model_names

# Table 1 - T2
save(loo, file = "ExtractedResults/looCapriva_up2018_Bayes.Rda")
 
# Extract (200) temperature trajectories and generated states (Out of sample)
# Figure 4-5 and S8-S11
outTraj <- lapply(1:length(TempTraj), function(x){
  lapply(1:length(TempTraj[[1]]), function(z) TempTraj[[x]][[z]][, idx])
})
save(outTraj, file = "ExtractedResults/TempTrajCapriva_from2019_to2024_Bayes_all.Rda")

# Extract temperature trajectories (last regime N, Out of Sample)
# Figure 3 - S4-S5-S6-S7
TempTraj <- TempTraj[idx_model]
save(TempTraj, file = "ExtractedResults/TempTrajCapriva_from2019_to2024_Bayes.Rda")

save(StatesOut, file = "ExtractedResults/StatesOutCapriva_from2019_to2024_Bayes.Rda")

# Extract (200) temperature trajectories and generated states (In sample)
# Figure 4-5 and S8-S11
outTraj <- lapply(1:length(TempPred), function(x){
  lapply(1:length(TempPred[[1]]), function(z) TempPred[[x]][[z]][, idx])
})

outStates <- lapply(1:length(TempPred), function(x){
  lapply(1:length(TempPred[[1]]), function(z) StatesIn[[x]][[z]][, idx])
})
 
save(outTraj, file = "ExtractedResults/TempPredCapriva_up2018_Bayes.Rda")
save(outStates, file = "ExtractedResults/StatesInCapriva_up2018_Bayes.Rda")

# Fig 2
HatProbIn <- HatProbIn[c(1,4,7,10)]
save(HatProbIn, file = "ExtractedResults/HatProbInCapriva_up2018_Bayes.Rda")
save(Exc, file = "ExtractedResults/ExcCapriva_up2018_Bayes.Rda")

# Figure 3-S7
save(Quantiles5, file = "ExtractedResults/Quantiles5Capriva_from2019_to2024_Bayes.Rda")
save(Quantiles, file = "ExtractedResults/QuantilesCapriva_from2019_to2024_Bayes.Rda" )

# Table T1 - Figure S2-S3 
saveRDS(extracted_chains, file = "ExtractedResults/extracted_samples.rds")

## Trieste----

location <- "Trieste"
stations <- "_T"

comm_seas <- TRUE
comm_trend <- FALSE
AR <- FALSE
comm_AR <- FALSE

load("data/dataFVG/datiTrieste.Rdata")

model_names <- paste0("fit_", as.vector(t(outer(label_reg, label_LASTreg, FUN = paste0))), stations, "_s")
model_names2 <- paste0("fit_", as.vector(t(outer(label_reg, label_LASTreg, FUN = paste0))), stations)
file_select <- paste0(model_names,  ".RData") 
fileName <- list.files(path = "ResultsGST", pattern = stations)
get_model <- which(file_select%in%fileName)

# Train and  Validation 
dataTrain <- dati[dati$Year <= year_end &  dati$Month >= 5 & dati$Month <= 9,]
dataVal <- dati[dati$Year > year_end & dati$Year <= year_end2 & dati$Month >= 5 & dati$Month <= 9,]

log_lik <- loo <- list()
TempPred <- list()
StatesIn <- list()
TempTraj <- list()
StatesOut <- list()
HatProbIn <- list()
Exc <- list()
Quantiles <- Quantiles5 <- list()

for(i in get_model){
  # Load STAN object 
  load(paste0("ResultsGST/",file_select[i]))
  fit <- mget(model_names2[i], envir = globalenv())[[1]] 
  rm(list=grep(model_names2[i],ls(),value=TRUE))
  
  # Get log-lik and loo
  log_lik[[i]] <- extract_log_lik(fit)
  loo[[i]] <- loo(log_lik[[i]])
  
  # Get Posterior
  Posterior <- ExtractPost(stan.obj = fit, K = label_K[i], last.regime = label_LASTreg2[i], 
                           comm_seas = comm_seas, comm_trend = comm_trend)
  
  # Fig 2 results 
  HatProbIn[[i]] <- as.data.frame(apply(exp(Posterior$pi),c(2,3), mean)[idxYearSel,])
  colnames(HatProbIn[[i]]) <- 1:label_K[i]
  HatProbIn[[i]]$doy <- 1:153
  
  Exc[[i]] <- get_excedance(K = label_K[i], Posterior$pi, year_start, year_end, last.regime = label_LASTreg2[i])
  
  # Fig 3 results + MCMC samples
  if(i == 10){ # Results of Model NNNNN
    QuantilesNEW <- extract_Quantiles(Posterior, data = dataVal, AR = AR, K = label_K[i],
                                      last.regime = "N", quantiles = c(0.025,0.25,0.5,0.75,0.975))
    
    QuantilesNEW <- sapply(1:label_K[i], function(x) QuantilesNEW[[x]][,label_K[i]])
    Quantiles5 <- QuantilesNEW[idx2024,]
    colnames(Quantiles5) <- paste0("qNEW5",c("0.025","0.25","0.5","0.75", "0.975"))
    
    predictedValues <- extract_PredSummary(Posterior, dataVal, AR = AR,
                                           K = label_K[i], last.regime = label_LASTreg2[i])
    
    
    Quantiles1 <- predictedValues$mean_reg_quantiles[,idx2024,]
    Quantiles <- t(apply(Quantiles1, 2, c))
    colnames(Quantiles) <- paste0("q",rep(1:label_K[i],each=5), c("0.025","0.25","0.5","0.75", "0.975"))
    
  }
  
  # Temperature trajectories and generated regimes (In Sample)
  getTempTraj_insample <- TempIn(Posterior, dataTrain, AR = AR, comm_AR = comm_AR,
                                 comm_seas = comm_seas, comm_trend = comm_trend,
                                 K = label_K[i], last.regime = label_LASTreg2[i], Nsim = 200, proj = FALSE)
  
  TempPred[[i]] <- lapply(1: (year_end - year_start + 1), function (x){
    getTempTraj_insample[[x]]$gen$TempTrajectories$y_pred
  } )
  
  StatesIn[[i]] <- lapply(1: (year_end - year_start + 1), function (x){
    getTempTraj_insample[[x]]$gen$TempTrajectories$hat_regime
  } )
  
  # Temperature trajectories and generated regimes (Out of Sample)
  getTempTraj_OOS <- TempOut(Posterior, dataTrain, dataTest = dataVal, AR = AR, comm_AR = comm_AR,
                             comm_seas = comm_seas, comm_trend = comm_trend,
                             K = label_K[i], last.regime = label_LASTreg2[i], Nsim = 200, proj = TRUE)
  TempTraj[[i]] <- lapply(1: (year_end2 - year_end), function (x){
    getTempTraj_OOS[[x]]$gen$TempTrajectories$y_pred
  } )
  
  StatesOut[[i]]<-lapply(1: (year_end2 - year_end), function (x){
    getTempTraj_OOS[[x]]$gen$TempTrajectories$hat_regime
  } )
  
  print(i)
}

names(loo) <- names(log_lik) <- model_names
names(TempTraj) <-  model_names
names(TempPred) <-  model_names
names(StatesIn) <-  model_names
names(StatesOut) <-  model_names
names(HatProbIn) <- model_names
names(Exc) <- model_names

# Table 1 - T2
save(loo, file = "ExtractedResults/looTrieste_up2018_Bayes.Rda")

# Extract (200) temperature trajectories and generated states (Out of sample)
# Figure 4-5 and S8-S11
outTraj <- lapply(1:length(TempTraj), function(x){
  lapply(1:length(TempTraj[[1]]), function(z) TempTraj[[x]][[z]][, idx])
})
save(outTraj, file = "ExtractedResults/TempTrajTrieste_from2019_to2024_Bayes_all.Rda")

# Extract temperature trajectories (last regime N, Out of Sample)
# Figure 3 - S4-S5-S6-S7
TempTraj <- TempTraj[idx_model]
save(TempTraj, file = "ExtractedResults/TempTrajTrieste_from2019_to2024_Bayes.Rda")

save(StatesOut, file = "ExtractedResults/StatesOutTrieste_from2019_to2024_Bayes.Rda")

# Extract (200) temperature trajectories and generated states (In sample)
# Figure 4-5 and S8-S11
outTraj <- lapply(1:length(TempPred), function(x){
  lapply(1:length(TempPred[[1]]), function(z) TempPred[[x]][[z]][, idx])
})

outStates <- lapply(1:length(TempPred), function(x){
  lapply(1:length(TempPred[[1]]), function(z) StatesIn[[x]][[z]][, idx])
})

save(outTraj, file = "ExtractedResults/TempPredTrieste_up2018_Bayes.Rda")
save(outStates, file = "ExtractedResults/StatesInTrieste_up2018_Bayes.Rda")

# Fig 2
HatProbIn <- HatProbIn[c(1,4,7,10)]
save(HatProbIn, file = "ExtractedResults/HatProbInTrieste_up2018_Bayes.Rda")
save(Exc, file = "ExtractedResults/ExcTrieste_up2018_Bayes.Rda")

# Figure 3-S7
save(Quantiles5, file = "ExtractedResults/Quantiles5Trieste_from2019_to2024_Bayes.Rda")
save(Quantiles, file = "ExtractedResults/QuantilesTrieste_from2019_to2024_Bayes.Rda" )

