# This code extracts the results analyzed in the paper for the fitted models:
# - loo 
# - predicted temperatures (in and out of sample)
# - predicted states (in and out of sample)
# - predicted probabilities (in and out of sample)

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

## Capriva----

location <- "Capriva"
stations <- "_C"

comm_seas <- TRUE
comm_trend <- FALSE
AR <- FALSE
comm_AR <- FALSE

load("dati/dataFVG/datiCapriva.Rdata")

year_seq <- 1995:2024

model_names <- paste0("fit_", as.vector(t(outer(label_reg, label_LASTreg, FUN = paste0))), stations, "_s")
model_names2 <- paste0("fit_", as.vector(t(outer(label_reg, label_LASTreg, FUN = paste0))), stations)
file_select <- paste0(model_names,  ".RData") # "_K", label_K,
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
ProbFore <- list()
ProbFilt <- list()

for(i in get_model){

  # Load Data
  load(paste0("ResultsGST/",file_select[i]))
  fit <- mget(model_names2[i], envir = globalenv())[[1]] 
  rm(list=grep(model_names2[i],ls(),value=TRUE))
  
  # Get log-lik and loo
  log_lik[[i]] <- extract_log_lik(fit)
  loo[[i]] <- loo(log_lik[[i]])
  
  # Get Posterior
  Posterior <- ExtractPost(stan.obj = fit, K = label_K[i], last.regime = label_LASTreg2[i], 
                           comm_seas = comm_seas, comm_trend = comm_trend)

  
  # Temperature trajectories on IS 
  getTempTraj_insample <- TempIn(Posterior, dataTrain, AR = AR, comm_AR = comm_AR, 
                                 comm_seas = comm_seas, comm_trend = comm_trend, 
                                 K = label_K[i], last.regime = label_LASTreg2[i], Nsim = 200, proj = FALSE)
  

  # Temperature trajectories on OOS 
  getTempTraj_OOS <- TempOut(Posterior, dataTrain, dataTest = dataVal, AR = AR, comm_AR = comm_AR, 
                            comm_seas = comm_seas, comm_trend = comm_trend, 
                            K = label_K[i], last.regime = label_LASTreg2[i], Nsim = 200, proj = TRUE)

  
  TempPred[[i]] <- lapply(1: (year_end - year_start + 1), function (x){
    getTempTraj_insample[[x]]$gen$TempTrajectories$y_pred
  } )
  
  StatesIn[[i]] <- lapply(1: (year_end - year_start + 1), function (x){
    getTempTraj_insample[[x]]$gen$TempTrajectories$hat_regime
  } )
  
  ProbFilt[[i]] <- lapply(1: (year_end - year_start + 1), function (x){
    getTempTraj_insample[[x]]$gen$TempTrajectories$ProbFor
  } )
  
  TempTraj[[i]] <- lapply(1: (year_end2 - year_end), function (x){
    getTempTraj_OOS[[x]]$gen$TempTrajectories$y_pred
  } )
  
  StatesOut[[i]]<-lapply(1: (year_end2 - year_end), function (x){
    getTempTraj_OOS[[x]]$gen$TempTrajectories$hat_regime
  } )

   ProbFore[[i]] <- lapply(1: (year_end2 - year_end), function (x){
     getTempTraj_OOS[[x]]$gen$TempTrajectories$ProbFor
  } )
  
  print(i)
}

names(loo) <- names(log_lik) <- model_names

names(TempTraj) <-  model_names
names(TempPred) <-  model_names
names(StatesIn) <-  model_names
names(StatesOut) <-  model_names
names(ProbFore) <-  model_names
names(ProbFilt) <-  model_names

save(loo, file = "ExtractedResults/looCapriva_up2018_Bayes.Rda")

save(TempTraj, file = "ExtractedResults/TempTrajCapriva_from2019_to2024_Bayes.Rda")
save(TempPred, file = "ExtractedResults/TempPredCapriva_up2018_Bayes.Rda")

save(StatesIn, file = "ExtractedResults/StatesInCapriva_up2018_Bayes.Rda")
save(StatesOut, file = "ExtractedResults/StatesOutCapriva_from2019_to2024_Bayes.Rda")

save(ProbFilt, file = "ExtractedResults/ProbFiltCapriva_up2018_Bayes.Rda")
save(ProbFore, file = "ExtractedResults/ProbForeCapriva_from2019_to2024_Bayes.Rda")

# Get loo_Capriva
load("ExtractedResults/looCapriva_up2018_Bayes.Rda")
loo_table_Capriva <- matrix(unlist(lapply(1:length(loo), function(x) loo[[x]]$estimates[3,1])), 
                            nrow = maxK - 1, ncol = length(label_LASTreg), byrow = TRUE)
round(loo_table_Capriva,1)
#        [,1]    [,2]    [,3]
# [1,] 18326.3 18316.1 18284.7
# [2,] 17659.3 17621.3 17617.8
# [3,] 17271.4 17255.6 17253.8
# [4,] 17117.6 17109.5 17110.3

# regime_comparison
idx <-lapply(1:(maxK-1), function(x) (3*(x-1) + 1):(3*x))
looComp <- lapply(1:length(idx), function(x) loo_compare(loo[idx[[x]]]))
looComp
# [[1]]
# elpd_diff se_diff
# fit_NrW_C_s   0.0       0.0  
# fit_NG_C_s  -15.7       5.4  
# fit_NN_C_s  -20.8       6.3  
# 
# [[2]]
# elpd_diff se_diff
# fit_NNrW_C_s   0.0       0.0  
# fit_NNG_C_s   -1.8       2.6  
# fit_NNN_C_s  -20.7       7.3  
# 
# [[3]]
# elpd_diff se_diff
# fit_NNNrW_C_s  0.0       0.0   
# fit_NNNG_C_s  -0.9       2.6   
# fit_NNNN_C_s  -8.8       4.8   
# 
# [[4]]
# elpd_diff se_diff
# fit_NNNNG_C_s   0.0       0.0   
# fit_NNNNrW_C_s -0.4       2.6   
# fit_NNNNN_C_s  -4.0       6.0   

# K comparison 
idx <- lapply(1:3, function(x) seq(x, 3*(maxK-1), 3))
looComp <- lapply(1:length(idx), function(x) loo_compare(loo[idx[[x]]]))
looComp
# [[1]]
# elpd_diff se_diff
# fit_NNNNN_C_s    0.0       0.0 
# fit_NNNN_C_s   -76.9      18.6 
# fit_NNN_C_s   -270.8      26.1 
# fit_NN_C_s    -604.3      38.7 
# 
# [[2]]
# elpd_diff se_diff
# fit_NNNNG_C_s    0.0       0.0 
# fit_NNNG_C_s   -73.0      18.1 
# fit_NNG_C_s   -255.9      25.7 
# fit_NG_C_s    -603.3      39.7 
# 
# [[3]]
# elpd_diff se_diff
# fit_NNNNrW_C_s    0.0       0.0 
# fit_NNNrW_C_s   -71.7      18.3 
# fit_NNrW_C_s   -253.7      25.4 
# fit_NrW_C_s    -587.2      38.8 

## Trieste----

location <- "Trieste"
stations <- "_T"

comm_seas <- TRUE
comm_trend <- FALSE
AR <- FALSE
comm_AR <- FALSE

load("dati/dataFVG/datiTrieste.Rdata")

year_seq <- 1995:2024

model_names <- paste0("fit_", as.vector(t(outer(label_reg, label_LASTreg, FUN = paste0))), stations, "_s")
model_names2 <- paste0("fit_", as.vector(t(outer(label_reg, label_LASTreg, FUN = paste0))), stations)
file_select <- paste0(model_names,  ".RData") 
fileName <- list.files(path = "../ResultsGST", pattern = stations)
get_model <- which(file_select%in%fileName)

# Train and  Validation 
dataTrain <- dati[dati$Year <= year_end &  dati$Month >= 5 & dati$Month <= 9,]
dataVal <- dati[dati$Year > year_end & dati$Year <= year_end2 & dati$Month >= 5 & dati$Month <= 9,]

log_lik <- loo <- list()

TempPred <- list()
StatesIn <- list()
TempTraj <- list()
StatesOut <- list()
ProbFore <- list()
ProbFilt <- list()

for(i in get_model){
  
  # Load Data
  load(paste0("ResultsGST/",file_select[i]))
  fit <- mget(model_names2[i], envir = globalenv())[[1]] 
  rm(list=grep(model_names2[i],ls(),value=TRUE))
  
  # Get log-lik and loo
  log_lik[[i]] <- extract_log_lik(fit)
  loo[[i]] <- loo(log_lik[[i]])
  
  # Get Posterior
  Posterior <- ExtractPost(stan.obj = fit, K = label_K[i], last.regime = label_LASTreg2[i], 
                           comm_seas = comm_seas, comm_trend = comm_trend)
  
  # Temperature trajectories on IS 
  getTempTraj_insample <- TempIn(Posterior, dataTrain, AR = AR, comm_AR = comm_AR, 
                                 comm_seas = comm_seas, comm_trend = comm_trend, 
                                 K = label_K[i], last.regime = label_LASTreg2[i], 
                                 Nsim = 200, proj = FALSE)

  # Temperature trajectories on OOS 
  getTempTraj_OOS <- TempOut(Posterior, dataTrain, dataTest = dataVal, AR = AR, comm_AR = comm_AR, 
                            comm_seas = comm_seas, comm_trend = comm_trend, 
                            K = label_K[i], last.regime = label_LASTreg2[i], Nsim = 200, proj = FALSE)
  
  
  TempPred[[i]] <- lapply(1: (year_end - year_start + 1), function (x){
    getTempTraj_insample[[x]]$gen$TempTrajectories$y_pred
  } )
  
  StatesIn[[i]] <- lapply(1: (year_end - year_start + 1), function (x){
    getTempTraj_insample[[x]]$gen$TempTrajectories$hat_regime
  } )
  
  ProbFilt[[i]] <- lapply(1: (year_end - year_start + 1), function (x){
    getTempTraj_insample[[x]]$gen$TempTrajectories$ProbFor
  } )
  
  TempTraj[[i]] <- lapply(1: (year_end2 - year_end), function (x){
    getTempTraj_OOS[[x]]$gen$TempTrajectories$y_pred
  } )
  
  StatesOut[[i]]<-lapply(1: (year_end2 - year_end), function (x){
    getTempTraj_OOS[[x]]$gen$TempTrajectories$hat_regime
  } )
  
  ProbFore[[i]] <- lapply(1: (year_end2 - year_end), function (x){
    getTempTraj_OOS[[x]]$gen$TempTrajectories$ProbFor
  } )
  
  print(i)
}

names(loo) <- names(log_lik) <- model_names

names(TempTraj) <-  model_names
names(TempPred) <-  model_names

names(StatesOut) <-  model_names
names(StatesIn) <-  model_names

names(ProbFore) <-  model_names
names(ProbFilt) <-  model_names

save(loo, file = "ExtractedResults/looTrieste_up2018_Bayes.Rda")

save(TempTraj, file = "ExtractedResults/TempTrajTrieste_from2019_to2024_Bayes.Rda")
save(TempPred, file = "ExtractedResults/TempPredTrieste_up2018_Bayes.Rda")

save(StatesIn, file = "ExtractedResults/StatesInTrieste_up2018_Bayes.Rda")
save(StatesOut, file = "ExtractedResults/StatesOutTrieste_from2019_to2024_Bayes.Rda")

save(ProbFilt, file = "ExtractedResults/ProbFiltTrieste_up2018_Bayes.Rda")
save(ProbFore, file = "ExtractedResults/ProbForeTrieste_from2019_to2024_Bayes.Rda")

# Get loo_Trieste
load("ExtractedResults/looTrieste_up2018_Bayes.Rda")
loo_table_Trieste <- matrix(unlist(lapply(1:length(loo), function(x) loo[[x]]$estimates[3,1])), 
                            nrow = maxK - 1, ncol = length(label_LASTreg), byrow = TRUE)

round(loo_table_Trieste,1)
# [,1]    [,2]    [,3]
# [1,] 16278.7 16205.5 16185.8
# [2,] 15546.3 15520.0 15518.3
# [3,] 15163.7 15167.2 15157.0
# [4,] 15065.1 15089.2 15059.1

# regime_comparison
idx <- lapply(1:(maxK-1), function(x) (3*(x-1) + 1):(3*x))
looComp <- lapply(1:length(idx), function(x) loo_compare(loo[idx[[x]]]))
looComp
# [[1]]
# elpd_diff se_diff
# fit_NrW_T_s   0.0       0.0  
# fit_NG_T_s   -9.9       4.0  
# fit_NN_T_s  -46.5      10.3  
# 
# [[2]]
# elpd_diff se_diff
# fit_NNrW_T_s   0.0       0.0  
# fit_NNG_T_s   -0.9       3.3  
# fit_NNN_T_s  -14.0      13.7  
# 
# [[3]]
# elpd_diff se_diff
# fit_NNNrW_T_s  0.0       0.0   
# fit_NNNN_T_s  -3.3       3.5   
# fit_NNNG_T_s  -5.1       3.0   
# 
# [[4]]
# elpd_diff se_diff
# fit_NNNNrW_T_s   0.0       0.0  
# fit_NNNNN_T_s   -3.0       3.4  
# fit_NNNNG_T_s  -15.0       3.7   

# K comparison 
idx <- lapply(1:3, function(x) seq(x, 3*(maxK-1), 3))
looComp <- lapply(1:length(idx), function(x) loo_compare(loo[idx[[x]]]))
looComp
# [[1]]
# elpd_diff se_diff
# fit_NNNNN_T_s    0.0       0.0 
# fit_NNNN_T_s   -49.3      12.2 
# fit_NNN_T_s   -240.6      24.9 
# fit_NN_T_s    -606.8      35.4 
# 
# [[2]]
# elpd_diff se_diff
# fit_NNNNG_T_s    0.0       0.0 
# fit_NNNG_T_s   -39.0      12.6 
# fit_NNG_T_s   -215.4      23.6 
# fit_NG_T_s    -558.2      34.6 
# 
# [[3]]
# elpd_diff se_diff
# fit_NNNNrW_T_s    0.0       0.0 
# fit_NNNrW_T_s   -48.9      12.4 
# fit_NNrW_T_s   -229.6      24.0 
# fit_NrW_T_s    -563.3      34.7 
