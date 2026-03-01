# This file fits the MS models 
# Last update: 20/01/2025

# Set the working directory before running the code
source("functions/load_Packages.R")
source("functions/fitSTAN.R")



# Model formulas

foo_N_N_K2 <- list(mu_1 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_1 ~ 1 | dist ~ Gaussian,
                   mu_2 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_2 ~ 1 | dist ~ Gaussian,
                   Tp_1.2 ~ 1, Tp_2.1 ~ 1)
foo_N_G_K2 <- list(mu_1 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_1 ~ 1 | dist ~ Gaussian,
                   mu_2 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_2 ~ 1 | dist ~ Gumbel,
                   Tp_1.2 ~ 1, Tp_2.1 ~ 1)
foo_N_rW_K2 <- list(mu_1 ~  GST_c + s(sc_doy, bs = "bs", k=5)|  sc_1 ~ 1 | dist ~ Gaussian,
                    mu_2 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_2 ~ 1 | sh_2 ~ 1 | dist ~ rWeibull,
                    Tp_1.2 ~ 1, Tp_2.1 ~ 1)

foo_N_N_K3 <- list(mu_1 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_1 ~ 1 | dist ~ Gaussian,
                   mu_2 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_2 ~ 1 | dist ~ Gaussian,
                   mu_3 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_3 ~ 1 | dist ~ Gaussian,
                   Tp_1.2 ~ 1, Tp_2.1 ~ 1, Tp_2.3 ~ 1, Tp_3.2 ~ 1)
foo_N_G_K3 <- list(mu_1 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_1 ~ 1 | dist ~ Gaussian,
                   mu_2 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_2 ~ 1 | dist ~ Gaussian,
                   mu_3 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_3 ~ 1 | dist ~ Gumbel,
                   Tp_1.2 ~ 1, Tp_2.1 ~ 1, Tp_2.3 ~ 1, Tp_3.2 ~ 1)
foo_N_rW_K3 <- list(mu_1 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_1 ~ 1 | dist ~ Gaussian,
                    mu_2 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_2 ~ 1 | dist ~ Gaussian,
                    mu_3 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_3 ~ 1 | sh_3 ~ 1 | dist ~ rWeibull,
                    Tp_1.2 ~ 1, Tp_2.1 ~ 1, Tp_2.3 ~ 1, Tp_3.2 ~ 1)

foo_N_N_K4 <- list(mu_1 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_1 ~ 1 | dist ~ Gaussian,
                   mu_2 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_2 ~ 1 | dist ~ Gaussian,
                   mu_3 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_3 ~ 1 | dist ~ Gaussian,
                   mu_4 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_4 ~ 1 | dist ~ Gaussian,
                   Tp_1.2 ~ 1, Tp_2.1 ~ 1, Tp_2.3 ~ 1, Tp_3.2 ~ 1, Tp_3.4 ~ 1, Tp_4.3 ~ 1)
foo_N_G_K4 <- list(mu_1 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_1 ~ 1 | dist ~ Gaussian,
                   mu_2 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_2 ~ 1 | dist ~ Gaussian,
                   mu_3 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_3 ~ 1 | dist ~ Gaussian,
                   mu_4 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_4 ~ 1 | dist ~ Gumbel,
                   Tp_1.2 ~ 1, Tp_2.1 ~ 1, Tp_2.3 ~ 1, Tp_3.2 ~ 1, Tp_3.4 ~ 1, Tp_4.3 ~ 1)
foo_N_rW_K4 <- list(mu_1 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_1 ~ 1 | dist ~ Gaussian,
                    mu_2 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_2 ~ 1 | dist ~ Gaussian,
                    mu_3 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_3 ~ 1 | dist ~ Gaussian,
                    mu_4 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_4 ~ 1 | sh_4 ~ 1 | dist ~ rWeibull,
                    Tp_1.2 ~ 1, Tp_2.1 ~ 1, Tp_2.3 ~ 1, Tp_3.2 ~ 1, Tp_3.4 ~ 1, Tp_4.3 ~ 1)

foo_N_N_K5 <- list(mu_1 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_1 ~ 1 | dist ~ Gaussian,
                   mu_2 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_2 ~ 1 | dist ~ Gaussian,
                   mu_3 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_3 ~ 1 | dist ~ Gaussian,
                   mu_4 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_4 ~ 1 | dist ~ Gaussian,
                   mu_5 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_5 ~ 1 | dist ~ Gaussian,
                   Tp_1.2 ~ 1, Tp_2.1 ~ 1, Tp_2.3 ~ 1, Tp_3.2 ~ 1, Tp_3.4 ~ 1, Tp_4.3 ~ 1, Tp_4.5 ~ 1, Tp_5.4 ~ 1)
foo_N_G_K5 <- list(mu_1 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_1 ~ 1 | dist ~ Gaussian,
                   mu_2 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_2 ~ 1 | dist ~ Gaussian,
                   mu_3 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_3 ~ 1 | dist ~ Gaussian,
                   mu_4 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_4 ~ 1 | dist ~ Gaussian,
                   mu_5 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_5 ~ 1 | dist ~ Gumbel,
                   Tp_1.2 ~ 1, Tp_2.1 ~ 1, Tp_2.3 ~ 1, Tp_3.2 ~ 1, Tp_3.4 ~ 1, Tp_4.3 ~ 1, Tp_4.5 ~ 1, Tp_5.4 ~ 1)
foo_N_rW_K5 <- list(mu_1 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_1 ~ 1 | dist ~ Gaussian,
                    mu_2 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_2 ~ 1 | dist ~ Gaussian,
                    mu_3 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_3 ~ 1 | dist ~ Gaussian,
                    mu_4 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_4 ~ 1 | dist ~ Gaussian,
                    mu_5 ~  GST_c + s(sc_doy, bs = "bs", k=5)| sc_5 ~ 1 | sh_5 ~ 1 | dist ~ rWeibull,
                    Tp_1.2 ~ 1, Tp_2.1 ~ 1, Tp_2.3 ~ 1, Tp_3.2 ~ 1, Tp_3.4 ~ 1, Tp_4.3 ~ 1, Tp_4.5 ~ 1, Tp_5.4 ~ 1)


## TRIESTE ----
location <- "Trieste"
load(paste0("data/dataFVG/dati", location, ".Rdata"))
dataT <- dati[dati$citta == location,] 

## K2 ---- 

# Normal - Normal
fit_NN_T <- model_fitting(data = dataT, modfoo = foo_N_N_K2, trainStart = 1995, trainEnd = 2018, regimes = "N-N",
                          K = 2, iter = 2000, initialise_param = "custom",
                          comm_seas = TRUE, comm_trend = F,
                          piStart = c(0.5,0.5))
save(fit_NN_T, file = "ResultsGST/fit_NN_T_s.RData")
rm(fit_NN_T)

# Normal - Gumbel
fit_NG_T <- model_fitting(data = dataT, modfoo = foo_N_G_K2, trainStart = 1995, trainEnd = 2018, regimes = "N-G",
                          K = 2, iter = 2000, initialise_param = "custom",
                          comm_seas = TRUE, comm_trend = F,
                          piStart = c(0.5,0.5))
save(fit_NG_T, file = "ResultsGST/fit_NG_T_s.RData")
rm(fit_NG_T)

# Normal - rWeibull
fit_NrW_T <- model_fitting(data = dataT, modfoo = foo_N_rW_K2, trainStart = 1995, trainEnd = 2018, regimes = "N-rW",
                           K = 2, iter = 2000, initialise_param = "custom",
                           comm_seas = TRUE, comm_trend = F,
                           piStart = c(0.5,0.5))
save(fit_NrW_T, file = "ResultsGST/fit_NrW_T_s.RData")
rm(fit_NrW_T)

## K3 ----
# Normal - Normal
fit_NNN_T <- model_fitting(data = dataT, modfoo = foo_N_N_K3, trainStart = 1995, trainEnd = 2018, regimes = "N-N",
                           K = 3, iter = 2000, initialise_param = "custom",
                           comm_seas = TRUE, comm_trend = F,
                           piStart = c(0.5,rep(0.5/(2),2)))
save(fit_NNN_T, file = "ResultsGST/fit_NNN_T_s.RData")
rm(fit_NNN_T)

# Normal - Gumbel
fit_NNG_T <-  model_fitting(data = dataT, modfoo = foo_N_G_K3,  trainStart = 1995, trainEnd = 2018,
                            regimes = "N-G", K = 3, iter = 2000, 
                            piStart = c(0.5,rep(0.5/(2),2)),
                            initialise_param = "custom", comm_seas = T, comm_trend = F)
save(fit_NNG_T, file = "ResultsGST/fit_NNG_T_s.RData")
rm(fit_NNG_T)

# Normal - rWeibull
fit_NNrW_T <- model_fitting(data = dataT, modfoo = foo_N_rW_K3, trainStart = 1995, trainEnd = 2018,
                            regimes = "N-rW", K = 3, iter = 2000, 
                            piStart = c(0.5,rep(0.5/(2),2)),
                            initialise_param = "custom", comm_seas = T, comm_trend = F)
save(fit_NNrW_T, file = "ResultsGST/fit_NNrW_T_s.RData")
rm(fit_NNrW_T)

## K4 ----
# Normal - Normal
fit_NNNN_T <- model_fitting(data = dataT, modfoo = foo_N_N_K4, trainStart = 1995, trainEnd = 2018, regimes = "N-N",
                           K = 4, iter = 2000, initialise_param = "custom",
                           comm_seas = TRUE, comm_trend = F,
                           piStart = c(0.5,rep(0.5/(3),3)))
save(fit_NNNN_T, file = "ResultsGST/fit_NNNN_T_s.RData")
rm(fit_NNNN_T)

# Normal - Gumbel
fit_NNNG_T <-  model_fitting(data = dataT, modfoo = foo_N_G_K4,  trainStart = 1995, trainEnd = 2018,
                            regimes = "N-G", K = 4, iter = 2000, 
                            piStart = c(0.5,rep(0.5/(3),3)),
                            initialise_param = "custom", comm_seas = T, comm_trend = F)
save(fit_NNNG_T, file = "ResultsGST/fit_NNNG_T_s.RData")
rm(fit_NNNG_T)

# Normal - rWeibull
fit_NNNrW_T <- model_fitting(data = dataT, modfoo = foo_N_rW_K4, trainStart = 1995, trainEnd = 2018,
                            regimes = "N-rW", K = 4, iter = 2000, 
                            piStart = c(0.5,rep(0.5/(3),3)),
                            initialise_param = "custom", comm_seas = T, comm_trend = F)
save(fit_NNNrW_T, file = "ResultsGST/fit_NNNrW_T_s.RData")
rm(fit_NNNrW_T)

## K5 ----
# Normal - Normal
fit_NNNNN_T <- model_fitting(data = dataT, modfoo = foo_N_N_K5, trainStart = 1995, trainEnd = 2018, regimes = "N-N", 
                             K = 5, iter = 2000, initialise_param = "custom", 
                             comm_seas = TRUE, comm_trend = F, 
                             piStart = c(0.5,rep(0.5/(4),4))) 
save(fit_NNNNN_T, file = "ResultsGST/fit_NNNNN_T_s.RData")
rm(fit_NNNNN_T)

fit_stan

# Normal - Gumbel 
fit_NNNNG_T <-  model_fitting(data = dataT, modfoo = foo_N_G_K5,  trainStart = 1995, trainEnd = 2018, 
                              regimes = "N-G", K = 5, iter = 2000, 
                              piStart = c(0.5,rep(0.5/(4),4)),
                              initialise_param = "MAP", comm_seas = T, comm_trend = F) 
save(fit_NNNNG_T, file = "ResultsGST/fit_NNNNG_T_s.RData")
rm(fit_NNNNG_T)

# Normal - rWeibull 
fit_NNNNrW_T <- model_fitting(data = dataT, modfoo = foo_N_rW_K5, trainStart = 1995, trainEnd = 2018, 
                              regimes = "N-rW", K = 5, iter = 2000, 
                              piStart = c(0.5,rep(0.5/(4),4)),
                              initialise_param = "custom", comm_seas = T, comm_trend = F)
save(fit_NNNNrW_T, file = "ResultsGST/fit_NNNNrW_T_s.RData")
rm(fit_NNNNrW_T)

## CAPRIVA ----
location <- "Capriva"
load(paste0("data/dataFVG/dati", location, ".Rdata"))
dataC <- dati[dati$citta == location,] 

## K2 ----
# Normal - Normal
fit_NN_C <- model_fitting(data = dataC, modfoo = foo_N_N_K2, trainStart = 1995, trainEnd = 2018, regimes = "N-N",
                          K = 2, iter = 2000, initialise_param = "custom",
                          comm_seas = TRUE, comm_trend = F,
                          piStart = c(0.5,0.5))
save(fit_NN_C, file = "ResultsGST/fit_NN_C_s.RData")
rm(fit_NN_C)

# Normal - Gumbel
fit_NG_C <- model_fitting(data = dataC, modfoo = foo_N_G_K2, trainStart = 1995, trainEnd = 2018, regimes = "N-G",
                          K = 2, iter = 2000, initialise_param = "custom",
                          comm_seas = TRUE, comm_trend = F,
                          piStart = c(0.5,0.5))
save(fit_NG_C, file = "ResultsGST/fit_NG_C_s.RData")
rm(fit_NG_C)

# Normal - rWeibull
fit_NrW_C <- model_fitting(data = dataC, modfoo = foo_N_rW_K2, trainStart = 1995, trainEnd = 2018, regimes = "N-rW",
                           K = 2, iter = 2000, initialise_param = "custom",
                           comm_seas = TRUE, comm_trend = F,
                           piStart = c(0.5,0.5))
save(fit_NrW_C, file = "ResultsGST/fit_NrW_C_s.RData")
rm(fit_NrW_C)

## K3 ----
# Normal - Normal
fit_NNN_C <- model_fitting(data = dataC, modfoo = foo_N_N_K3, trainStart = 1995, trainEnd = 2018, regimes = "N-N",
                           K = 3, iter = 2000, initialise_param = "custom",
                           comm_seas = TRUE, comm_trend = F,
                           piStart = c(0.5,rep(0.5/(2),2)))
save(fit_NNN_C, file = "ResultsGST/fit_NNN_C_s.RData")
rm(fit_NNN_C)

# Normal - Gumbel
fit_NNG_C <-  model_fitting(data = dataC, modfoo = foo_N_G_K3,  trainStart = 1995, trainEnd = 2018,
                            regimes = "N-G", K = 3, iter = 2000, piStart = c(0.5,rep(0.5/(2),2)),
                            initialise_param = "custom", comm_seas = T, comm_trend = F)
save(fit_NNG_C, file = "ResultsGST/fit_NNG_C_s.RData")
rm(fit_NNG_C)

# Normal - rWeibull
fit_NNrW_C <- model_fitting(data = dataC, modfoo = foo_N_rW_K3, trainStart = 1995, trainEnd = 2018,
                            regimes = "N-rW", K = 3, iter = 2000, piStart = c(0.5,rep(0.5/(2),2)),
                            initialise_param = "custom", comm_seas = T, comm_trend = F)
save(fit_NNrW_C, file = "ResultsGST/fit_NNrW_C_s.RData")
rm(fit_NNrW_C)

## K4 ----
# Normal - Normal
fit_NNNN_C <- model_fitting(data = dataC, modfoo = foo_N_N_K4, trainStart = 1995, trainEnd = 2018, regimes = "N-N",
                            K = 4, iter = 2000, initialise_param = "custom",
                            comm_seas = TRUE, comm_trend = F,
                            piStart = c(0.5,rep(0.5/(3),3)))
save(fit_NNNN_C, file = "ResultsGST/fit_NNNN_C_s.RData")
rm(fit_NNNN_C)

# Normal - Gumbel
fit_NNNG_C <-  model_fitting(data = dataC, modfoo = foo_N_G_K4,  trainStart = 1995, trainEnd = 2018,
                             regimes = "N-G", K = 4, iter = 2000, piStart = c(0.5,rep(0.5/(3),3)),
                             initialise_param = "custom", comm_seas = T, comm_trend = F)
save(fit_NNNG_C, file = "ResultsGST/fit_NNNG_C_s.RData")
rm(fit_NNNG_C)

# Normal - rWeibull
fit_NNNrW_C <- model_fitting(data = dataC, modfoo = foo_N_rW_K4, trainStart = 1995, trainEnd = 2018,
                             regimes = "N-rW", K = 4, iter = 2000, piStart = c(0.5,rep(0.5/(3),3)),
                             initialise_param = "custom", comm_seas = T, comm_trend = F)
save(fit_NNNrW_C, file = "ResultsGST/fit_NNNrW_C_s.RData")
rm(fit_NNNrW_C)

## K5 ----
# Normal - Normal
fit_NNNNN_C <- model_fitting(data = dataC, modfoo = foo_N_N_K5, trainStart = 1995, trainEnd = 2018, regimes = "N-N",
                             K = 5, iter = 2000, initialise_param = "custom",
                             comm_seas = TRUE, comm_trend = F,
                             piStart = c(0.5,rep(0.5/(4),4)))
save(fit_NNNNN_C, file = "ResultsGST/fit_NNNNN_C_s.RData")
rm(fit_NNNNN_C)

# Normal - Gumbel
fit_NNNNG_C <-  model_fitting(data = dataC, modfoo = foo_N_G_K5,  trainStart = 1995, trainEnd = 2018,
                              regimes = "N-G", K = 5, iter = 2000, piStart = c(0.5,rep(0.5/(4),4)),
                              initialise_param = "MAP", comm_seas = T, comm_trend = F)
save(fit_NNNNG_C, file = "ResultsGST/fit_NNNNG_C_s.RData")
rm(fit_NNNNG_C)

# Normal - rWeibull
fit_NNNNrW_C <- model_fitting(data = dataC, modfoo = foo_N_rW_K5, trainStart = 1995, trainEnd = 2018,
                              regimes = "N-rW", K = 5, iter = 2000, piStart = c(0.5,rep(0.5/(4),4)),
                              initialise_param = "custom", comm_seas = T, comm_trend = F)
save(fit_NNNNrW_C, file = "ResultsGST/fit_NNNNrW_C_s.RData")
rm(fit_NNNNrW_C)

