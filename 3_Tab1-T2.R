### Table 1 and T2 - LOOIC
maxK <- 5
label_LASTreg <- c("N", "G", "rW")

source("functions/load_Packages.R")

# Capriva ----
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

# Trieste ----
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

