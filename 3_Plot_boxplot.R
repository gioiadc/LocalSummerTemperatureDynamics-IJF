# This code reproduce the boxplot presented in the main manuscript and in
# Supplementary Material

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

load("HeatwaveR/extract_HWsw_summary_C.Rdata")
load("HeatwaveR/extract_HWsw_summary_Cwith2024.Rdata")

load("HeatwaveR/extract_HWsw_summary_T.Rdata")
load("HeatwaveR/extract_HWsw_summary_Twith2024.Rdata")

load("ExtractedResults/database_Hwr.Rdata")

comm_seas <- TRUE
comm_trend <- FALSE
comm_AR <- FALSE
AR <- FALSE

idx <- sample(1:4000, 200, replace = F)

# Out of sample ----
## Capriva ----
stations <- "_C"
model_names2 <- paste0("fit_", as.vector(t(outer(label_reg, label_LASTreg, FUN = paste0))), stations)

load("ExtractedResults/TempTrajCapriva_from2019_to2024_Bayes.Rda")
load("ExtractedResults/StatesOutCapriva_from2019_to2024_Bayes.Rda")

outTraj <- lapply(1:length(TempTraj), function(x){
  lapply(1:length(TempTraj[[1]]), function(z) TempTraj[[x]][[z]][, idx])
})

outStates <- lapply(1:length(TempTraj), function(x){
  lapply(1:length(TempTraj[[1]]), function(z) StatesOut[[x]][[z]][, idx])
})


foreHWdist_Capriva <- getForeHWdist2(trajectories = outTraj, 
                                     states =  outStates,
                                     year_start = 2019, 
                                     year_end = 2024, 
                                     AR = AR, 
                                     model_names = model_names2,
                                     K = label_K, 
                                     last.regime = label_LASTreg2,
                                     site = "Capriva", 
                                     Nsim = 200,
                                     year_thresh = 2019) 
HWdist_CaprivaS <- do.call("rbind",lapply(1:length(TempTraj), function(x) do.call("rbind", foreHWdist_Capriva[[x]])))
HWdist_CaprivaS$method <- "State"


foreHWdist_Capriva <- getForeHWdist(trajectories = outTraj, 
                                     dataHWr = extract_HWsw_summary_C3$res, 
                                     year_start = 2019, 
                                     year_end = 2024, 
                                     AR = AR, 
                                     model_names = model_names2,
                                     K = label_K, 
                                     last.regime = label_LASTreg2,
                                     site = "Capriva", 
                                     Nsim = 200,
                                     year_thresh = 2019) 
HWdist_CaprivaT <- do.call("rbind",lapply(1:length(TempTraj), function(x) do.call("rbind", foreHWdist_Capriva[[x]])))
HWdist_CaprivaT$method <- "Perk"

HWdist_Capriva <- rbind(HWdist_CaprivaS, HWdist_CaprivaT)

## Trieste ----
stations <- "_T"
model_names2 <- paste0("fit_", as.vector(t(outer(label_reg, label_LASTreg, FUN = paste0))), stations)

load("ExtractedResults/TempTrajTrieste_from2019_to2024_Bayes.Rda")
load("ExtractedResults/StatesOutTrieste_from2019_to2024_Bayes.Rda")

outTraj <- lapply(1:length(TempTraj), function(x){
  lapply(1:length(TempTraj[[1]]), function(z) TempTraj[[x]][[z]][, idx])
})

outStates <- lapply(1:length(TempTraj), function(x){
  lapply(1:length(TempTraj[[1]]), function(z) StatesOut[[x]][[z]][, idx])
})


foreHWdist_Trieste <- getForeHWdist2(trajectories = outTraj, 
                                     states =  outStates,
                                     year_start = 2019, 
                                     year_end = 2024, 
                                     AR = AR, 
                                     model_names = model_names2,
                                     K = label_K, 
                                     last.regime = label_LASTreg2,
                                     site = "Trieste", 
                                     Nsim = 200,
                                     year_thresh = 2019) 
HWdist_TriesteS <- do.call("rbind",lapply(1:length(TempTraj), function(x) do.call("rbind", foreHWdist_Trieste[[x]])))
HWdist_TriesteS$method <- "State"

foreHWdist_Trieste <- getForeHWdist(trajectories = outTraj, 
                                    dataHWr = extract_HWsw_summary_T3$res, 
                                    year_start = 2019, 
                                    year_end = 2024, 
                                    AR = AR, 
                                    model_names = model_names2,
                                    K = label_K, 
                                    last.regime = label_LASTreg2,
                                    site = "Trieste", 
                                    Nsim = 200,
                                    year_thresh = 2019) 
HWdist_TriesteT <- do.call("rbind",lapply(1:length(TempTraj), function(x) do.call("rbind", foreHWdist_Trieste[[x]])))
HWdist_TriesteT$method <- "Perk"

HWdist_Trieste <- rbind(HWdist_TriesteS, HWdist_TriesteT)

HWdist_out <- rbind.data.frame(HWdist_Capriva, HWdist_Trieste)



# In sample ----
## Capriva ----
load("ExtractedResults/TempPredCapriva_up2018_Bayes.Rda")
load("ExtractedResults/StatesInCapriva_up2018_Bayes.Rda")

outTraj <- lapply(1:length(TempPred), function(x){
  lapply(1:length(TempPred[[1]]), function(z) TempPred[[x]][[z]][, idx])
})

outStates <- lapply(1:length(TempPred), function(x){
  lapply(1:length(TempPred[[1]]), function(z) StatesIn[[x]][[z]][, idx])
})


foreHWdist_Capriva <- getForeHWdist2(trajectories = outTraj, 
                                     states = outStates,
                                     year_start = 1995, 
                                     year_end = 2018, 
                                     AR = AR, 
                                     model_names = model_names2,
                                     K = label_K, 
                                     last.regime = label_LASTreg2,
                                     site = "Capriva", 
                                     Nsim = 200,
                                     year_thresh = 1995) 
HWdist_CaprivaS<- do.call("rbind",lapply(1:length(TempPred), function(x) do.call("rbind", foreHWdist_Capriva[[x]]  )))
HWdist_CaprivaS$method <- "State"

foreHWdist_Capriva <- getForeHWdist(trajectories = outTraj, 
                                    dataHWr = extract_HWsw_summary_C$res, 
                                    year_start = 1995, 
                                    year_end = 2018, 
                                    AR = AR, 
                                    model_names = model_names2,
                                    K = label_K, 
                                    last.regime = label_LASTreg2,
                                    site = "Capriva", 
                                    Nsim = 200,
                                    year_thresh = 1995) 
HWdist_CaprivaT <- do.call("rbind",lapply(1:length(TempPred), function(x) do.call("rbind", foreHWdist_Capriva[[x]])))
HWdist_CaprivaT$method <- "Perk"

HWdist_Capriva <- rbind(HWdist_CaprivaS, HWdist_CaprivaT)

## Trieste ----
load("ExtractedResults/TempPredTrieste_up2018_Bayes.Rda")
load("ExtractedResults/StatesInTrieste_up2018_Bayes.Rda")

outTraj <- lapply(1:length(TempPred), function(x){
  lapply(1:length(TempPred[[1]]), function(z) TempPred[[x]][[z]][, idx])
})

outStates <- lapply(1:length(TempPred), function(x){
  lapply(1:length(TempPred[[1]]), function(z) StatesIn[[x]][[z]][, idx])
})

foreHWdist_Trieste <- getForeHWdist2(trajectories = outTraj, 
                                     states = outStates,
                                     year_start = 1995, 
                                     year_end = 2018, 
                                     AR = AR, 
                                     model_names = model_names2,
                                     K = label_K, 
                                     last.regime = label_LASTreg2,
                                     site = "Trieste", 
                                     Nsim = 200,
                                     year_thresh = 1995) 
HWdist_TriesteS<- do.call("rbind",lapply(1:length(TempPred), function(x) do.call("rbind", foreHWdist_Trieste[[x]]  )))
HWdist_TriesteS$method <- "State"


foreHWdist_Trieste <- getForeHWdist(trajectories = outTraj, 
                                    dataHWr = extract_HWsw_summary_T$res, 
                                    year_start = 1995, 
                                    year_end = 2018, 
                                    AR = AR, 
                                    model_names = model_names2,
                                    K = label_K, 
                                    last.regime = label_LASTreg2,
                                    site = "Trieste", 
                                    Nsim = 200,
                                    year_thresh = 1995) 
HWdist_TriesteT <- do.call("rbind",lapply(1:length(TempPred), function(x) do.call("rbind", foreHWdist_Trieste[[x]])))
HWdist_TriesteT$method <- "Perk"

HWdist_Trieste <- rbind(HWdist_TriesteS, HWdist_TriesteT)

HWdist_in <- rbind.data.frame(HWdist_Capriva, HWdist_Trieste)

# Plot Main ----

HWdist_out$sample <- "out"
HWdist_in$sample <- "in"

HWdist_all <- rbind.data.frame(HWdist_in, HWdist_out)

load("data/dataFVG/datiCapriva.Rdata") # to extract GST
dataTrain <- dati[dati$Year <= 2018 &  dati$Month >= 5 & dati$Month <= 9,]
dataVal <- dati[dati$Year > 2018 & dati$Year <= 2024 & dati$Month >= 5 & dati$Month <= 9,]

GST <- cbind.data.frame(GST = c(unique(dataTrain$GST_c),unique(dataVal$GST_c)),
                        year = 1995:2024, 
                        GSTclass = cut(c(unique(dataTrain$GST_c)+0.76,unique(dataVal$GST_c)+0.76),
                                       breaks = seq(-.41,.59,.2)+0.76)) #GSTmean 0.76

HWdist_all <- merge(HWdist_all, GST, by = "year")
HWrData <- merge(database_Hwr, GST, by = "year")

HWdist_all$last.regime <- factor(HWdist_all$last.regime, levels = c("N", "G", "rW"))
levels(HWdist_all$last.regime) <- c("Normal", "Gumbel", "rWeibull")

HWrData$sample <- "out"
HWrData$sample[HWrData$year<2019] <- "in" 


HWrData$K <- "HWr"
HWrData$last.regime <- "HWr"
HWrData$method <- "HWr"
colnames(HWdist_all)[3] <- "Type"
HWdist_all_v2 <- rbind.data.frame(HWdist_all[,-c(2,4)], 
                                  HWrData[,c(1,2,10,11,4,3,5,12,9,7,8)])

#### Capriva ----
p1 <- ggplot(HWdist_all_v2[HWdist_all_v2$site == "Capriva"&HWdist_all_v2$measure == "HWd"&HWdist_all_v2$K!="HWr"
                     &HWdist_all_v2$method == "State",], 
       aes(x = as.factor(K), y = Freq, fill = last.regime)) +
  geom_boxplot(col = "gray40", outlier.size = .5, alpha = 0.5) +
  geom_jitter(data = HWdist_all_v2[HWdist_all_v2$site == "Capriva"&HWdist_all_v2$measure == "HWd"&HWdist_all_v2$K=="HWr"
                                   &HWdist_all_v2$method == "HWr",], 
              alpha = 0.5, col = "black", size = 1, width = 0.2, height = 0)+
  facet_grid(cols = vars(GSTclass), rows = vars(sample), 
             scales = "fixed", labeller = "label_value")+
  labs(x = "K", y = "Days in heat wave") +
  theme(legend.position="top",
        legend.text = element_text(size=11),
        text = element_text(family = "serif", size = 11))+ 
  scale_y_continuous(breaks = breaks_pretty())+
  scale_fill_manual(
    name = "Last regime distribution",
    values = c(
      "Normal" = "#666666", "Gumbel" = "#CE1126", "rWeibull" = "#007A3D")
  )+
  scale_x_discrete(guide = guide_axis(n.dodge = 1))

p2 <- ggplot(HWdist_all_v2[HWdist_all_v2$site == "Capriva"&HWdist_all_v2$measure == "HWn"&HWdist_all_v2$K!="HWr"
                     &HWdist_all_v2$method == "State",], 
       aes(x = as.factor(K), y = Freq, fill = last.regime)) +
  geom_boxplot(col = "gray40", outlier.size = .5, alpha = 0.5) +
  geom_jitter(data = HWdist_all_v2[HWdist_all_v2$site == "Capriva"&HWdist_all_v2$measure == "HWn"&HWdist_all_v2$K=="HWr"
                                   &HWdist_all_v2$method == "HWr",], 
              alpha = 0.5, col = "black", size = 1, width = 0.2, height = 0)+
  facet_grid(cols = vars(GSTclass), rows = vars(sample), 
             scales = "fixed", labeller = "label_value")+
  labs(x = "K", y = "Number of heat waves") +
  theme(legend.position="top",
        legend.text = element_text(size=11),
        text = element_text(family = "serif", size = 11))+ 
  scale_y_continuous(breaks = breaks_pretty())+
  scale_fill_manual(
    name = "Last regime distribution",
    values = c(
      "Normal" = "#666666", "Gumbel" = "#CE1126", "rWeibull" = "#007A3D")
  )+
  scale_x_discrete(guide = guide_axis(n.dodge = 1))

ggarrange(p1,p2, common.legend = T, ncol = 1)

ggsave("Plots/HW_C_GST_review_state.pdf",   
       width = 174, height = 190, units = "mm", dpi = 600)

### Trieste ----
p1 <- ggplot(HWdist_all_v2[HWdist_all_v2$site == "Trieste"&HWdist_all_v2$measure == "HWd"&HWdist_all_v2$K!="HWr"
                     &HWdist_all_v2$method == "State",], 
       aes(x = as.factor(K), y = Freq, fill = last.regime)) +
  geom_boxplot(col = "gray40", outlier.size = .5, alpha = 0.5) +
  geom_jitter(data = HWdist_all_v2[HWdist_all_v2$site == "Trieste"&HWdist_all_v2$measure == "HWd"&HWdist_all_v2$K=="HWr"
                                   &HWdist_all_v2$method == "HWr",], 
              alpha = 0.5, col = "black", size = 1, width = 0.2, height = 0)+
  facet_grid(cols = vars(GSTclass), rows = vars(sample), 
             scales = "fixed", labeller = "label_value")+
  labs(x = "K", y = "Days in heat waves") +
  theme(legend.position="top",
        legend.text = element_text(size=11),
        text = element_text(family = "serif", size = 11))+ 
  scale_y_continuous(breaks = breaks_pretty())+
  scale_fill_manual(
    name = "Last regime distribution",
    values = c(
      "Normal" = "#666666", "Gumbel" = "#CE1126", "rWeibull" = "#007A3D")
  )+
  scale_x_discrete(guide = guide_axis(n.dodge = 1))

p2 <- ggplot(HWdist_all_v2[HWdist_all_v2$site == "Trieste"&HWdist_all_v2$measure == "HWn"&HWdist_all_v2$K!="HWr"
                     &HWdist_all_v2$method == "State",], 
       aes(x = as.factor(K), y = Freq, fill = last.regime)) +
  geom_boxplot(col = "gray40", outlier.size = .5, alpha = 0.5) +
  geom_jitter(data = HWdist_all_v2[HWdist_all_v2$site == "Trieste"&HWdist_all_v2$measure == "HWn"&HWdist_all_v2$K=="HWr"
                                   &HWdist_all_v2$method == "HWr",], 
              alpha = 0.5, col = "black", size = 1, width = 0.2, height = 0)+
  facet_grid(cols = vars(GSTclass), rows = vars(sample), 
             scales = "fixed", labeller = "label_value")+
  labs(x = "K", y = "Number of heat waves") +
  theme(legend.position="top",
        legend.text = element_text(size=11),
        text = element_text(family = "serif", size = 11))+ 
  scale_y_continuous(breaks = breaks_pretty())+
  scale_fill_manual(
    name = "Last regime distribution",
    values = c(
      "Normal" = "#666666", "Gumbel" = "#CE1126", "rWeibull" = "#007A3D")
  )+
  scale_x_discrete(guide = guide_axis(n.dodge = 1))

ggarrange(p1,p2, common.legend = T, ncol = 1)

ggsave("Plots/HW_T_GST_review_state.pdf",   
       width = 174, height = 190, units = "mm", dpi = 600)



# Plot Supplementary ----
### Capriva ----
p1 <- ggplot(HWdist_all_v2[HWdist_all_v2$site == "Capriva"&HWdist_all_v2$measure == "HWd"&HWdist_all_v2$K!="HWr"
                           &HWdist_all_v2$method == "Perk",], 
             aes(x = as.factor(K), y = Freq, fill = last.regime)) +
  geom_boxplot(col = "gray40", outlier.size = .5, alpha = 0.5) +
  geom_jitter(data = HWdist_all_v2[HWdist_all_v2$site == "Capriva"&HWdist_all_v2$measure == "HWd"&HWdist_all_v2$K=="HWr"
                                   &HWdist_all_v2$method == "HWr",], 
              alpha = 0.5, col = "black", size = 1, width = 0.2, height = 0)+
  facet_grid(cols = vars(GSTclass), rows = vars(sample), 
             scales = "fixed", labeller = "label_value")+
  labs(x = "K", y = "Days in heat wave") +
  theme(legend.position="top",
        legend.text = element_text(size=11),
        text = element_text(family = "serif", size = 11))+ 
  scale_y_continuous(breaks = breaks_pretty())+
  scale_fill_manual(
    name = "Last regime distribution",
    values = c(     
      "Normal" = "#666666", "Gumbel" = "#CE1126", "rWeibull" = "#007A3D")
  )+
  scale_x_discrete(guide = guide_axis(n.dodge = 1))

p2 <- ggplot(HWdist_all_v2[HWdist_all_v2$site == "Capriva"&HWdist_all_v2$measure == "HWn"&HWdist_all_v2$K!="HWr"
                           &HWdist_all_v2$method == "Perk",], 
             aes(x = as.factor(K), y = Freq, fill = last.regime)) +
  geom_boxplot(col = "gray40", outlier.size = .5, alpha = 0.5) +
  geom_jitter(data = HWdist_all_v2[HWdist_all_v2$site == "Capriva"&HWdist_all_v2$measure == "HWn"&HWdist_all_v2$K=="HWr"
                                   &HWdist_all_v2$method == "HWr",], 
              alpha = 0.5, col = "black", size = 1, width = 0.2, height = 0)+
  facet_grid(cols = vars(GSTclass), rows = vars(sample), 
             scales = "fixed", labeller = "label_value")+
  labs(x = "K", y = "Number of heat waves") +
  theme(legend.position="top",
        legend.text = element_text(size=11),
        text = element_text(family = "serif", size = 11))+ 
  scale_y_continuous(breaks = breaks_pretty())+
  scale_fill_manual(
    name = "Last regime distribution",
    values = c(
      "Normal" = "#666666", "Gumbel" = "#CE1126", "rWeibull" = "#007A3D")
  )+
  scale_x_discrete(guide = guide_axis(n.dodge = 1))

ggarrange(p1,p2, common.legend = T, ncol = 1)

ggsave("Plots/HW_C_GST_review_HWr.pdf",   
       width = 174, height = 190, units = "mm", dpi = 600)

### Trieste ----
p1 <- ggplot(HWdist_all_v2[HWdist_all_v2$site == "Trieste"&HWdist_all_v2$measure == "HWd"&HWdist_all_v2$K!="HWr"
                           &HWdist_all_v2$method == "Perk",], 
             aes(x = as.factor(K), y = Freq, fill = last.regime)) +
  geom_boxplot(col = "gray40", outlier.size = .5, alpha = 0.5) +
  geom_jitter(data = HWdist_all_v2[HWdist_all_v2$site == "Trieste"&HWdist_all_v2$measure == "HWd"&HWdist_all_v2$K=="HWr"
                                   &HWdist_all_v2$method == "HWr",], 
              alpha = 0.5, col = "black", size = 1, width = 0.2, height = 0)+
  facet_grid(cols = vars(GSTclass), rows = vars(sample), 
             scales = "fixed", labeller = "label_value")+
  labs(x = "K", y = "Days in heat waves") +
  theme(legend.position="top",
        legend.text = element_text(size=11),
        text = element_text(family = "serif", size = 11))+ 
  scale_y_continuous(breaks = breaks_pretty())+
  scale_fill_manual(
    name = "Last regime distribution",
    values = c(
      "Normal" = "#666666", "Gumbel" = "#CE1126", "rWeibull" = "#007A3D")
  )+
  scale_x_discrete(guide = guide_axis(n.dodge = 1))

p2 <- ggplot(HWdist_all_v2[HWdist_all_v2$site == "Trieste"&HWdist_all_v2$measure == "HWn"&HWdist_all_v2$K!="HWr"
                           &HWdist_all_v2$method == "Perk",], 
             aes(x = as.factor(K), y = Freq, fill = last.regime)) +
  geom_boxplot(col = "gray40", outlier.size = .5, alpha = 0.5) +
  geom_jitter(data = HWdist_all_v2[HWdist_all_v2$site == "Trieste"&HWdist_all_v2$measure == "HWn"&HWdist_all_v2$K=="HWr"
                                   &HWdist_all_v2$method == "HWr",], 
              alpha = 0.5, col = "black", size = 1, width = 0.2, height = 0)+
  facet_grid(cols = vars(GSTclass), rows = vars(sample), 
             scales = "fixed", labeller = "label_value")+
  labs(x = "K", y = "Number of heat waves") +
  theme(legend.position="top",
        legend.text = element_text(size=11),
        text = element_text(family = "serif", size = 11))+ 
  scale_y_continuous(breaks = breaks_pretty())+
  scale_fill_manual(
    name = "Last regime distribution",
    values = c(
      "Normal" = "#666666", "Gumbel" = "#CE1126", "rWeibull" = "#007A3D")
  )+
  scale_x_discrete(guide = guide_axis(n.dodge = 1))

ggarrange(p1,p2, common.legend = T, ncol = 1)

ggsave("Plots/HW_T_GST_review_HWr.pdf",   
       width = 174, height = 190,  units = "mm", dpi = 600)

# Plot Supplementary by Year----

HWdist_all_v2$last.regime <- factor(HWdist_all_v2$last.regime, 
                                    levels= c("Normal", "Gumbel", "rWeibull", "HWr"))
HWrData$K <- 2
HWrData_v2 <- rbind(HWrData,HWrData, HWrData, HWrData)
HWrData_v2$K <- rep(2:5, each =120)
HWrData_v2$last.regime <- "Normal"
HWrData_v2 <- rbind(HWrData_v2,HWrData_v2, HWrData_v2)
HWrData_v2$last.regime <- rep(c("Normal", "Gumbel", "rWeibull"), each =480)
HWrData_v2$last.regime <- factor(HWrData_v2$last.regime, 
                                    levels= c("Normal", "Gumbel", "rWeibull", "HWr"))
ggplot(HWdist_all_v2[HWdist_all_v2$site == "Capriva"&HWdist_all_v2$sample=="out"&HWdist_all_v2$last.regime!="HWr"&HWdist_all_v2$method=="Perk",], aes(x = as.factor(year), y = Freq)) +
  geom_boxplot(col = "gray40", outlier.size = 1) +
  geom_point(data = HWrData_v2[HWrData_v2$site == "Capriva" & HWrData_v2$year>2018,], aes(x = as.factor(year)), col = "darkred", size = 1)+
  facet_grid(cols = vars(K), rows = vars(last.regime,measure), scales = "free", labeller = "label_value")+
  labs(x = "Year", y = "") +
  theme(legend.position="top",
        legend.text = element_text(size=9),
        text = element_text(family = "serif", size = 9))+ 
  scale_y_continuous(breaks = breaks_pretty())+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))
ggsave("Plots/HW_C_Year.pdf", dpi = 1000, width = 200, height = 230, units = "mm")

ggplot(HWdist_all_v2[HWdist_all_v2$site == "Trieste"&HWdist_all_v2$sample=="out"&HWdist_all_v2$last.regime!="HWr"&HWdist_all_v2$method=="Perk",], aes(x = as.factor(year), y = Freq)) +
  geom_boxplot(col = "gray40", outlier.size = 1) +
  geom_point(data = HWrData_v2[HWrData_v2$site == "Trieste" & HWrData_v2$year>2018,], aes(x = as.factor(year)), col = "darkred", size = 1)+
  facet_grid(cols = vars(K), rows = vars(last.regime,measure), scales = "free", labeller = "label_value")+
  labs(x = "Year", y = "") +
  theme(legend.position="top",
        legend.text = element_text(size=9),
        text = element_text(family = "serif", size = 9))+ 
  scale_y_continuous(breaks = breaks_pretty())+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))
ggsave("Plots/HW_T_Year.pdf", dpi = 1000, width = 200, height = 230, units = "mm")
