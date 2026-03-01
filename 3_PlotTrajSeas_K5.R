# This code reproduce figure 3 in the main manuscript 

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

location <- "Capriva"
stations <- "_C"

comm_seas <- TRUE
comm_trend <- FALSE
AR <- FALSE
comm_AR <- FALSE

load("dati/dataFVG/datiCapriva.Rdata")
dataTrain <- dati[dati$Year <= year_end &  dati$Month >= 5 & dati$Month <= 9,]
dataTest <- dati[dati$Year > year_end &  dati$Month >= 5 & dati$Month <= 9,]


model_names <- paste0("fit_", as.vector(t(outer(label_reg, label_LASTreg, FUN = paste0))), stations, "_s")
model_names2 <- paste0("fit_", as.vector(t(outer(label_reg, label_LASTreg, FUN = paste0))), stations)
file_select <- paste0(model_names,  ".RData") 
fileName <- list.files(path = "ResultsGST", pattern = stations)
get_model <- which(file_select%in%fileName)

i <- 10  
K <- label_K[i]

load(paste0("ResultsGST/",file_select[i]))
fit <- mget(model_names2[i], envir = globalenv())[[1]] 
rm(list=grep(model_names2[i],ls(),value=TRUE))

Posterior <- ExtractPost(stan.obj = fit, K = label_K[i], last.regime = label_LASTreg2[i], 
                         comm_seas = comm_seas, comm_trend = comm_trend)
QuantilesNEW <- extract_Quantiles(Posterior, data = dataTest, AR = AR, K = 5, 
                               last.regime = "N", quantiles = c(0.025,0.25,0.5,0.75,0.975))

Quantiles5 <- sapply(1:5, function(x) QuantilesNEW[[x]][,5])
Quantiles5 <- Quantiles5[766:918,]
colnames(Quantiles5) <- paste0("qNEW5",c("0.025","0.25","0.5","0.75", "0.975"))

predictedValues <- extract_PredSummary(Posterior, dataTest, AR = AR,
                                       K = label_K[i], last.regime = label_LASTreg2[i])

idx2024 <-   ((2024 - 2019) * 153 + 1) : ((2024 - 2019 + 1) * 153)   

Quantiles1 <- predictedValues$mean_reg_quantiles[,idx2024,]
Quantiles <- t(apply(Quantiles1, 2, c))


colnames(Quantiles) <- paste0("q",rep(1:K,each=5), c("0.025","0.25","0.5","0.75", "0.975"))
dati <- dataTest[idx2024,c("doy", "MaxT")] 

plotSeasDataC <- cbind(dati, Quantiles, Quantiles5,
                       site = location) 


## Trieste----

location <- "Trieste"
stations <- "_T"

comm_seas <- TRUE
comm_trend <- FALSE
AR <- FALSE
comm_AR <- FALSE

load("dati/dataFVG/datiTrieste.Rdata")
dataTrain <- dati[dati$Year <= year_end &  dati$Month >= 5 & dati$Month <= 9,]


model_names <- paste0("fit_", as.vector(t(outer(label_reg, label_LASTreg, FUN = paste0))), stations, "_s")
model_names2 <- paste0("fit_", as.vector(t(outer(label_reg, label_LASTreg, FUN = paste0))), stations)
file_select <- paste0(model_names,  ".RData")
fileName <- list.files(path = "ResultsGST", pattern = stations)
get_model <- which(file_select%in%fileName)

i <- 10 
K <- label_K[i]

load(paste0("ResultsGST/",file_select[i]))
fit <- mget(model_names2[i], envir = globalenv())[[1]] 
rm(list=grep(model_names2[i],ls(),value=TRUE))


Posterior <- ExtractPost(stan.obj = fit, K = label_K[i], last.regime = label_LASTreg2[i], 
                         comm_seas = comm_seas, comm_trend = comm_trend)
QuantilesNEW <- extract_Quantiles(Posterior, data = dataTest, AR = AR, K = 5, 
                                  last.regime = "N", quantiles = c(0.025,0.25,0.5,0.75,0.975))
Quantiles5 <- sapply(1:5, function(x) QuantilesNEW[[x]][,5])
Quantiles5 <- Quantiles5[766:918,]
colnames(Quantiles5) <- paste0("qNEW5",c("0.025","0.25","0.5","0.75", "0.975"))

predictedValues <- extract_PredSummary(Posterior, dataTest, AR = AR,
                                       K = label_K[i], last.regime = label_LASTreg2[i])

idx2024 <-   ((2024 - 2019) * 153 + 1) : ((2024 - 2019 + 1) * 153)   


Quantiles1 <- predictedValues$mean_reg_quantiles[,idx2024,]
Quantiles <- t(apply(Quantiles1, 2, c))

colnames(Quantiles) <- paste0("q",rep(1:K,each=5), c("0.025","0.25","0.5","0.75", "0.975"))
dati <- dataTrain[idx2024,c("doy", "MaxT")] 

plotSeasDataT <- cbind(dati, Quantiles, Quantiles5,#Prob, 
                       site = location) 

dataSeas <- rbind(plotSeasDataC, plotSeasDataT)


seq_year <- cumsum(c(0,152))+1
label_year <- 1995:2024

idx_origin <- (120:272 - 120)/152
idx <- sample(1:4000,200, replace = FALSE)

load("ExtractedResults/TempTrajCapriva_from2019_to2024_Bayes.Rda")
TempTrajC <- lapply(1:length(TempTraj), function(x){
  lapply(1:length(TempTraj[[1]]), function(z) TempTraj[[x]][[z]][, idx])
})
load("ExtractedResults/TempTrajTrieste_from2019_to2024_Bayes.Rda")
TempTrajT <- lapply(1:length(TempTraj), function(x){
  lapply(1:length(TempTraj[[1]]), function(z) TempTraj[[x]][[z]][, idx])
})


TempTraj <- rbind.data.frame(cbind.data.frame(TempTrajC[[10]][[6]], site = "Capriva"),
                             cbind.data.frame(TempTrajT[[10]][[6]], site = "Trieste"))
traj <- reshape2::melt(TempTraj)
traj$doy <- 1:153


load("HeatwaveR/extract_HWsw_summary_C.Rdata")
load("HeatwaveR/extract_HWsw_summary_Cwith2024.Rdata")

load("HeatwaveR/extract_HWsw_summary_T.Rdata")
load("HeatwaveR/extract_HWsw_summary_Twith2024.Rdata")




HWrTrhC <- extract_HWsw_summary_C3$res[extract_HWsw_summary_C3$res$years==2024,]
HWrTrhT <- extract_HWsw_summary_T3$res[extract_HWsw_summary_T3$res$years==2024,]
HWrTrh <- rbind.data.frame(cbind.data.frame(HWrTrhC, site = "Capriva"),
                           cbind.data.frame(HWrTrhT, site = "Trieste"))

## plot paper -----
ggplot(traj, aes(x = doy))+
  geom_line(aes(group = variable,y = value),alpha = 0.2, linewidth = 0.2, col = "grey")+
  geom_line(data = traj[traj$variable=="147",], aes(y = value, group = variable), linewidth = 0.2, col = "grey20", linetype = 2)+
  geom_line(data = HWrTrh, aes(x = doy, y = MaxT), col = "grey20", linetype = 1, linewidth = 0.2)+
  geom_line(data = dataSeas, aes(y = q10.5), linewidth = .3, color = "darkgreen")+
  geom_line(data = dataSeas, aes(y = q20.5), linewidth = .3, color = "gold2")+
  geom_line(data = dataSeas, aes(y = q30.5), linewidth = .3, color = "darkorange")+
  geom_line(data = dataSeas, aes(y = q40.5), linewidth = .3, color = "red")+
  geom_line(data = dataSeas, aes(y = q50.5), linewidth = .3, color = "darkred")+
  facet_wrap(~site,nrow = 2,  scales = "free")+
  xlab("Year 2024") + 
  theme_default()+
  labs(y = "Maximum daily temperature (Celsius)") +
  theme(text = element_text(family = "serif", size = 9),
        legend.key.width  = unit(.3, "lines"),
        legend.key.height = unit(4, "lines"),
        legend.position="right",
        legend.text=element_text(size=9))+
  scale_x_continuous(breaks=label_month,
                     labels=c("May", "June", "July", "August", "September"))+ 
  scale_y_continuous(breaks = breaks_pretty())

ggsave("Plots/Traj2024_NNNNN.pdf", dpi = 1000, width = 190, height = 160, units = "mm")

ggplot(traj, aes(x = doy))+
  geom_line(data = traj, aes(group = variable,y = value),alpha = 0.3, linewidth = 0.2, col = "grey")+
  geom_ribbon(data = dataSeas, aes(ymin = qNEW50.025, ymax = qNEW50.975), fill = "darkred", alpha = 0.3) +
  geom_ribbon(data = dataSeas, aes(ymin = qNEW50.25, ymax = qNEW50.75), fill = "darkred", alpha = 0.2) +
 geom_line(data = dataSeas, aes(y = q10.5), linewidth = .3, color = "darkgreen")+
  geom_line(data = dataSeas, aes(y = q20.5), linewidth = .3, color = "gold2")+
  geom_line(data = dataSeas, aes(y = q30.5), linewidth = .3, color = "darkorange")+
  geom_line(data = dataSeas, aes(y = q40.5), linewidth = .3, color = "red")+
  geom_line(data = dataSeas, aes(y = qNEW50.5), linewidth = .3, color = "darkred")+
  facet_wrap(~site,nrow = 2,  scales = "free")+
  xlab("Year 2024") + 
  theme_default()+
  labs(y = "Maximum daily temperature (Celsius)") +
  theme(text = element_text(family = "serif", size = 9),
        legend.key.width  = unit(.3, "lines"),
        legend.key.height = unit(4, "lines"),
        legend.position="right",
        legend.text=element_text(size=9))+
  scale_x_continuous(breaks=label_month,
                     labels=c("May", "June", "July", "August", "September"))+ 
  scale_y_continuous(breaks = breaks_pretty())
ggsave("Plots/Traj2024_NNNNN_SM.pdf", dpi = 1000, width = 190, height = 160, units = "mm")
