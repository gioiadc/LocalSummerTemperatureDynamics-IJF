# This code reproduce figure 3 in the main manuscript 
# and Figure S7 in the Supplementary Material 

source("functions/ExtractPosterior.R")
source("functions/load_Packages.R")

seq_month <- rep(0:153)
label_month <- seq_month[c(1, 32, 62, 93, 124)] 

idx2024 <- ((2024 - 2019) * 153 + 1) : ((2024 - 2019 + 1) * 153)   
K <- 5

year_end <- 2018

#set.seed(23)
idx <- sample(1:4000,200, replace = FALSE)

# Capriva ----
location <- "Capriva"
load("data/dataFVG/datiCapriva.Rdata")
dataTrain <- dati[dati$Year <= year_end &  dati$Month >= 5 & dati$Month <= 9,]
dataTest <- dati[dati$Year > year_end &  dati$Month >= 5 & dati$Month <= 9,]


load("ExtractedResults/QuantilesCapriva_from2019_to2024_Bayes.Rda")
load("ExtractedResults/Quantiles5Capriva_from2019_to2024_Bayes.Rda")
dati <- dataTest[idx2024,c("doy", "MaxT")] 

plotSeasDataC <- cbind(dati, 
                       Quantiles, 
                       Quantiles5,
                       site = location) 

# Trieste ----
location <- "Trieste"
load("data/dataFVG/datiTrieste.Rdata")
dataTrain <- dati[dati$Year <= year_end &  dati$Month >= 5 & dati$Month <= 9,]
dataTest <- dati[dati$Year > year_end &  dati$Month >= 5 & dati$Month <= 9,]


load("ExtractedResults/QuantilesTrieste_from2019_to2024_Bayes.Rda")
load("ExtractedResults/Quantiles5Trieste_from2019_to2024_Bayes.Rda")

dati <- dataTest[idx2024,c("doy", "MaxT")] 

plotSeasDataT <- cbind(dati, 
                       Quantiles, 
                       Quantiles5,
                       site = location) 


dataSeas <- rbind(plotSeasDataC, plotSeasDataT)

seq_year <- cumsum(c(0,152))+1
label_year <- 1995:2024

idx_origin <- (120:272 - 120)/152

load("ExtractedResults/TempTrajCapriva_from2019_to2024_Bayes_1.Rda")
load("ExtractedResults/TempTrajCapriva_from2019_to2024_Bayes_2.Rda")
TempTraj <- c(TempTraj1, TempTraj2)
rm(TempTraj1, TempTraj2)

TempTrajC <- lapply(1:length(TempTraj), function(x){
  lapply(1:length(TempTraj[[1]]), function(z) TempTraj[[x]][[z]][, idx])
})

load("ExtractedResults/TempTrajTrieste_from2019_to2024_Bayes_1.Rda")
load("ExtractedResults/TempTrajTrieste_from2019_to2024_Bayes_2.Rda")
TempTraj <- c(TempTraj1, TempTraj2)
rm(TempTraj1, TempTraj2)

TempTrajT <- lapply(1:length(TempTraj), function(x){
   lapply(1:length(TempTraj[[1]]), function(z) TempTraj[[x]][[z]][, idx])
 })

TempTraj <- rbind.data.frame(cbind.data.frame(TempTrajC[[4]][[6]], site = "Capriva"),
                             cbind.data.frame(TempTrajT[[4]][[6]], site = "Trieste"))
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
  geom_line(data = dataSeas, aes(y = q10.5), linewidth = .3, color = "#032163")+
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
  geom_line(data = dataSeas, aes(y = q10.5), linewidth = .3, color = "#032163")+
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

