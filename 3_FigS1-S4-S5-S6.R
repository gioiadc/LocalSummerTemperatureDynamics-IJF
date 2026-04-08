# This file reproduce results in the Supplementary Material:
# - the maximum temperature plot for 2023 and 2024
# - the extremal dependence measures (chi and theta)

source("functions/load_Packages.R")

# Maximum temperature plot for 2023 and 2024 years with the empirical based threshold ----

load("data/dataFVG/datiCapriva.Rdata")
datiC <- dati[dati$Year>2022,c("MaxT", "Year", "doy")]
datiC$Site <- "Capriva"
load("data/dataFVG/datiTrieste.Rdata")
datiT <- dati[dati$Year>2022,c("MaxT", "Year", "doy")]
datiT$Site <- "Trieste"

load("HeatwaveR/extract_HWsw_summary_Cwith2024.Rdata")
load("HeatwaveR/extract_HWsw_summary_Twith2024.Rdata")
datiC$thresh <- extract_HWsw_summary_C3$res$tresh[1:306]
datiT$thresh <- extract_HWsw_summary_T3$res$tresh[1:306]

datiPlot <- rbind(datiC, datiT)

seq_month <- rep(0:153)
label_month <- seq_month[c(1, 32, 62, 93, 124)]

ggplot(data = datiPlot, aes(x = doy, group = Site)) +
  geom_line(data = datiPlot, aes(y = thresh), col= darker_red, lty = 3, linewidth = .3) +
  geom_point(data = datiPlot, aes(y = MaxT), col = "grey20", size = 0.5, alpha = 0.5)+
  geom_line(data = datiPlot, aes(y = MaxT), alpha = 0.3, col = "gray60")+
  xlab("Time") + 
  theme_default()+
  labs(y = "Maximum daily temperature (Celsius)") +
  theme(legend.position="none",
        #legend.text = element_text(size=9),
        text = element_text(family = "serif", size = 10))+
  scale_x_continuous(breaks=label_month,
                     labels=c("May", "June", "July", "Aug", "Sept"))+
  facet_grid(cols = vars(Year), rows = vars(Site), scales = "free")

ggsave(file="Plots/FigS1.pdf", width = 214, height = 120, 
       dpi = 600, units = "mm")



# Extremal dependence index for Capriva, last regime Normal ----

lags <- 1 # Set lags1 to obtain Fig S4, then set 5 to obtain Fig S5
u <- seq(0.5,0.98,0.01)

# empirical result

load("data/dataFVG/datiCapriva.Rdata")
x <- dati$MaxT[dati$Year>2018]

X <- cbind(head(x,length(x)-lags), x[-c(1:lags)])

chi_vals <- sapply(lags, function(h) chi(X, nq = length(u), qlim = range(u)))
chi_vals <- chi_vals[[1]]

# MS results

load("ExtractedResults/TempTrajCapriva_from2019_to2024_Bayes_1.Rda")
load("ExtractedResults/TempTrajCapriva_from2019_to2024_Bayes_2.Rda")

TempTraj <- c(TempTraj1, TempTraj2)
rm(TempTraj1, TempTraj2)

X <- do.call("rbind",TempTraj[[1]])
matK2N <- list()
for(i in 1:4000){
  Xs <- cbind(head(X[,i],length(X[,i])-lags), X[-c(1:lags),i])
  temp <- try(sapply(lags, function(h) chi(Xs, nq = length(u), qlim = range(u))), TRUE)
  if(!is.character(temp)){
    matK2N[[i]] <- temp[[1]]
  } else {
    matK2N[[i]] <- matrix(NA, nrow = length(u), ncol = 3)
  }
  print(i)
}
matK2N <- sapply(matK2N, function(mat) mat[, 2])

X <- do.call("rbind",TempTraj[[2]])
matK3N <- list()
for(i in 1:4000){
  Xs <- cbind(head(X[,i],length(X[,i])-lags), X[-c(1:lags),i])
  temp <- try(sapply(lags, function(h) chi(Xs, nq = length(u), qlim = range(u))), TRUE)
  if(!is.character(temp)){
    matK3N[[i]] <- temp[[1]]
  } else {
    matK3N[[i]] <- matrix(NA, nrow = length(u), ncol = 3)
  }
  print(i)
}
matK3N <- sapply(matK3N, function(mat) mat[, 2])

X <- do.call("rbind",TempTraj[[3]])
matK4N <- list()
for(i in 1:4000){
  Xs <- cbind(head(X[,i],length(X[,i])-lags), X[-c(1:lags),i])
  temp <- try(sapply(lags, function(h) chi(Xs, nq = length(u), qlim = range(u))), TRUE)
  if(!is.character(temp)){
    matK4N[[i]] <- temp[[1]]
  } else {
    matK4N[[i]] <- matrix(NA, nrow = length(u), ncol = 3)
  }
  print(i)
}
matK4N <- sapply(matK4N, function(mat) mat[, 2])

X <- do.call("rbind",TempTraj[[4]])
matK5N <- list()
for(i in 1:4000){
  Xs <- cbind(head(X[,i],length(X[,i])-lags), X[-c(1:lags),i])
  temp <- try(sapply(lags, function(h) chi(Xs, nq = length(u), qlim = range(u))), TRUE)
  if(!is.character(temp)){
    matK5N[[i]] <- temp[[1]]
  } else {
    matK5N[[i]] <- matrix(NA, nrow = length(u), ncol = 3)
  }
  print(i)
}
matK5N <- sapply(matK5N, function(mat) mat[, 2])


K2 <- melt(matK2N)
K2$K <- 2
K3 <- melt(matK3N)
K3$K <- 3
K4 <- melt(matK4N)
K4$K <- 4
K5 <- melt(matK5N)
K5$K <- 5

matPlot <- rbind.data.frame(K2, K3, K4, K5)
matPlot$u <- u

if(lags == 1){
  chi_valsPlot <- cbind.data.frame(u = u, value = chi_vals[1:length(u),2],
                                   valueLB = chi_vals[1:length(u),1],
                                   valueUB = chi_vals[1:length(u),3],
                                   K = rep(2:5, each = length(u)))
} else {
  chi_valsPlot <- cbind.data.frame(u = u, value = chi_vals[1:length(u),2],
                                   valueLB = pmax(0,chi_vals[1:length(u),1]),
                                   valueUB = chi_vals[1:length(u),3],
                                   K = rep(2:5, each = length(u)))
}


ggplot(matPlot)+
  geom_line(aes(x = u, y = value, group = X2), col = "darkgrey", alpha = 0.2)+ 
  geom_line(data = chi_valsPlot, aes(x = u, y = value), col = "steelblue")+
  geom_ribbon(data = chi_valsPlot, aes(x = u, ymin = valueLB, ymax = valueUB), 
              fill = "steelblue", col = "steelblue", alpha = 0.5)+
  geom_line(data = chi_valsPlot, aes(x = u, y = value))+
  ylim(c(-0.15,1))+
  ylab(expression(paste(chi,"(u)")))+
  theme_light()+
  facet_wrap(vars(K), ncol = 1)
ggsave(paste0("Plots/chi_N_lag",lags,"_OUT_98.pdf"), width = 15, height = 15, units = "cm")

## Index theta for Capriva, last regime Normal ----

# empirical result
load("data/dataFVG/datiCapriva.Rdata")
x <- dati$MaxT[dati$Year>2018]

load("HeatwaveR/extract_HWsw_summary_Cwith2024.Rdata")
Perk <- extract_HWsw_summary_C3$res$tresh
thetaPerk <- extremalindex(x, Perk, method = c("runs"), run.length = 1)[1]
thetaPerk <-  cbind(theta = thetaPerk, K = 2:5)

# MS results

X_list <- list(
  do.call("rbind", TempTraj[[1]]),
  do.call("rbind", TempTraj[[2]]),
  do.call("rbind", TempTraj[[3]]),
  do.call("rbind", TempTraj[[4]])
)

thetaK2 <- numeric(4000)
thetaK3 <- numeric(4000)
thetaK4 <- numeric(4000)
thetaK5 <- numeric(4000)

for (i in 1:4000) {
  thetaK2[i] <- extremalindex(X_list[[1]][, i], Perk, method = "runs", run.length = 1)[1]
  thetaK3[i] <- extremalindex(X_list[[2]][, i], Perk, method = "runs", run.length = 1)[1]
  thetaK4[i] <- extremalindex(X_list[[3]][, i], Perk, method = "runs", run.length = 1)[1]
  thetaK5[i] <- extremalindex(X_list[[4]][, i], Perk, method = "runs", run.length = 1)[1]
  }

theta <- melt(cbind(thetaK2, thetaK3, thetaK4, thetaK5))
levels(theta$X2) <- c("2", "3", "4", "5")
colnames(theta)[2:3] <- c("K","theta")

ggplot(theta) + 
  geom_histogram(aes(x = theta)) + 
  geom_vline(data = thetaPerk, aes(xintercept = theta))+
  ylab("Freq")+
  xlab(expression(paste(theta)))+
  theme_light()+
  facet_wrap(vars(K), ncol = 1)
ggsave(paste0("Plots/theta_OUT_Perk.pdf"), width = 15, height = 15, units = "cm")

