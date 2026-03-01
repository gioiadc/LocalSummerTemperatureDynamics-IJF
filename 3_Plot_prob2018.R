
source("functions/ExtractPosterior.R")
source("functions/load_Packages.R")
source("functions/fitSTAN.R")

minK <- 2
maxK <- 5
year_start <- 1995
year_end <- 2018
year_end2 <- 2024

seq_month <- rep(0:153)
label_month <- seq_month[c(1, 32, 62, 93, 124)] 
label_LASTreg <- c("N")
label_LASTreg2 <- rep(label_LASTreg, length(minK:maxK))
label_reg <- c("N", "NN", "NNN", "NNNN")
label_K <- rep(minK:maxK, each = length(label_LASTreg))

aaa <- eee <- list()
dataSeas <- list()
HWsegm_all <- list()
yearSel <- c(2018)


comm_seas <- TRUE
comm_trend <- FALSE
AR <- FALSE
comm_AR <- FALSE
locationVect <- c("Capriva", "Trieste")
stationsVect <- c("_C", "_T")

yearSel <- 2018
percVect <- c("", "95")
for(site in 1:2){
  location <- locationVect[site]
  stations <- stationsVect[site]
  ProbFore <- list()
  for(K in 2:5){
  
    Kpos <- K-1
 
    load(paste0("dati/dataFVG/dati",location,".Rdata"))
    dataTrain <- dati[dati$Year <= year_end &  dati$Month >= 5 & dati$Month <= 9,]
    
    model_names <- paste0("fit_", as.vector(t(outer(label_reg, label_LASTreg, FUN = paste0))), stations,"_s")
    model_names2 <- paste0("fit_", as.vector(t(outer(label_reg, label_LASTreg, FUN = paste0))), stations)
    file_select <- paste0(model_names,  ".RData") 
    fileName <- list.files(path = "ResultsGST", pattern = stations)
    get_model <- which(file_select%in%fileName)
    
    load(paste0("ResultsGST/",file_select[Kpos]))
    fit <- mget(model_names2[Kpos], envir = globalenv())[[1]] 
    rm(list=grep(model_names2[Kpos],ls(),value=TRUE))
    
    Posterior <- ExtractPost(stan.obj = fit, K = label_K[Kpos], last.regime = label_LASTreg2[Kpos], 
                             comm_seas = comm_seas, comm_trend = comm_trend)
    
    predictedValues <- extract_PredSummary(Posterior, dataTrain, AR = AR, 
                                           K = label_K[Kpos], last.regime = label_LASTreg2[Kpos])
    
    idxYearSel <-   ((yearSel - 1995) * 153 + 1) : ((yearSel - 1995 + 1) * 153)  
    hat_prob_in <- as.data.frame(apply(exp(Posterior$pi),c(2,3), mean)[idxYearSel,])

    colnames(hat_prob_in) <- 1:K
    hat_prob_in$doy <- 1:153
    
    Quantiles1 <- predictedValues$mean_reg_quantiles[,idxYearSel,]
    Quantiles <- t(apply(Quantiles1, 2, c))
    colnames(Quantiles) <- paste0("q",rep(1:K,each=5), c("0.025","0.25","0.5","0.75", "0.975"))
    
    dati <- dataTrain[idxYearSel, c("doy", "MaxT")] 
    Prob <- hat_prob_in[, K]
    
    plotSeasData <- cbind(dati, Quantiles, Prob, site = location) 
    dataSeas[[site]] <- plotSeasData
    
    # HW detect 
    HWsegmTemp <- list()
    for(perc in c(1:2)){
      load(paste0("HeatwaveR/extract_HWsw_summary",stations,percVect[perc],".Rdata"))
      fileTemp <- mget(paste0("extract_HWsw_summary",stations), envir = globalenv())[[1]] 
      
      Exc <- get_excedance(K = K, Posterior$pi, year_start, year_end, last.regime = label_LASTreg2[Kpos])
      Concordance <- list()
      Concordance$Overall <- cbind(fileTemp$res[,c("years","doy", "MaxT", "tresh", "excedance")], Exc$HWs)
      
      HWs_summary_hwr <- numberingHW(fileTemp$res[fileTemp$res$year == yearSel & 
                                                    fileTemp$res$flag_excedance >= 3,])
      HWs_summary_fit <- numberingHW(Concordance$Overall[Concordance$Overall$year == yearSel & 
                                                           Concordance$Overall[,6] == K & 
                                                           Concordance$Overall$exc >2,])
      if(nrow(HWs_summary_hwr)!=0){
        HWs_summary_hwr$site  <- location
        HWs_summary_hwr$method <- "HWr"
      } else {
        HWs_summary_hwr$site <- numeric(0)
        HWs_summary_hwr$method <- numeric(0)
      }
      
      if(nrow(HWs_summary_fit)!=0){
        HWs_summary_fit$site <- location
        HWs_summary_fit$method <- "Prob"
      } else{
        HWs_summary_fit$site <- numeric(0)
        HWs_summary_fit$method <- numeric(0)
      }  
      
      if(nrow(HWs_summary_fit)!=0 & nrow(HWs_summary_hwr)!=0){
        HWs_summary <- rbind(HWs_summary_fit[,c(2,3,9:11)], HWs_summary_hwr[,c(1,4,9:11)])
      } else if(nrow(HWs_summary_fit)!=0){
        HWs_summary <- rbind(HWs_summary_fit[,c(2,3,9:11)])
      } else if(nrow(HWs_summary_hwr)!=0){
        HWs_summary <- rbind(HWs_summary_hwr[,c(1,4,9:11)])
      }
      
      HWs_summary$id <- paste0(HWs_summary$site,
                                    HWs_summary$method,
                                    HWs_summary$progNumbHW)
      
      HWsegm <- cbind.data.frame(from = tapply(HWs_summary$doy, HWs_summary$id, min),
                                 to = tapply(HWs_summary$doy, HWs_summary$id, max))
      
      HWsegm$site <- factor(substring(rownames(HWsegm),1,1), labels = location)
      
      if(nrow(HWs_summary_fit)!=0 & nrow(HWs_summary_hwr)!=0){
        HWsegm$method <- factor(substring(rownames(HWsegm),nchar(rownames(HWsegm))-3, nchar(rownames(HWsegm))-1), labels = c("HWr","Prob"))
      } else if(nrow(HWs_summary_fit)!=0){
        HWsegm$method <- factor(substring(rownames(HWsegm),nchar(rownames(HWsegm))-3, nchar(rownames(HWsegm))-1), labels = c("Prob"))
      } else if(nrow(HWs_summary_hwr)!=0){
        HWsegm$method <- factor(substring(rownames(HWsegm),nchar(rownames(HWsegm))-3, nchar(rownames(HWsegm))-1), labels = c("HWr"))
      }
      HWsegm$nHW <- substring(rownames(HWsegm),nchar(rownames(HWsegm)), nchar(rownames(HWsegm)))
      
      
      HWsegm$ystart <- 13.1
      
      HWsegm$ystart[HWsegm$method == "Prob"] <- HWsegm$ystart[HWsegm$method == "Prob"] -13.15
      HWsegm$K <- K
      HWsegm$perc <- perc
        
      HWsegmTemp[[perc]]<- HWsegm
    }
    
    HWsegm_all[[Kpos]] <- do.call("rbind", HWsegmTemp)

    ProbFore[[Kpos]] <- melt(hat_prob_in, id.vars = "doy")
    ProbFore[[Kpos]]$K <- K
  }
  
  aaa[[site]] <- do.call("rbind", ProbFore)
  aaa[[site]]$K <- factor(aaa[[site]]$K)
  aaa[[site]]$year <- yearSel
  aaa[[site]]$site <- location
  
  eee[[site]] <- do.call("rbind", HWsegm_all)
  eee[[site]]$year <- yearSel
  eee[[site]]$site <- location
  
  print(site)
}
  
bbb <-  do.call("rbind", aaa)
bbb$col <- as.character(bbb$variable)
bbb$col[bbb$K==2&bbb$col==2] <- 5
bbb$col[bbb$K==3&bbb$col==3] <- 5
bbb$col[bbb$K==4&bbb$col==4] <- 5

ccc <- do.call("rbind", dataSeas)
ccc$year <- rep(2018,each = 153)

ddd <- do.call("rbind", eee)

plotProb <- ggplot(bbb)+
  geom_area(aes(x = doy, y = value, group = variable, fill = col), alpha = 0.7, size=.001, colour="black")+
  scale_fill_manual(values = c("darkgreen", "gold2", "darkorange2","orangered3","darkred"), guide = "none")+  
  theme_default()+
  facet_grid(rows = vars(K), cols = vars(site))+
  theme(text = element_text(family = "serif", size = 11),
        legend.key.width  = unit(1, "lines"),
        legend.key.height = unit(.2, "lines"),
        legend.text=element_text(size=11),
        strip.text.x = element_blank(),
        axis.title.y = element_text(family = "sans",size = 11),
        legend.position="none")+
  coord_cartesian(xlim = c(0, 153))+
  scale_x_continuous(breaks=label_month, limits = c(0,153),
                     labels=c("May", "June", "July", "Aug", "Sept"))+
  scale_y_continuous(breaks=c(0,0.5,1), limits = c(-0.1,1.01))+#,
  xlab("2018")+
  ylab(expression(pi[k]^t))+
  geom_segment(data = droplevels(ddd[ddd$method == "Prob",]), 
               aes(x = from, xend =  to, y = ystart, yend = ystart, linetype =  method), 
               col= "darkred", lwd = 0.4) +
  geom_segment(data = ddd[ddd$method == "Prob",], 
               aes(x = to-0.1, xend =  to, y = ystart, yend = ystart), 
               arrow = arrow(angle = 90, length = unit(0.02, "inches"), ends = "last"), 
               col= "darkred", lty = 1, lwd = 0.4)+
  geom_segment(data = ddd[ddd$method == "Prob",], 
               aes(x = from, xend =  from+0.1, y = ystart, yend = ystart),
               arrow = arrow(angle = 90, length = unit(0.02, "inches"), ends = "first"),  
               col= "darkred", lty = 1, lwd = 0.4) +
  scale_linetype_manual(labels = c(""), 
                        values = c("Prob" = 2), guide = "legend", name ="Max Prob") 

load(paste0("HeatwaveR/extract_HWsw_summary_C.Rdata"))
load(paste0("HeatwaveR/extract_HWsw_summary_T.Rdata"))
Thre <- rbind(cbind(extract_HWsw_summary_C$res[extract_HWsw_summary_C$res$years == 1995,], site="Capriva"),
              cbind(extract_HWsw_summary_T$res[extract_HWsw_summary_T$res$years == 1995,], site = "Trieste"))
load(paste0("HeatwaveR/extract_HWsw_summary_C95.Rdata"))
load(paste0("HeatwaveR/extract_HWsw_summary_T95.Rdata"))
Thre95 <- data.frame(doy = rep(1:153,2),
                         thresh = c(extract_HWsw_summary_C$res[extract_HWsw_summary_C$res$years == 1995, ]$tresh,
                                    extract_HWsw_summary_T$res[extract_HWsw_summary_T$res$years == 1995, ]$tresh),
                         site = rep(c("Capriva", "Trieste"), each = 153))
Thre95 <- cbind(Thre95, Thre[,-c(1,6,9)])

ddd$perc <- as.factor(ddd$perc)

plotThr <- ggplot(data = ccc, aes(x = doy))+
  geom_point(aes(y = MaxT), col = "gray40", fill = "gray40", size = 0.3, alpha = 0.9)+
  geom_line(aes(y = MaxT), alpha = 0.3, col = "gray60")+  
  geom_line(aes(y = MaxT), alpha = 0.3, col = "gray60")+  
  geom_line(data = Thre, aes(y = tresh), alpha = 0.7, col = "darkred", linetype = 3)+  
  geom_line(data = Thre95, aes(y = thresh), alpha = 0.7, col = "darkred", linetype = 1)+
  theme_default()+
  ylim(15,40)+
  theme(text = element_text(family = "serif", size = 11),
        legend.key.width  = unit(1, "lines"),
        legend.key.height = unit(.2, "lines"),
        legend.text=element_text(family = "serif",size=11),
        legend.position="none",
        strip.text.y = element_blank()
  )+
  scale_x_continuous(breaks=label_month,limits = c(0,153),
                     labels=c("May", "June", "July", "Aug", "Sept"))+
  xlab("")+
  coord_cartesian(xlim = c(0, 153))+
  facet_grid(~site)+
  labs(y = "Max daily temperature")+
  geom_segment(data = ddd[ddd$method == "HWr"&ddd$perc==1,], 
               aes(x = from, xend =  to, y = ystart+2, yend = ystart+2, linetype =  perc), 
               col= "darkred", lwd = 0.4) +
  geom_segment(data = ddd[ddd$method == "HWr"&ddd$perc==1,], 
               aes(x = to-0.1, xend =  to, y = ystart+2, yend = ystart+2), 
               arrow = arrow(angle = 90, length = unit(0.02, "inches"), ends = "last"), 
               col= "darkred", lty = 1, lwd = 0.4)+
  geom_segment(data = ddd[ddd$method == "HWr"&ddd$perc==1,], 
               aes(x = from, xend =  from+0.1, y = ystart+2, yend = ystart+2),
               arrow = arrow(angle = 90, length = unit(0.02, "inches"), ends = "first"),  
               col= "darkred", lty = 1, lwd = 0.4) +
  geom_segment(data = ddd[ddd$method == "HWr"&ddd$perc==2,], 
               aes(x = from, xend =  to, y = ystart+3, yend = ystart+3, linetype =  perc), 
               col= "darkred", lwd = 0.4) +
  geom_segment(data = ddd[ddd$method == "HWr"&ddd$perc==2,], 
               aes(x = to-0.1, xend =  to, y = ystart+3, yend = ystart+3), 
               arrow = arrow(angle = 90, length = unit(0.02, "inches"), ends = "last"), 
               col= "darkred", lty = 1, lwd = 0.4)+
  geom_segment(data = ddd[ddd$method == "HWr"&ddd$perc==2,], 
               aes(x = from, xend =  from+0.1, y = ystart+3, yend = ystart+3),
               arrow = arrow(angle = 90, length = unit(0.02, "inches"), ends = "first"),  
               col= "darkred", lty = 1, lwd = 0.4) +
  scale_linetype_manual(labels = c("90","95"), 
                        values = c("1" =3, "2" = 1), guide = "legend", name ="HWr") 

plotThr / plotProb + plot_layout(heights = c(1, 2), guides = "collect") &
  theme(legend.position = "top")
ggsave(filename = paste0("Plots/probPlot_2018.png"), dpi = 900, width = 190, height = 190, units = "mm")

