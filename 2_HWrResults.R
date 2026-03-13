# This function  build the database related to the emprirical-based approach  

source("functions/load_Packages.R")

# This function takes in argument:                                  
# obj: the resulting output of the call to ts2clm() function        
# year_start: first year of the observed time series                
# year_end: last year of the observed time series                   
# gap_days: to select to nember of consecutive days along which     
#           the observed temperature exceed the estimated threshold 
# summer_month: to select the summer months                         

extract_HWsw_summary <- function(obj, year_start, year_end, 
                                 gap_days = 3, summer_month = "MJJAS"){
  
  if(summer_month == "MJJAS"){ 
    ndays <- 153
    start <- paste0(seq(year_start, year_end, 1), "-05-01")
    end <- paste0(seq(year_start, year_end, 1), "-09-30")
  }
  nyears <- year_end - year_start + 1
  
  list_period <- data.frame(start = as.Date(start, format = "%Y-%m-%d"),
                            end  = as.Date(end, format = "%Y-%m-%d"))  
  
  data_smoothwindow <- data.frame(doy = rep(1:ndays, nyears), 
                                  years = rep(year_start : year_end, each = ndays),
                                  date = rep(NA, ndays * nyears),
                                  MaxT = rep(NA, ndays * nyears),
                                  seas = rep(NA, ndays * nyears),
                                  thresh = rep(NA, ndays * nyears),
                                  excedance = rep(NA, ndays *nyears),
                                  flag_excedance = rep(NA, ndays *nyears))
  
  for(j in 1 : nyears){
    data_smoothwindow$date[((j - 1) * ndays + 1) : (j * ndays)] <- obj$date[obj$date>=list_period$start[j] & obj$date <= list_period$end[j]]  
    #data_smoothwindow[((j - 1) * ndays + 1) : (j * ndays), c(4,5,6)] <-  obj[obj$date>=list_period$start[j] & obj$date <= list_period$end[j], 3:5]
    data_smoothwindow[((j - 1) * ndays + 1) : (j * ndays), c("MaxT","seas", "thresh")] <-  obj[obj$date>=list_period$start[j] & obj$date <= list_period$end[j], c("MaxT", "seas", "thresh") ]
  }
  colnames(data_smoothwindow)[which(colnames(data_smoothwindow) == "thresh")] <- "tresh"
  
  data_smoothwindow$date <- as.Date(data_smoothwindow$date)
  
  years_lab <- unique(data_smoothwindow$years)
  HW_year <- list()
  count <- 1
  for(j in years_lab){ 
    idx_year <- which(data_smoothwindow$years == j)
    data_smoothwindow$excedance[idx_year] <- data_smoothwindow$MaxT[idx_year] > data_smoothwindow$tresh[idx_year]
    tmp <- rle(data_smoothwindow$excedance[idx_year])
    data_smoothwindow$flag_excedance[idx_year] <- with(tmp, rep(lengths * values, lengths))
    HW_year[[count]] <- data_smoothwindow[idx_year,][ which(data_smoothwindow$flag_excedance[idx_year] >= gap_days), ]
    count <- count + 1
  }
  
  HW<- do.call("rbind", HW_year)
  
  HW_summary <- list()
  flag_OK <- TRUE
  count <- 1
  iter <- 1 
  if(nrow(HW) > 0){
    while(flag_OK){
      value <- HW$flag_excedance[count] 
      
      HW_summary[[iter]] <- data.frame(year = HW$years[count],
                                       doy_start = HW$doy[count],
                                       doy_end = HW$doy[count + value - 1], 
                                       ndays = value) 
      count <- count + value
      iter <- iter + 1
      if(count > nrow(HW) ){
        flag_OK <- FALSE
      }
    }
    out <- do.call("rbind", HW_summary)  
    
    
    nHWpy <- data.frame(year = year_start:year_end,
                        nHW = rep(0, year_end - year_start + 1)) 
    nHWpy$nHW[which(nHWpy$year %in% as.numeric(names(table(out$year))))] <-   table(out$year)
    nHWdpy <- data.frame(year = year_start:year_end,
                         nHWdays = rep(0, year_end - year_start + 1))
    
    extrNdays <- aggregate(ndays ~ year, data = out, FUN=cumsum)[,2]
    nHWdpy$nHWdays[which(nHWpy$year %in% as.numeric(names(table(out$year))))] <- unlist(lapply(1:length(extrNdays), function(x) extrNdays[[x]][length(extrNdays[[x]])]))  
    
    return(list(res = data_smoothwindow,
                summaryHW = out,
                noHWperyear = nHWpy,
                noHWdaysperyear = nHWdpy))
  } else {
    return(list(res = data_smoothwindow,
                summaryHW = data.frame(year = year_start:year_end, doy_start = NA, doy_end = NA, ndays = NA),
                noHWperyear = data.frame(year = year_start:year_end, nHW = 0),
                noHWdaysperyear = data.frame(year = year_start:year_end, nHWdays = 0)))
    
  }
}


# train set 1995 - 2018 ----

# Capriva
location = "Capriva"
load(paste0("data/dataFVG/dati_ALL", location, ".Rdata"))
dataC <- dati[dati$citta == location,]

res_C_up18 <- ts2clm3(data= dataC, x = date, y = MaxT,
                     climatologyPeriod = c("1995-01-01", "2018-12-31"),
                     windowHalfWidth = 7, smoothPercentile = FALSE, pctile = 90)
extract_HWsw_summary_C <- extract_HWsw_summary(obj = res_C_up18, 
                                               year_start = 1995, 
                                               year_end = 2018)
save(extract_HWsw_summary_C, file = "HeatwaveR/extract_HWsw_summary_C.Rdata")

res_C_up24 <- ts2clm3(data= dataC, x = date, y = MaxT,
                     climatologyPeriod = c("1995-01-01", "2018-10-07"),
                     windowHalfWidth = 7, smoothPercentile = FALSE, pctile = 90)

extract_HWsw_summary_C3 <- extract_HWsw_summary(obj = res_C_up24, 
                                                year_start = 2019, 
                                                year_end = 2024)
save(extract_HWsw_summary_C3, file = "HeatwaveR/extract_HWsw_summary_Cwith2024.Rdata")



# Trieste 
location = "Trieste"
load(paste0("data/dataFVG/dati_ALL", location, ".Rdata"))
dataT <- dati[dati$citta == location,]

res_T_up18 <- ts2clm3(data= dataT, x = date, y = MaxT,
                     climatologyPeriod = c("1995-01-01", "2018-12-31"),
                     windowHalfWidth = 7, smoothPercentile = FALSE, pctile = 90)


extract_HWsw_summary_T <- extract_HWsw_summary(obj = res_T_up18, 
                                               year_start = 1995, 
                                               year_end = 2018)
save(extract_HWsw_summary_T, file = "HeatwaveR/extract_HWsw_summary_T.Rdata")

res_T_up24 <- ts2clm3(data= dataT, x = date, y = MaxT,
                     climatologyPeriod = c("1995-01-01", "2018-10-07"),
                     windowHalfWidth = 7, smoothPercentile = FALSE, pctile = 90)

extract_HWsw_summary_T3 <- extract_HWsw_summary(obj = res_T_up24, 
                                                year_start = 2019, 
                                                year_end = 2024)
save(extract_HWsw_summary_T3, file = "HeatwaveR/extract_HWsw_summary_Twith2024.Rdata")

# Capriva
location = "Capriva"
load(paste0("data/dataFVG/dati_ALL", location, ".Rdata"))
dataC <- dati[dati$citta == location,]

res_C_up18 <- ts2clm3(data= dataC, x = date, y = MaxT,
                     climatologyPeriod = c("1995-01-01", "2018-12-31"),
                     windowHalfWidth = 7, smoothPercentile = FALSE, pctile = 95)
extract_HWsw_summary_C <- extract_HWsw_summary(obj = res_C_up18, year_start = 1995, year_end = 2018)
save(extract_HWsw_summary_C, file = "HeatwaveR/extract_HWsw_summary_C95.Rdata")


# Trieste 
location = "Trieste"
load(paste0("data/dataFVG/dati_ALL", location, ".Rdata"))
dataT <- dati[dati$citta == location,]

res_T_up18 <- ts2clm3(data= dataT, x = date, y = MaxT,
                     climatologyPeriod = c("1995-01-01", "2018-12-31"),
                     windowHalfWidth = 7, smoothPercentile = FALSE, pctile = 95)


extract_HWsw_summary_T <- extract_HWsw_summary(obj = res_T_up18, year_start = 1995, year_end = 2018)
save(extract_HWsw_summary_T, file = "HeatwaveR/extract_HWsw_summary_T95.Rdata")


year_start <- 1995
year_end <- 2018
year_end2 <- 2024


## heatwaveR ----
load("HeatwaveR/extract_HWsw_summary_C.Rdata")
load("HeatwaveR/extract_HWsw_summary_Cwith2024.Rdata")

load("HeatwaveR/extract_HWsw_summary_T.Rdata")
load("HeatwaveR/extract_HWsw_summary_Twith2024.Rdata")

database_Hwr_out <- data.frame(Type = "Hwr",
                               year = rep(((year_end+1):year_end2), times = 4),
                               measure = rep(c("HWn", "HWd"), each = 2*length(((year_end+1):year_end2))),
                               site = rep(c("Capriva", "Trieste", "Capriva", "Trieste"), each = length(((year_end+1):year_end2))),
                               Freq = c(extract_HWsw_summary_C3$noHWperyear$nHW[extract_HWsw_summary_C3$noHWperyear$year%in%((year_end+1):year_end2)],
                                        extract_HWsw_summary_T3$noHWperyear$nHW[extract_HWsw_summary_T3$noHWperyear$year%in%((year_end+1):year_end2)],
                                        extract_HWsw_summary_C3$noHWdaysperyear$nHWdays[extract_HWsw_summary_C3$noHWdaysperyear$year%in%((year_end+1):year_end2)],
                                        extract_HWsw_summary_T3$noHWdaysperyear$nHWdays[extract_HWsw_summary_T3$noHWdaysperyear$year%in%((year_end+1):year_end2)]),
                               period = "out")
database_Hwr_in <- data.frame(Type = "Hwr",
                              year = rep(((year_start):year_end), times = 4),
                              measure = rep(c("HWn", "HWd"), each = 2*length(((year_start):year_end))),
                              site = rep(c("Capriva", "Trieste", "Capriva", "Trieste"), each = length(((year_start):year_end))),
                              Freq = c(extract_HWsw_summary_C$noHWperyear$nHW,
                                       extract_HWsw_summary_T$noHWperyear$nHW,
                                       extract_HWsw_summary_C$noHWdaysperyear$nHWdays,
                                       extract_HWsw_summary_T$noHWdaysperyear$nHWdays),
                              period= "in")
database_Hwr <- rbind.data.frame(database_Hwr_in, database_Hwr_out)
save(database_Hwr, file = "HeatwaveR/database_Hwr.Rdata")

