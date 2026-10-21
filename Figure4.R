########figure4

library(ggalt)
library(dplyr)
library(tidyr)
library(aplot)
library(GGally)
library(scales)
library(ggbreak)
library(ggplot2)
library(apyramid)
library(showtext)
library(patchwork)
font_add("arial",'arial.ttf',bold = 'arialbd.ttf',italic = 'ariali.ttf',bolditalic = 'arialbi.ttf')
font_add("times",'times.ttf',bold = 'timesbd.ttf',italic = 'timesi.ttf',bolditalic = 'timesbi.ttf')
showtext_auto(enable=TRUE)
line_size <- 1/.pt
text_pt <- 8
dpi <- 300
text_size <- dpi/72 * text_pt
family <- "times"

find_interval <- function(min, max, 
                          n_range = c(0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,3,4,5,10,20,50,100,200,300,400,500,600)) {
  ideal_step <- (max - min) / 6.5 
  best_step <- min(n_range[which(n_range >= ideal_step)]) 
  return(best_step) 
}
####

save_path <-"./Putian_influenza_project/operation_files/DLNM_results_20251110/"

library(dlnm)
library(ggsci) 
library(tidyr)
library(tsModel)
library(ggplot2)
library(splines) 
library(reshape2)
library(grDevices)
library(tidyverse)
library(gridExtra)
library(patchwork)
library(RColorBrewer)
library(directlabels)

#dir.create(paste0(save_path,"DLNM_plot/original_data/df_selection/"), recursive = TRUE)
#dir.create("./Putian_influenza_project/operation_files/df_selection", recursive = TRUE)
dir.create(paste0(save_path,"DLNM_plot4/original_data/df_selection/"), recursive = TRUE)

modeldf_names <- c("respond","var","df_cb","df_time","qAIC","cb_P<0.1","times_P<0.1","5%_value","95%_value")
model_all <- data.frame(matrix(nrow=0,ncol=length(modeldf_names)))
colnames(model_all) <- modeldf_names; rm(modeldf_names)

peakdf_names <- c("respond","var","p_num","peak_x_range_min","peak_x_range_max","peak_x","peak_y","valley_x_range_min","valley_x_range_max","valley_x","valley_y")
picture_all <- data.frame(matrix(nrow=0,ncol=length(peakdf_names)))
colnames(picture_all) <- peakdf_names; rm(peakdf_names)

fqaic <- function(model){
  loglik <- sum(dpois(model$y,model$fitted.values,log = TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)    
}

cut_fc <- function(posi, nega, num) {
  if (nrow(posi)!=0){
    for (i in 1:nrow(posi)) {  posi$dd[i] <- round(posi[i+1,1]-posi[i,1],1)  }
    posi_cut <- which(posi$dd!=min(posi$dd,na.rm=TRUE),arr.ind = TRUE)
    if (length(posi_cut)==0) { 
      picture_all[nrow(picture_all)+1,] <- c(respon, var, num, range(posi[,1]),
                                             posi[,1][which(posi[,2]==max(posi[,2],na.rm=TRUE),arr.ind = TRUE)],
                                             round(posi[,2][which(posi[,2]==max(posi[,2],na.rm=TRUE),arr.ind = TRUE)],3),
                                             "-","-","-","-") 
    } else { 
      posi_cut <- c(0,posi_cut,nrow(posi))
      for (i in 1:c(length(posi_cut)-1)) { 
        picture_all[nrow(picture_all)+1,] <- c(respon, var, num,
                                               paste0(range(posi[,1][c(posi_cut[i]+1):posi_cut[i+1]])),
                                               posi[,1][c(posi_cut[i]+1):posi_cut[i+1]][which(
                                                 posi[,2][c(posi_cut[i]+1):posi_cut[i+1]]==max(posi[,2][c(posi_cut[i]+1):posi_cut[i+1]],na.rm=TRUE),
                                                 arr.ind = TRUE)],
                                               round(posi[,2][c(posi_cut[i]+1):posi_cut[i+1]][which(
                                                 posi[,2][c(posi_cut[i]+1):posi_cut[i+1]]==max(posi[,2][c(posi_cut[i]+1):posi_cut[i+1]],na.rm=TRUE),
                                                 arr.ind = TRUE)] ,3),
                                               "-","-","-","-")  
      } 
    } 
  } else {   
    picture_all[nrow(picture_all)+1,] <- c(respon,var,num,"NA","NA","-","-","-","-","-","-") 
  }
  if (nrow(nega)!=0){
    for (i in 1:nrow(nega)) {  nega$dd[i] <- round(nega[i+1,1]-nega[i,1],1 ) }
    nega_cut <- which(nega$dd > min(nega$dd,na.rm=TRUE),arr.ind = TRUE)
    if (length(nega_cut)==0) {
      picture_all[nrow(picture_all)+1,] <- c(respon,var,num,"-","-","-","-", range(nega[,1]),
                                             nega[,1][which(nega[,2]==min(nega[,2],na.rm=TRUE),arr.ind = TRUE)],
                                             round(nega[,2][which(nega[,2]==min(nega[,2],na.rm=TRUE),arr.ind = TRUE)],3 ))
    } else { 
      nega_cut <- c(0,nega_cut,nrow(nega))
      for (i in 1:c(length(nega_cut)-1)) { 
        picture_all[nrow(picture_all)+1,] <- c(respon,var,num,"-","-","-","-",
                                               paste0(range(nega[,1][c(nega_cut[i]+1):nega_cut[i+1]])),
                                               nega[,1][c(nega_cut[i]+1):nega_cut[i+1]][which(
                                                 nega[,2][c(nega_cut[i]+1):nega_cut[i+1]]==min(nega[,2][c(nega_cut[i]+1):nega_cut[i+1]],na.rm=TRUE),
                                                 arr.ind = TRUE)],
                                               round(nega[,2][c(nega_cut[i]+1):nega_cut[i+1]][which(
                                                 nega[,2][c(nega_cut[i]+1):nega_cut[i+1]]==min(nega[,2][c(nega_cut[i]+1):nega_cut[i+1]],na.rm=TRUE),
                                                 arr.ind = TRUE)] ,3))  
      }  
    }  
  } else { 
    picture_all[nrow(picture_all)+1,] <- c(respon,var,num,"-","-","-","-","NA","NA","-","-")    
  }  
}



library(dlnm)
library(readr)    # write_excel_csv
library(stringr)  # str_detect

# setwd("./Putian_influenza_project/operation_files/df_selection1")
wtcases <- read.csv("df_merge1_joined.csv", fileEncoding = "GBK")

## ---- DLNM analysis ----
wtcases$winddir <- as.numeric(wtcases$winddir)
wtcases <- wtcases[order(wtcases$Date),]
wtcases$time <- c(1:nrow(wtcases))
wtcases$winddir_use <- wtcases$winddir
wtcases$winddir_use[which(wtcases$winddir>337.5)] <- wtcases$winddir[which(wtcases$winddir>337.5)]-360

flu_types <- c("typeA", "typeB")
vars_list <- c("temp", "tempdiff", "humidity", "precip", "windspeed", "uvindex", "winddir_use", "solarradiation")

for (respon in flu_types) {
  for (var in vars_list) {
    M <- median(as.numeric(wtcases[,var]))
    lag <- 15
    qaic_model <- as.data.frame(matrix(nrow = 4, ncol = 5))
    
    for (i in 2:5) {
      for (n in 6:10) {
        try({
          cb_test = crossbasis(wtcases[,var], lag = lag,
                               argvar = list(fun="ns", df=i),
                               arglag = list(fun="poly", degree = 4))
          model_test <- glm(wtcases[,respon] ~ cb_test + ns(time,n*5) + log_weekly_detection + DOW + facial_coverings,
                            family = quasipoisson(), data = wtcases, na.action = na.omit)
          qaic_model[i-1,n-5] <- as.numeric(paste(fqaic(model_test)))
        }, silent = TRUE)
      }
    }
    colnames(qaic_model) <- c(6:10)
    rownames(qaic_model) <- c(2:5)
    
    write.csv(qaic_model,
              file = paste0(save_path, "DLNM_plot4/original_data/df_selection/",
                            respon, "-", var, "_df_selection_", Sys.Date(), ".csv"))
    
    df_cb  <- which(qaic_model == min(qaic_model, na.rm = TRUE), arr.ind = TRUE)[,1]+1
    df_time<- which(qaic_model == min(qaic_model, na.rm = TRUE), arr.ind = TRUE)[,2]+5
    
    cb <- crossbasis(wtcases[,var], lag = lag,
                     argvar = list(fun="ns", df=df_cb),
                     arglag = list(fun="poly", degree = 4))
    model <- glm(wtcases[,respon] ~ cb + ns(time,df_time*5) + log_weekly_detection + DOW + facial_coverings,
                 family = quasipoisson(), data = wtcases, na.action = na.omit)
    
    coefficients <- as.data.frame(summary(model)$coefficients)
    list_p0.1 <- rownames(coefficients[which(coefficients$`Pr(>|t|)`<0.1),])
    
    model_all[nrow(model_all)+1,] <- c(respon, var, df_cb, df_time,
                                       round(fqaic(model),3),
                                       sum(str_detect(list_p0.1,"cb")),
                                       sum(str_detect(list_p0.1,"time")), "", "")
    
    pred <- crosspred(cb, model, cen = M, bylag=0.2)
    rm(cb); rm(cb_test); rm(coefficients); rm(list_p0.1); rm(model); rm(model_test); rm(qaic_model); rm(i); rm(n)
    
    RRfit <- as.data.frame(pred$matRRfit)
    breaks_lines <- seq(min(RRfit), max(RRfit), by=((max(RRfit)-min(RRfit))/10))
    RRfit$var <- rownames(RRfit)
    map <- melt(RRfit)
    colnames(map) <- c("var","lag_day","value_")
    map$var <- as.numeric(map$var)
    map$lag_day <- as.numeric(gsub('[lag]', "", map$lag_day))
    
    RRall <- as.data.frame(pred$allRRfit)
    RRall$var <- as.numeric(rownames(RRall))
    colnames(RRall) <- c("RR","var")
    RRall.low <- as.data.frame(pred$allRRlow)
    RRall.low$var <- as.numeric(rownames(RRall.low))
    colnames(RRall.low) <- c("RR.low","var")
    RRall.high <- as.data.frame(pred$allRRhigh)
    RRall.high$var <- as.numeric(rownames(RRall.high))
    colnames(RRall.high) <- c("RR.high","var")
    comRR <- merge(merge(RRall, RRall.high, by = "var"), RRall.low, by = "var")
    rm(RRall); rm(RRall.low); rm(RRall.high)
    
    row_min = which.min(abs(quantile(wtcases[,var],0.05) - as.numeric(rownames(RRfit))))
    model_all[nrow(model_all),8] <- rownames(as.data.frame(pred$matRRfit)[row_min,])
    low.RR <- as.data.frame(t(as.data.frame(pred$matRRfit)[row_min,]))
    low.RR[,2] <- as.numeric(substr(rownames(low.RR),4,8))
    colnames(low.RR) <- c("RR","lag")
    low.RR.low <- as.data.frame(t(as.data.frame(pred$matRRlow)[row_min,]))
    low.RR.low[,2] <- as.numeric(substr(rownames(low.RR.low),4,8))
    colnames(low.RR.low) <- c("RR.low","lag")
    low.RR.up <- as.data.frame(t(as.data.frame(pred$matRRhigh)[row_min,]))
    low.RR.up[,2] <- as.numeric(substr(rownames(low.RR.up),4,8))
    colnames(low.RR.up) <- c("RR.up","lag")
    low <- merge(merge(low.RR, low.RR.up, by="lag"), low.RR.low, by="lag")
    
    row_max = which.min(abs(quantile(wtcases[,var],0.95) - as.numeric(rownames(RRfit))))
    model_all[nrow(model_all),9] <- rownames(as.data.frame(pred$matRRfit)[row_max,])
    high.RR <- as.data.frame(t(as.data.frame(pred$matRRfit)[row_max,]))
    high.RR[,2] <- as.numeric(substr(rownames(high.RR),4,8))
    colnames(high.RR) <- c("RR","lag")
    high.RR.low <- as.data.frame(t(as.data.frame(pred$matRRlow)[row_max,]))
    high.RR.low[,2] <- as.numeric(substr(rownames(high.RR.low),4,8))
    colnames(high.RR.low) <- c("RR.low","lag")
    high.RR.up <- as.data.frame(t(as.data.frame(pred$matRRhigh)[row_max,]))
    high.RR.up[,2] <- as.numeric(substr(rownames(high.RR.up),4,8))
    colnames(high.RR.up) <- c("RR.up","lag")
    high <- merge(merge(high.RR, high.RR.up, by="lag"), high.RR.low, by="lag")
    
    rm(low.RR, low.RR.up, low.RR.low, high.RR, high.RR.up, high.RR.low, df_cb, df_time, row_min, row_max)
    
    write.csv(comRR, file = paste0(save_path, "DLNM_plot4/original_data/", respon, "_", var, "_P2_", Sys.Date(), ".csv"), row.names = FALSE)
    write.csv(low,   file = paste0(save_path, "DLNM_plot4/original_data/", respon, "_", var, "_P3_", Sys.Date(), ".csv"), row.names = FALSE)
    write.csv(high,  file = paste0(save_path, "DLNM_plot4/original_data/", respon, "_", var, "_p4_", Sys.Date(), ".csv"), row.names = FALSE)
    
    posi <- comRR[which(comRR[,2]>1 & comRR[,3]>1 & comRR[,4]>1),]
    nega <- comRR[which(comRR[,2]<1 & comRR[,3]<1 & comRR[,4]<1),]
    num <- 2
    cut_fc(posi = posi, nega = nega, num = num)
    rm(posi, nega, num)
    
    posi <- low[which(low[,2]>1 & low[,3]>1 & low[,4]>1),]
    nega <- low[which(low[,2]<1 & low[,3]<1 & low[,4]<1),]
    num <- 3
    cut_fc(posi = posi, nega = nega, num = num)
    rm(posi, nega, num)
    
    posi <- high[which(high[,2]>1 & high[,3]>1 & high[,4]>1),]
    nega <- high[which(high[,2]<1 & high[,3]<1 & high[,4]<1),]
    num <- 4
    cut_fc(posi = posi, nega = nega, num = num)
    rm(posi, nega, num)
    
    rm(comRR, low, high)
  }
}

write_excel_csv(picture_all, file = paste0(save_path, "DLNM_plot4/original_data/DLNM_P2-P4_keypoints_summary_", Sys.Date(), ".csv"))
write_excel_csv(model_all,   file = paste0(save_path, "DLNM_plot4/original_data/DLNM_model_summary_", Sys.Date(), ".csv"))

# rm(M, lag, respon, var, vars_list, flu_types, map, pred, RRfit, breaks_lines, fqaic, cut_fc, model_all, picture_all)


df_path <- "./Putian_influenza_project/operation_files/DLNM_results_20251110/DLNM_plot4/"
### Fig4.5 DLNM typeA & typeB same color ####
theme <- theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(), 
        axis.line = element_line(color = "black", linewidth = line_size),
        axis.title.x = element_text(size=text_size,color="black",family = family),
        axis.title.y = element_text(size=text_size,color="black",family = family),
        axis.text  = element_text(size=text_size,color="black",family = family),
        plot.title = element_text(size = text_size, color = "black", family = family, hjust = 0.5,lineheight = 0.2),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt") )
DLNM_P1 <- function() {
  if (is.null(unit)) {  ggplot(data = datap2, aes(x = var, y = RR)) + 
      geom_hline(aes(yintercept = 1), colour="grey70", linetype="dashed") +
      geom_ribbon(aes(ymin = RR.low,  ymax = RR.high), alpha = 0.2, fill = "grey80") + 
      geom_line(linewidth = 0.8, color = "#612C1B") +
      labs(x=paste0(var_name), y="Cumulative relative risk")  +
      ggtitle(title_1) +
      scale_x_continuous(expand = c(0,0), limits = xlimt1, breaks = xbreak, labels = xlabs ) + 
      scale_y_continuous(expand = c(0,0), limits = ylimt1, breaks = ybreak1 ) +
      theme +
      theme(plot.title=element_text(margin=margin(b=-0.5, unit = "cm"))) } else {
        ggplot(data = datap2, aes(x = var, y = RR)) + 
          geom_hline(aes(yintercept = 1), colour="grey70", linetype="dashed") +
          geom_ribbon(aes(ymin = RR.low,  ymax = RR.high), alpha = 0.2, fill = "grey80") +
          geom_line(linewidth = 0.8, color = "#612C1B") +
          labs(x=paste0(var_name," (",unit,")"), y="Cumulative relative risk")  +
          ggtitle(title_1) +
          scale_x_continuous(expand = c(0,0), limits = xlimt1, breaks = xbreak, labels = xlabs ) + 
          scale_y_continuous(expand = c(0,0), limits = ylimt1, breaks = ybreak1 ) +
          theme +
          theme(plot.title=element_text(margin=margin(b=-0.5, unit = "cm"))) } }
DLNM_P2 <- function() {
  ggplot()+
    geom_hline(aes(yintercept = 1), colour="grey70", linetype="dashed") +
    geom_line(data = datapp, aes(x=lag, y=RR, group=group, color=group), lwd=1)+
    geom_ribbon(data = datapp, aes(ymin=RR.low, ymax=RR.up, group=group, fill=group, x=lag), alpha = 0.2)+
    theme_classic()+
    labs(x="Lag times (Days)",y="Relative risk", fill = NULL, color = NULL)+
    ggtitle(title_2) +
    scale_color_manual(values = c('#759CC7','#F28E2B'), guide = guide_legend(by.row = T, nrow = 1))+
    scale_fill_manual(values = c('#C9D7E5','#E7D1C5'))+
    scale_x_continuous(expand = c(0,0), limits = c(0, 15), breaks = seq(0, 15, 3) ) + 
    scale_y_continuous(expand = c(0,0), limits = ylimit2, breaks = ybreak2 ) +
    theme +
    theme(legend.position = c(0.5,1),
          legend.justification = c(0.5,1),
          legend.text = element_text(size = text_size, family=family, colour = "black", hjust = 0.5),
          legend.title = element_text(size = text_size, family=family, colour = "black", ),
          legend.key.height = unit(0.5, "cm"),
          legend.key.width = unit(1, "cm"),
          legend.background=element_rect(fill=rgb(1,1,1,alpha=0.001),colour=NA)) }
plot_list_1A <- list()
plot_list_1B <- list()
plot_list_2A <- list()
plot_list_2B <- list()



##################### Cumulative + Extreme (Same color for typeA & typeB)
#### Temperature typeA ####
type <- "typeA"; var <- "temp" ;var_name <- "Average temperature" ; unit <- "°C"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 24.1",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("14", unit) ; datap4$group <- paste0("30.5", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("14", unit),paste0("30.5", unit)))
datap2$RR.high[which(datap2$RR.high > 100)] <- 100
xlimt1 <- c(8, 32.5) ; xbreak <- seq(8, 32, 2); xlabs <- xbreak
ylimt1 <- c(0,100) ; ybreak1 <- seq(0,100, 10)
ylimit2 <- c(0.5, 6) ;ybreak2 <- seq(0.5, 6, 0.5)
plot_list_1A[[1]] <- DLNM_P1()
plot_list_2A[[1]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)

#### Temperature typeB ####
type <- "typeB"; var <- "temp" ;var_name <- "Average temperature" ; unit <- "°C"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 24.1",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("14", unit) ; datap4$group <- paste0("30.5", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("14", unit),paste0("30.5", unit)))
datap2$RR.high[which(datap2$RR.high > 44)] <- 44
xlimt1 <- c(8, 32.5) ; xbreak <- seq(8, 32, 2); xlabs <- xbreak
ylimt1 <- c(0,44) ; ybreak1 <- seq(0,44, 4)
ylimit2 <- c(0.6, 2) ;ybreak2 <- seq(0.6, 2, 0.2)
plot_list_1B[[1]] <- DLNM_P1()
plot_list_2B[[1]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)

#### Humidity typeA ####
type <- "typeA"; var <- "humidity" ;var_name <- "Humidity" ; unit <- "%"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 75.8",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("58", unit) ; datap4$group <- paste0("94", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("58", unit),paste0("94", unit)))
datap2$RR.high[which(datap2$RR.high > 30)] <- 30
#xlimt1 <- c(31, 99) ; xbreak <- c(31, 40, 50, 60, 70, 80, 90, 99); xlabs <- xbreak
xlimt1 <- c(41, 99) ; xbreak <- c(41, 50, 60, 70, 80, 90, 99); xlabs <- xbreak
ylimt1 <- c(0,30) ; ybreak1 <- seq(0,30, 3)
ylimit2 <- c(0.6, 1.8) ;ybreak2 <- seq(0.6, 1.8, 0.2)
plot_list_1A[[2]] <- DLNM_P1()
plot_list_2A[[2]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)

#### Humidity typeB ####
type <- "typeB"; var <- "humidity" ;var_name <- "Humidity" ; unit <- "%"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 75.8",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("58", unit) ; datap4$group <- paste0("94", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("58", unit),paste0("94", unit)))
datap2$RR.high[which(datap2$RR.high > 30)] <- 30
#xlimt1 <- c(31, 99) ; xbreak <- c(31, 40, 50, 60, 70, 80, 90, 99); xlabs <- xbreak
xlimt1 <- c(41, 99) ; xbreak <- c(41, 50, 60, 70, 80, 90, 99); xlabs <- xbreak
ylimt1 <- c(0,30) ; ybreak1 <- seq(0,30, 3)
ylimit2 <- c(0.5, 1.3) ;ybreak2 <- seq(0.5, 1.3, 0.1)
plot_list_1B[[2]] <- DLNM_P1()
plot_list_2B[[2]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)

#### Precipitation typeA ####
type <- "typeA"; var <- "precip" ;var_name <- "Precipitation" ; unit <- "mm"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 0 ",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("5 ", unit) ; datap4$group <- paste0("25 ", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("5 ", unit),paste0("25 ", unit)))
xlimt1 <- c(0,55) ; xbreak <- seq(0, 55, 5); xlabs <- xbreak
#ylimt1 <- c(0,5) ; ybreak1 <- seq(0,5, 0.5)
ylimt1 <- c(0,5) ; ybreak1 <- seq(0,5, 1)
ylimit2 <- c(0.5,1.4) ;ybreak2 <- seq(0.5, 1.4, 0.1)
plot_list_1A[[3]] <- DLNM_P1()
plot_list_2A[[3]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)

#### Precipitation typeB ####
type <- "typeB"; var <- "precip" ;var_name <- "Precipitation" ; unit <- "mm"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 0 ",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("5 ", unit) ; datap4$group <- paste0("25 ", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("5 ", unit),paste0("25 ", unit)))
#datap2$RR.high[which(datap2$RR.high > 1100)] <- 1100
datap2$RR.high[which(datap2$RR.high > 10)] <- 10
xlimt1 <- c(0,55) ; xbreak <- seq(0, 55, 5); xlabs <- xbreak
#ylimt1 <- c(0,1100) ; ybreak1 <- seq(0,1100, 100)
ylimt1 <- c(0,10) ; ybreak1 <- seq(0,10, 2)
ylimit2 <- c(0.6,2.8) ;ybreak2 <- seq(0.6, 2.8, 0.2)
plot_list_1B[[3]] <- DLNM_P1()
plot_list_2B[[3]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)

#### Solar radiation typeA ####
type <- "typeA"; var <- "solarradiation" ;var_name <- "Solar radiation" ; unit <- expression(W/m^2)
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T) #y = expression("Solar radiation (W/m"^2*")")
title_1 <- expression("Reference: 205.3 W/m"^2) 
title_2 <- expression("Solar radiation (Reference: 205.3 W/m"^2*")")
datap3$group <- "low" 
datap4$group <- "up"
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c("low","up"))
xlimt1 <- c(10, 340) ; xbreak <- seq(10,340, 30); xlabs <- xbreak
#ylimt1 <- c(0,35) ; ybreak1 <- seq(0,35, 5)
ylimt1 <- c(0,13) ; ybreak1 <- seq(0,13, 2)
ylimit2 <- c(0.5, 1.4) ;ybreak2 <- seq(0.5,1.4,0.1)
plot_list_1A[[4]] <- DLNM_P1() + labs(x=expression("Solar radiation (W/m"^2*")"))
plot_list_2A[[4]] <- DLNM_P2() + 
  scale_color_manual(values = c('#759CC7','#F28E2B'), guide = guide_legend(by.row = T, nrow = 1),
                     labels = c("low" = expression("75 W/m"^2),"up" = expression("305 W/m"^2)) ) +
  scale_fill_manual(values = c('#C9D7E5','#E7D1C5'),
                    labels = c("low" = expression("75 W/m"^2),"up" = expression("305 W/m"^2)) )
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)

#### Solar radiation typeB ####
type <- "typeB"; var <- "solarradiation" ;var_name <- "Solar radiation" ;  unit <- expression(W/m^2)
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <- expression("Reference: 205.3 W/m"^2) 
title_2 <- expression("Solar radiation ( Reference: 205.3 W/m"^2*")")
datap3$group <- "low" 
datap4$group <- "up"
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c("low","up"))
xlimt1 <- c(10, 340) ; xbreak <- seq(10,340, 30); xlabs <- xbreak
#ylimt1 <- c(0,7) ; ybreak1 <- seq(0,7,1)
ylimt1 <- c(0,13) ; ybreak1 <- seq(0,13,2)
ylimit2 <- c(0.5, 1.4) ;ybreak2 <- seq(0.5,1.4,0.1)
plot_list_1B[[4]] <- DLNM_P1() + labs(x=expression("Solar radiation (W/m"^2*")"))
plot_list_2B[[4]] <- DLNM_P2() + 
  scale_color_manual(values = c('#759CC7','#F28E2B'), guide = guide_legend(by.row = T, nrow = 1),
                     labels = c("low" = expression("75 W/m"^2),"up" = expression("305 W/m"^2)) ) +
  scale_fill_manual(values = c('#C9D7E5','#E7D1C5'),
                    labels = c("low" = expression("75 W/m"^2),"up" = expression("305 W/m"^2)) )
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)

#### Daily temperature range typeA ####
type <- "typeA"; var <- "tempdiff" ;var_name <- "Daily temperature ranges" ; unit <- "°C"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 8",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("3", unit) ; datap4$group <- paste0("13", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("3", unit),paste0("13", unit)))
##
datap2$RR.high[which(datap2$RR.high > 12)] <- 12
##
xlimt1 <- c(1, 16) ; xbreak <- seq(1, 16, 1); xlabs <- xbreak
#ylimt1 <- c(0,6) ; ybreak1 <- seq(0,6, 0.5)
ylimt1 <- c(0,12) ; ybreak1 <- seq(0,12, 2)
ylimit2 <- c(0.6, 1.3) ;ybreak2 <- seq(0.6, 1.3, 0.1)
plot_list_1A[[5]] <- DLNM_P1()
plot_list_2A[[5]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)

#### Daily temperature range typeB ####
type <- "typeB"; var <- "tempdiff" ;var_name <- "Daily temperature ranges" ; unit <- "°C"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 8",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("3", unit) ; datap4$group <- paste0("13", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("3", unit),paste0("13", unit)))
#datap2$RR.high[which(datap2$RR.high > 80)] <- 80
datap2$RR.high[which(datap2$RR.high > 18)] <- 18
xlimt1 <- c(1, 16) ; xbreak <- seq(1, 16, 1); xlabs <- xbreak
#ylimt1 <- c(0,80) ; ybreak1 <- seq(0,80, 10)
ylimt1 <- c(0,18) ; ybreak1 <- seq(0,18, 2)
ylimit2 <- c(0.6, 2.4) ;ybreak2 <- seq(0.6, 2.4, 0.2)
plot_list_1B[[5]] <- DLNM_P1()
plot_list_2B[[5]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)

#### UV index typeA ####
type <- "typeA"; var <- "uvindex" ;var_name <- "Ultraviolet index" ; unit <- NULL
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 7",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("3 ", unit) ; datap4$group <- paste0("10 ", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("3 ", unit),paste0("10 ", unit)))
xlimt1 <- c(0,10) ; xbreak <- seq(0, 10, 1); xlabs <- xbreak
#ylimt1 <- c(0,12) ; ybreak1 <- seq(0, 12, 1)
ylimt1 <- c(0,5) ; ybreak1 <- seq(0, 5, 1)
ylimit2 <- c(0.4, 1.5) ;ybreak2 <- seq(0.4, 1.5, 0.1)
plot_list_1A[[6]] <- DLNM_P1()
plot_list_2A[[6]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)

#### UV index typeB ####
type <- "typeB"; var <- "uvindex" ;var_name <- "Ultraviolet index" ; unit <- NULL
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 7",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("3 ", unit) ; datap4$group <- paste0("10 ", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("3 ", unit),paste0("10 ", unit)))
#datap2$RR.high[which(datap2$RR.high > 100)] <- 100
datap2$RR.high[which(datap2$RR.high > 40)] <- 40
xlimt1 <- c(0,10) ; xbreak <- seq(0, 10, 1); xlabs <- xbreak
#ylimt1 <- c(0,100) ; ybreak1 <- seq(0, 100, 10)
ylimt1 <- c(0,40) ; ybreak1 <- seq(0, 40, 5)
ylimit2 <- c(0.4, 1.5) ;ybreak2 <- seq(0.4, 1.5, 0.1)
plot_list_1B[[6]] <- DLNM_P1()
plot_list_2B[[6]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)

#### Wind speed typeA ####
type <- "typeA"; var <- "windspeed" ;var_name <- "Wind speed" ; unit <- "km/h"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 18 ",unit) 
title_2 <- paste0(var_name, " (",title_1,")")
datap3$group <- paste0("10 ", unit) ; datap4$group <- paste0("26 ", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("10 ", unit),paste0("26 ", unit)))
##
datap2$RR.high[which(datap2$RR.high > 25)] <- 25
##
xlimt1 <- c(8, 34) ; xbreak <- seq(8, 34, 2); xlabs <- xbreak
ylimt1 <- c(0,25) ; ybreak1 <- seq(0,25, 5)
ylimit2 <- c(0.4, 1.8) ;ybreak2 <- seq(0.4, 1.8, 0.2)
plot_list_1A[[7]] <- DLNM_P1()
plot_list_2A[[7]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)

#### Wind speed typeB ####
type <- "typeB"; var <- "windspeed" ;var_name <- "Wind speed" ; unit <- "km/h"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 18 ",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("10 ", unit) ; datap4$group <- paste0("26 ", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("10 ", unit),paste0("26 ", unit)))
xlimt1 <- c(8, 34) ; xbreak <- seq(8, 34, 2); xlabs <- xbreak
ylimt1 <- c(0,4) ; ybreak1 <- seq(0,4, 0.5)
ylimit2 <- c(0.4, 1.8) ;ybreak2 <- seq(0.4, 1.8, 0.2)
plot_list_1B[[7]] <- DLNM_P1()
plot_list_2B[[7]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)

#### Wind direction typeA ####
type <- "typeA"; var <- "winddir_use" ;var_name <- "Wind direction" ; unit <- NULL
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: SW",unit) 
title_2 <- paste0("Wind direction (Reference: SW",unit,")") 
xlimt1 <- c(-20,330) ; xbreak <- c(0,45,90,135,180,225,270,315); xlabs <- c("N","NE","E","SE","S","SW","W","NW")
#ylimt1 <- c(0,6) ; ybreak1 <- seq(0,6,1)
ylimt1 <- c(0,16) ; ybreak1 <- seq(0,16,4)
datap3$group <- "NW" ; datap4$group <- "E"
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c("NW","E"))
ylimit2 <- c(0.2, 1.8) ;ybreak2 <- seq(0.2, 1.8, 0.2)
plot_list_1A[[8]] <- DLNM_P1()
plot_list_2A[[8]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)

#### Wind direction typeB ####
type <- "typeB"; var <- "winddir_use" ;var_name <- "Wind direction" ; unit <- NULL
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: SW",unit) 
title_2 <- paste0("Wind direction (Reference: SW",unit,")") 
xlimt1 <- c(-20,330) ; xbreak <- c(0,45,90,135,180,225,270,315); xlabs <- c("N","NE","E","SE","S","SW","W","NW")
#ylimt1 <- c(0,500) ; ybreak1 <- seq(0,500,50)
ylimt1 <- c(0,15) ; ybreak1 <- seq(0,15,3)
#datap2$RR.high[which(datap2$RR.high > 500)] <- 500
datap2$RR.high[which(datap2$RR.high > 15)] <- 15
datap3$group <- "NW" ; datap4$group <- "E"
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c("NW","E"))
ylimit2 <- c(0, 2.4) ;ybreak2 <- seq(0, 2.4, 0.2)
plot_list_1B[[8]] <- DLNM_P1()
plot_list_2B[[8]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)

#### Combine plots ####
wrap_plots(plot_list_1A) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  theme(panel.spacing = unit(0, "mm"),
        legend.position = "bottom",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        plot.tag = element_text(size = text_size, colour = "black", family = family),
        plot.margin = margin(t = 1, r = 5, b = 1, l = 1, unit = "pt"))
ggsave(paste0(df_path,"DLNM_cumulative_risk_typeA.tiff"), height = 9, width = 7,units = "in",bg = NULL, dpi=dpi)

wrap_plots(plot_list_1B) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  theme(panel.spacing = unit(0, "mm"),
        legend.position = "bottom",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        plot.tag = element_text(size = text_size, colour = "black", family = family),
        plot.margin = margin(t = 1, r = 5, b = 1, l = 1, unit = "pt"))
ggsave(paste0(df_path,"DLNM_cumulative_risk_typeB.tiff"), height = 9, width = 7,units = "in",bg = NULL, dpi=dpi)

wrap_plots(plot_list_2A) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A") & 
  theme(panel.spacing = unit(0, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        plot.tag = element_text(size = text_size, colour = "black", family = family),
        plot.margin = margin(t = 1, r = 5, b = 1, l = 1, unit = "pt"))
ggsave(paste0(df_path,"DLNM_extreme_meteorology_typeA.tiff"), height = 9, width = 7,units = "in",bg = NULL, dpi=dpi)

wrap_plots(plot_list_2B) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A") & 
  theme(panel.spacing = unit(0, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        plot.tag = element_text(size = text_size, colour = "black", family = family),
        plot.margin = margin(t = 1, r = 5, b = 1, l = 1, unit = "pt"))
ggsave(paste0(df_path,"DLNM_extreme_meteorology_typeB.tiff"), height = 9, width = 7,units = "in",bg = NULL, dpi=dpi)

wrap_plots(plot_list_1A) +
  plot_layout(ncol = 2, guides = "collect") +
  theme(panel.spacing = unit(0, "mm"),
        legend.position = "bottom",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        plot.tag = element_text(size = text_size, colour = "black", family = family),
        plot.margin = margin(t = 1, r = 5, b = 1, l = 1, unit = "pt"))
ggsave(paste0(df_path,"DLNM_cumulative_risk_typeA_no_tag.tiff"), height = 9, width = 7,units = "in",bg = NULL, dpi=dpi)

wrap_plots(plot_list_1B) +
  plot_layout(ncol = 2, guides = "collect") +
  theme(panel.spacing = unit(0, "mm"),
        legend.position = "bottom",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        plot.tag = element_text(size = text_size, colour = "black", family = family),
        plot.margin = margin(t = 1, r = 5, b = 1, l = 1, unit = "pt"))
ggsave(paste0(df_path,"DLNM_cumulative_risk_typeB_no_tag.tiff"), height = 9, width = 7,units = "in",bg = NULL, dpi=dpi)

wrap_plots(plot_list_2A) +
  plot_layout(ncol = 2) +
  theme(panel.spacing = unit(0, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        plot.tag = element_text(size = text_size, colour = "black", family = family),
        plot.margin = margin(t = 1, r = 5, b = 1, l = 1, unit = "pt"))
ggsave(paste0(df_path,"DLNM_extreme_meteorology_typeA_no_tag.tiff"), height = 9, width = 7,units = "in",bg = NULL, dpi=dpi)

wrap_plots(plot_list_2B) +
  plot_layout(ncol = 2) +
  theme(panel.spacing = unit(0, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        plot.tag = element_text(size = text_size, colour = "black", family = family),
        plot.margin = margin(t = 1, r = 5, b = 1, l = 1, unit = "pt"))
ggsave(paste0(df_path,"DLNM_extreme_meteorology_typeB_no_tag.tiff"), height = 9, width = 7,units = "in",bg = NULL, dpi=dpi)
#rm(theme, plot_list_1A,plot_list_1B,plot_list_2A,plot_list_2B,DLNM_P1, DLNM_P2)

### Fig4.5 DLNM typeA & typeB different color ####
theme <- theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(), 
        axis.line = element_line(color = "black", linewidth = line_size),
        axis.title.x = element_text(size=text_size,color="black",family = family),
        axis.title.y = element_text(size=text_size,color="black",family = family),
        axis.text  = element_text(size=text_size,color="black",family = family),
        plot.title = element_text(size = text_size, color = "black", family = family, hjust = 0.5,lineheight = 0.2),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt") )

DLNM_P1 <- function() {
  if (type == "typeA") {color <- "#A32A31"} else {color <- "#551C6E"}
  if (is.null(unit)) {  ggplot(data = datap2, aes(x = var, y = RR)) + 
      geom_hline(aes(yintercept = 1), colour="grey70", linetype="dashed") +
      geom_ribbon(aes(ymin = RR.low,  ymax = RR.high), alpha = 0.2, fill = "grey80") +
      geom_line(linewidth = 0.8, color = color) +
      labs(x=paste0(var_name), y="Cumulative relative risk")  +
      ggtitle(title_1) +
      scale_x_continuous(expand = c(0,0), limits = xlimt1, breaks = xbreak, labels = xlabs ) + 
      scale_y_continuous(expand = c(0,0), limits = ylimt1, breaks = ybreak1 ) +
      theme +
      theme(plot.title=element_text(margin=margin(b=-0.5, unit = "cm"))) } else {
        ggplot(data = datap2, aes(x = var, y = RR)) + 
          geom_hline(aes(yintercept = 1), colour="grey70", linetype="dashed") +
          geom_ribbon(aes(ymin = RR.low,  ymax = RR.high), alpha = 0.2, fill = "grey80") +
          geom_line(linewidth = 0.8, color = color) +
          labs(x=paste0(var_name," (",unit,")"), y="Cumulative relative risk")  +
          ggtitle(title_1) +
          scale_x_continuous(expand = c(0,0), limits = xlimt1, breaks = xbreak, labels = xlabs ) + 
          scale_y_continuous(expand = c(0,0), limits = ylimt1, breaks = ybreak1 ) +
          theme +
          theme(plot.title=element_text(margin=margin(b=-0.5, unit = "cm")))
      }
}

DLNM_P2 <- function() {
  if (type == "typeA") {
    color1 <- c('#759CC7','#F28E2B')
    color2 <- c('#C9D7E5','#E7D1C5')
  } else {
    color1 <- c('#62B197', '#E1C855' )
    color2 <- c('#CDE0C7', '#F7EBC6' )
  }
  ggplot()+
    geom_hline(aes(yintercept = 1), colour="grey70", linetype="dashed") +
    geom_line(data = datapp, aes(x=lag, y=RR, group=group, color=group), lwd=1)+
    geom_ribbon(data = datapp, aes(ymin=RR.low, ymax=RR.up, group=group, fill=group, x=lag), alpha = 0.2)+
    theme_classic()+
    labs(x="Lag times (Days)",y="Relative risk", fill = NULL, color = NULL)+
    ggtitle(title_2) +
    scale_color_manual(values = color1, guide = guide_legend(by.row = T, nrow = 1))+
    scale_fill_manual(values = color2 )+
    scale_x_continuous(expand = c(0,0), limits = c(0, 15), breaks = seq(0, 15, 3) ) + 
    scale_y_continuous(expand = c(0,0), limits = ylimit2, breaks = ybreak2 ) +
    theme +
    theme(legend.position = c(0.5,1),
          legend.justification = c(0.5,1),
          legend.text = element_text(size = text_size, family=family, colour = "black", hjust = 0.5),
          legend.title = element_text(size = text_size, family=family, colour = "black", ),
          legend.key.height = unit(0.5, "cm"),
          legend.key.width = unit(1, "cm"),
          legend.background=element_rect(fill=rgb(1,1,1,alpha=0.001),colour=NA))
}

plot_list_1A <- list()
plot_list_1B <- list()
plot_list_2A <- list()
plot_list_2B <- list()

################################

##################Cumulative + Extreme Stratification
###Cumulative Stratification
####Temperature A####
type <- "typeA"; var <- "temp" ;var_name <- "Average temperature" ; unit <- "°C"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 24.1",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("14", unit) ; datap4$group <- paste0("30.5", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("14", unit),paste0("30.5", unit)))
datap2$RR.high[which(datap2$RR.high > 100)] <- 100
xlimt1 <- c(8, 32.5) ; xbreak <- seq(8, 32, 2); xlabs <- xbreak
ylimt1 <- c(0,100) ; ybreak1 <- seq(0,100, 10)
ylimit2 <- c(0.5, 6) ;ybreak2 <- seq(0.5, 6, 0.5)
plot_list_1A[[1]] <- DLNM_P1()
plot_list_2A[[1]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)
####Temperature B####
type <- "typeB"; var <- "temp" ;var_name <- "Average temperature" ; unit <- "°C"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 24.1",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("14", unit) ; datap4$group <- paste0("30.5", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("14", unit),paste0("30.5", unit)))
datap2$RR.high[which(datap2$RR.high > 44)] <- 44
xlimt1 <- c(8, 32.5) ; xbreak <- seq(8, 32, 2); xlabs <- xbreak
ylimt1 <- c(0,44) ; ybreak1 <- seq(0,44, 4)
ylimit2 <- c(0.6, 2) ;ybreak2 <- seq(0.6, 2, 0.2)
plot_list_1B[[1]] <- DLNM_P1()
plot_list_2B[[1]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)
####Humidity A####
type <- "typeA"; var <- "humidity" ;var_name <- "Humidity" ; unit <- "%"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 75.8",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("58", unit) ; datap4$group <- paste0("94", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("58", unit),paste0("94", unit)))
datap2$RR.high[which(datap2$RR.high > 30)] <- 30
#xlimt1 <- c(31, 99) ; xbreak <- c(31, 40, 50, 60, 70, 80, 90, 99); xlabs <- xbreak
xlimt1 <- c(41, 99) ; xbreak <- c(41, 50, 60, 70, 80, 90, 99); xlabs <- xbreak
ylimt1 <- c(0,30) ; ybreak1 <- seq(0,30, 3)
ylimit2 <- c(0.6, 1.8) ;ybreak2 <- seq(0.6, 1.8, 0.2)
plot_list_1A[[2]] <- DLNM_P1()
plot_list_2A[[2]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)
####Humidity B####
type <- "typeB"; var <- "humidity" ;var_name <- "Humidity" ; unit <- "%"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 75.8",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("58", unit) ; datap4$group <- paste0("94", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("58", unit),paste0("94", unit)))
datap2$RR.high[which(datap2$RR.high > 30)] <- 30
#xlimt1 <- c(31, 99) ; xbreak <- c(31, 40, 50, 60, 70, 80, 90, 99); xlabs <- xbreak
xlimt1 <- c(41, 99) ; xbreak <- c(41, 50, 60, 70, 80, 90, 99); xlabs <- xbreak
ylimt1 <- c(0,30) ; ybreak1 <- seq(0,30, 3)
ylimit2 <- c(0.5, 1.3) ;ybreak2 <- seq(0.5, 1.3, 0.1)
plot_list_1B[[2]] <- DLNM_P1()
plot_list_2B[[2]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)
####Precipitation A####
type <- "typeA"; var <- "precip" ;var_name <- "Precipitation" ; unit <- "mm"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 0 ",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("5 ", unit) ; datap4$group <- paste0("25 ", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("5 ", unit),paste0("25 ", unit)))
xlimt1 <- c(0,55) ; xbreak <- seq(0, 55, 5); xlabs <- xbreak
#ylimt1 <- c(0,5) ; ybreak1 <- seq(0,5, 0.5)
ylimt1 <- c(0,5) ; ybreak1 <- seq(0,5, 1)
ylimit2 <- c(0.5,1.4) ;ybreak2 <- seq(0.5, 1.4, 0.1)
plot_list_1A[[3]] <- DLNM_P1()
plot_list_2A[[3]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)
####Precipitation B####
type <- "typeB"; var <- "precip" ;var_name <- "Precipitation" ; unit <- "mm"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 0 ",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("5 ", unit) ; datap4$group <- paste0("25 ", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("5 ", unit),paste0("25 ", unit)))
#datap2$RR.high[which(datap2$RR.high > 1100)] <- 1100
datap2$RR.high[which(datap2$RR.high > 10)] <- 10
xlimt1 <- c(0,55) ; xbreak <- seq(0, 55, 5); xlabs <- xbreak
#ylimt1 <- c(0,1100) ; ybreak1 <- seq(0,1100, 100)
ylimt1 <- c(0,10) ; ybreak1 <- seq(0,10, 2)
ylimit2 <- c(0.6,2.8) ;ybreak2 <- seq(0.6, 2.8, 0.2)
plot_list_1B[[3]] <- DLNM_P1()
plot_list_2B[[3]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)
####Solar Radiation A####
type <- "typeA"; var <- "solarradiation" ;var_name <- "Solar radiation" ; unit <- expression(W/m^2)
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T) #y = expression("Solar radiation (W/m"^2*")")
title_1 <- expression("Reference: 205.3 W/m"^2) 
title_2 <- expression("Solar radiation (Reference: 205.3 W/m"^2*")")
datap3$group <- "low" 
datap4$group <- "up"
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c("low","up"))
xlimt1 <- c(10, 340) ; xbreak <- seq(10,340, 30); xlabs <- xbreak
#ylimt1 <- c(0,35) ; ybreak1 <- seq(0,35, 5)
ylimt1 <- c(0,13) ; ybreak1 <- seq(0,13, 2)
ylimit2 <- c(0.5, 1.4) ;ybreak2 <- seq(0.5,1.4,0.1)
plot_list_1A[[4]] <- DLNM_P1() + labs(x=expression("Solar radiation (W/m"^2*")"))
plot_list_2A[[4]] <- DLNM_P2() + 
  scale_color_manual(values = c('#759CC7','#F28E2B'), guide = guide_legend(by.row = T, nrow = 1),
                     labels = c("low" = expression("75 W/m"^2),"up" = expression("305 W/m"^2)) ) +
  scale_fill_manual(values = c('#C9D7E5','#E7D1C5'),
                    labels = c("low" = expression("75 W/m"^2),"up" = expression("305 W/m"^2)) )
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)
####Solar Radiation B####
type <- "typeB"; var <- "solarradiation" ;var_name <- "Solar radiation" ;  unit <- expression(W/m^2)
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <- expression("Reference: 205.3 W/m"^2) 
title_2 <- expression("Solar radiation ( Reference: 205.3 W/m"^2*")")
datap3$group <- "low" 
datap4$group <- "up"
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c("low","up"))
xlimt1 <- c(10, 340) ; xbreak <- seq(10,340, 30); xlabs <- xbreak
#ylimt1 <- c(0,7) ; ybreak1 <- seq(0,7,1)
ylimt1 <- c(0,13) ; ybreak1 <- seq(0,13,2)
ylimit2 <- c(0.5, 1.4) ;ybreak2 <- seq(0.5,1.4,0.1)
plot_list_1B[[4]] <- DLNM_P1() + labs(x=expression("Solar radiation (W/m"^2*")"))
plot_list_2B[[4]] <- DLNM_P2() + 
  scale_color_manual(values = c('#759CC7','#F28E2B'), guide = guide_legend(by.row = T, nrow = 1),
                     labels = c("low" = expression("75 W/m"^2),"up" = expression("305 W/m"^2)) ) +
  scale_fill_manual(values = c('#C9D7E5','#E7D1C5'),
                    labels = c("low" = expression("75 W/m"^2),"up" = expression("305 W/m"^2)) )
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)
####Diurnal Temperature Range A####
type <- "typeA"; var <- "tempdiff" ;var_name <- "Daily temperature ranges" ; unit <- "°C"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 8",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("3", unit) ; datap4$group <- paste0("13", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("3", unit),paste0("13", unit)))
##
datap2$RR.high[which(datap2$RR.high > 12)] <- 12
##
xlimt1 <- c(1, 16) ; xbreak <- seq(1, 16, 1); xlabs <- xbreak
#ylimt1 <- c(0,6) ; ybreak1 <- seq(0,6, 0.5)
ylimt1 <- c(0,12) ; ybreak1 <- seq(0,12, 2)
ylimit2 <- c(0.6, 1.3) ;ybreak2 <- seq(0.6, 1.3, 0.1)
plot_list_1A[[5]] <- DLNM_P1()
plot_list_2A[[5]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)
####Diurnal Temperature Range B####
type <- "typeB"; var <- "tempdiff" ;var_name <- "Daily temperature ranges" ; unit <- "°C"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 8",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("3", unit) ; datap4$group <- paste0("13", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("3", unit),paste0("13", unit)))
#datap2$RR.high[which(datap2$RR.high > 80)] <- 80
datap2$RR.high[which(datap2$RR.high > 18)] <- 18
xlimt1 <- c(1, 16) ; xbreak <- seq(1, 16, 1); xlabs <- xbreak
#ylimt1 <- c(0,80) ; ybreak1 <- seq(0,80, 10)
ylimt1 <- c(0,18) ; ybreak1 <- seq(0,18, 2)
ylimit2 <- c(0.6, 2.4) ;ybreak2 <- seq(0.6, 2.4, 0.2)
plot_list_1B[[5]] <- DLNM_P1()
plot_list_2B[[5]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)
####UV Index A####
type <- "typeA"; var <- "uvindex" ;var_name <- "Ultraviolet index" ; unit <- NULL
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 7",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("3 ", unit) ; datap4$group <- paste0("10 ", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("3 ", unit),paste0("10 ", unit)))
xlimt1 <- c(0,10) ; xbreak <- seq(0, 10, 1); xlabs <- xbreak
#ylimt1 <- c(0,12) ; ybreak1 <- seq(0, 12, 1)
ylimt1 <- c(0,5) ; ybreak1 <- seq(0, 5, 1)
ylimit2 <- c(0.4, 1.5) ;ybreak2 <- seq(0.4, 1.5, 0.1)
plot_list_1A[[6]] <- DLNM_P1()
plot_list_2A[[6]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)
####UV Index B####
type <- "typeB"; var <- "uvindex" ;var_name <- "Ultraviolet index" ; unit <- NULL
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 7",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("3 ", unit) ; datap4$group <- paste0("10 ", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("3 ", unit),paste0("10 ", unit)))
#datap2$RR.high[which(datap2$RR.high > 100)] <- 100
datap2$RR.high[which(datap2$RR.high > 40)] <- 40
xlimt1 <- c(0,10) ; xbreak <- seq(0, 10, 1); xlabs <- xbreak
#ylimt1 <- c(0,100) ; ybreak1 <- seq(0, 100, 10)
ylimt1 <- c(0,40) ; ybreak1 <- seq(0, 40, 5)
ylimit2 <- c(0.4, 1.5) ;ybreak2 <- seq(0.4, 1.5, 0.1)
plot_list_1B[[6]] <- DLNM_P1()
plot_list_2B[[6]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)
####Wind Speed A####
type <- "typeA"; var <- "windspeed" ;var_name <- "Wind speed" ; unit <- "km/h"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 18 ",unit) 
title_2 <- paste0(var_name, " (",title_1,")")
datap3$group <- paste0("10 ", unit) ; datap4$group <- paste0("26 ", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("10 ", unit),paste0("26 ", unit)))
##
datap2$RR.high[which(datap2$RR.high > 25)] <- 25
##
xlimt1 <- c(8, 34) ; xbreak <- seq(8, 34, 2); xlabs <- xbreak
ylimt1 <- c(0,25) ; ybreak1 <- seq(0,25, 5)
ylimit2 <- c(0.4, 1.8) ;ybreak2 <- seq(0.4, 1.8, 0.2)
plot_list_1A[[7]] <- DLNM_P1()
plot_list_2A[[7]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)
####Wind Speed B####
type <- "typeB"; var <- "windspeed" ;var_name <- "Wind speed" ; unit <- "km/h"
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: 18 ",unit) 
title_2 <- paste0(var_name, " (",title_1,")") 
datap3$group <- paste0("10 ", unit) ; datap4$group <- paste0("26 ", unit)
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c(paste0("10 ", unit),paste0("26 ", unit)))
xlimt1 <- c(8, 34) ; xbreak <- seq(8, 34, 2); xlabs <- xbreak
ylimt1 <- c(0,4) ; ybreak1 <- seq(0,4, 0.5)
ylimit2 <- c(0.4, 1.8) ;ybreak2 <- seq(0.4, 1.8, 0.2)
plot_list_1B[[7]] <- DLNM_P1()
plot_list_2B[[7]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)
####Wind Direction A####
type <- "typeA"; var <- "winddir_use" ;var_name <- "Wind direction" ; unit <- NULL
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: SW",unit) 
title_2 <- paste0("Wind direction (Reference: SW",unit,")") 
xlimt1 <- c(-20,330) ; xbreak <- c(0,45,90,135,180,225,270,315); xlabs <- c("N","NE","E","SE","S","SW","W","NW")
#ylimt1 <- c(0,6) ; ybreak1 <- seq(0,6,1)
ylimt1 <- c(0,16) ; ybreak1 <- seq(0,16,4)
datap3$group <- "NW" ; datap4$group <- "E"
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c("NW","E"))
ylimit2 <- c(0.2, 1.8) ;ybreak2 <- seq(0.2, 1.8, 0.2)
plot_list_1A[[8]] <- DLNM_P1()
plot_list_2A[[8]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)
####Wind Direction B####
type <- "typeB"; var <- "winddir_use" ;var_name <- "Wind direction" ; unit <- NULL
datap2 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P2"), full.names = T), header = T)
datap3 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_P3"), full.names = T), header = T)
datap4 <- read.csv(list.files(paste0(df_path,"raw_data/"), pattern = paste0(type,"_",var,"_p4"), full.names = T), header = T)
title_1 <-  paste0("Reference: SW",unit) 
title_2 <- paste0("Wind direction (Reference: SW",unit,")") 
xlimt1 <- c(-20,330) ; xbreak <- c(0,45,90,135,180,225,270,315); xlabs <- c("N","NE","E","SE","S","SW","W","NW")
#ylimt1 <- c(0,500) ; ybreak1 <- seq(0,500,50)
ylimt1 <- c(0,15) ; ybreak1 <- seq(0,15,3)
#datap2$RR.high[which(datap2$RR.high > 500)] <- 500
datap2$RR.high[which(datap2$RR.high > 15)] <- 15
datap3$group <- "NW" ; datap4$group <- "E"
datapp <- rbind(datap3, datap4)
datapp$group <- factor(datapp$group, levels = c("NW","E"))
ylimit2 <- c(0, 2.4) ;ybreak2 <- seq(0, 2.4, 0.2)
plot_list_1B[[8]] <- DLNM_P1()
plot_list_2B[[8]] <- DLNM_P2()
rm(type);rm(var);rm(var_name);rm(unit);rm(datap2);rm(datap3);rm(datap4);rm(title_1);rm(title_2)
rm(xlimt1);rm(xbreak);rm(xlabs);rm(ylimt1);rm(ybreak1);rm(datapp);rm(ylimit2);rm(ybreak2)
####Combine Plots####
wrap_plots(plot_list_1A) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  theme(panel.spacing = unit(0, "mm"),
        legend.position = "bottom",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        plot.tag = element_text(size = text_size, colour = "black", family = family),
        plot.margin = margin(t = 1, r = 5, b = 1, l = 1, unit = "pt"))
ggsave(paste0(df_path,"DLNM_cumulative_risk_InfluenzaA.tiff"), height = 9, width = 7,units = "in",bg = NULL, dpi=dpi)
wrap_plots(plot_list_1B) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  theme(panel.spacing = unit(0, "mm"),
        legend.position = "bottom",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        plot.tag = element_text(size = text_size, colour = "black", family = family),
        plot.margin = margin(t = 1, r = 5, b = 1, l = 1, unit = "pt"))
ggsave(paste0(df_path,"DLNM_cumulative_risk_InfluenzaB.tiff"), height = 9, width = 7,units = "in",bg = NULL, dpi=dpi)
wrap_plots(plot_list_2A) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A") & 
  theme(panel.spacing = unit(0, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        plot.tag = element_text(size = text_size, colour = "black", family = family),
        plot.margin = margin(t = 1, r = 5, b = 1, l = 1, unit = "pt"))

ggsave(paste0(df_path,"DLNM_cumulative_risk_InfluenzaA_no_label.tiff"), height = 9, width = 7,units = "in",bg = NULL, dpi=dpi)
wrap_plots(plot_list_1B) +
  plot_layout(ncol = 2, guides = "collect") +
  theme(panel.spacing = unit(0, "mm"),
        legend.position = "bottom",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        plot.tag = element_text(size = text_size, colour = "black", family = family),
        plot.margin = margin(t = 1, r = 5, b = 1, l = 1, unit = "pt"))
ggsave(paste0(df_path,"DLNM_cumulative_risk_InfluenzaB_no_label.tiff"), height = 9, width = 7,units = "in",bg = NULL, dpi=dpi)
wrap_plots(plot_list_2A) +
  plot_layout(ncol = 2) +
  theme(panel.spacing = unit(0, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        plot.tag = element_text(size = text_size, colour = "black", family = family),
        plot.margin = margin(t = 1, r = 5, b = 1, l = 1, unit = "pt"))




