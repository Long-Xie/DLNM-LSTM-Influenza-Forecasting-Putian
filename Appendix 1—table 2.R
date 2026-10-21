
save_path <-".\\Daangene\\1.Putian Influenza Project\\"
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

dir.create(paste0(save_path,"DLNM_Plots/RawData/df_selection/"), recursive = TRUE)
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
  return(qaic)    }
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
                                               "-","-","-","-")  } } 
  } else {   picture_all[nrow(picture_all)+1,] <- c(respon,var,num,"NULL","NULL","-","-","-","-","-","-") }
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
                                                 arr.ind = TRUE)] ,3))  }  }  
  } else { picture_all[nrow(picture_all)+1,] <- c(respon,var,num,"-","-","-","-","NULL","NULL","-","-")    }  }






wtcases <- read.csv("wtcases_0226.csv", fileEncoding = "GBK")
##DLNM####
wtcases$winddir <- as.numeric(wtcases$winddir)
wtcases <- wtcases[order(wtcases$Date),]
wtcases$time <- c(1:nrow(wtcases)) 
wtcases$winddir_use <- wtcases$winddir
wtcases$winddir_use[which(wtcases$winddir>337.5)] <- wtcases$winddir[which(wtcases$winddir>337.5)]-360
flu_types <- c("typeA", "typeB")
vars_list <- c("temp", "tempdiff", "humidity", "precip",  "windspeed", "uvindex", "winddir_use", "solarradiation")
for (respon in flu_types) {
  for (var in vars_list) {
    M <- median(as.numeric(wtcases[,var]))
    lag <- 15
    qaic_model <- as.data.frame(matrix(nrow = 4, ncol = 5))
    for (i in 2:5){ 
      for (n in 6:10) {  
        try({cb_test = crossbasis(wtcases[,var],lag = lag,argvar =list(fun="ns",df=i), arglag = list(fun="poly", degree = 4))
        model_test <-  glm(wtcases[,respon] ~ cb_test + ns(time,n*5), family=quasipoisson(), data = wtcases, na.action = na.omit)
        qaic_model[i-1,n-5] <- as.numeric(paste(fqaic(model_test)))  }, silent = TRUE)
      }  }
    colnames(qaic_model) <- c(6:10)
    rownames(qaic_model) <- c(2:5)
    write.csv(qaic_model,file = paste0(save_path,"DLNM_Plots/RawData/df_selection/",respon,"-",var,"df_selection",Sys.Date(),".csv"))
    df_cb <- which(qaic_model == min(qaic_model, na.rm = TRUE), arr.ind = TRUE)[,1]+1
    df_time <- which(qaic_model == min(qaic_model, na.rm = TRUE), arr.ind = TRUE)[,2]+5
    cb <- crossbasis(wtcases[,var], lag = lag, argvar = list(fun="ns",df=df_cb), arglag = list(fun="poly", degree = 4))
    model <-  glm(wtcases[,respon] ~ cb + ns(time,df_time*5), family=quasipoisson(), data = wtcases, na.action = na.omit)
    coefficients <- as.data.frame(summary(model)$coefficients)
    list_p0.1 <- rownames(coefficients[which(coefficients$`Pr(>|t|)`<0.1),])
    model_all[nrow(model_all)+1,] <- c(respon, var, df_cb, df_time,
                                       round(fqaic(model),3),sum(str_detect(list_p0.1,"cb")),
                                       sum(str_detect(list_p0.1,"time")),"","")
    pred <- crosspred(cb, model, cen = M, bylag=0.2)
    rm(cb);rm(cb_test);rm(coefficients);rm(list_p0.1);rm(model);rm(model_test);rm(qaic_model);rm(i);rm(n);
    RRfit <- as.data.frame(pred$matRRfit)
    breaks_lines <- seq(min(RRfit),max(RRfit),by=((max(RRfit)-min(RRfit))/10)) 
    RRfit$var <- rownames(RRfit) 
    map <- melt(RRfit) 
    colnames(map)<-c("var","lag_day","value_")
    map$var <- as.numeric(map$var)
    map$lag_day <- as.numeric(gsub('[lag]',"",map$lag_day))
    RRall <- as.data.frame(pred$allRRfit)
    RRall$var <- as.numeric(rownames(RRall))
    colnames(RRall) <- c("RR","var")
    RRall.low <- as.data.frame(pred$allRRlow)
    RRall.low$var <- as.numeric(rownames(RRall.low))
    colnames(RRall.low) <- c("RR.low","var")
    RRall.high <- as.data.frame(pred$allRRhigh)
    RRall.high$var <- as.numeric(rownames(RRall.high))
    colnames(RRall.high) <- c("RR.high","var")
    comRR <- merge(merge(RRall,RRall.high,by = "var"),RRall.low,by = "var")
    rm(RRall);rm(RRall.low);rm(RRall.high)
    row_min = which.min(abs(quantile(wtcases[,var],0.05) - as.numeric(rownames(RRfit))))
    model_all[nrow(model_all),8] <- rownames(as.data.frame(pred$matRRfit)[row_min,])
    low.RR <-  as.data.frame(t(as.data.frame(pred$matRRfit)[row_min,]))
    low.RR[,2] <- as.numeric(substr(rownames(low.RR),4,8))
    colnames(low.RR) <- c("RR","lag")
    low.RR.low <-  as.data.frame(t(as.data.frame(pred$matRRlow)[row_min,]))
    low.RR.low[,2] <- as.numeric(substr(rownames(low.RR.low),4,8))
    colnames(low.RR.low) <- c("RR.low","lag")
    low.RR.up <-  as.data.frame(t(as.data.frame(pred$matRRhigh)[row_min,]))
    low.RR.up[,2] <- as.numeric(substr(rownames(low.RR.up),4,8))
    colnames(low.RR.up) <- c("RR.up","lag")
    low <- merge(merge(low.RR,low.RR.up,by="lag"),low.RR.low,by="lag")
    row_max = which.min(abs(quantile(wtcases[,var],0.95) - as.numeric(rownames(RRfit))))
    model_all[nrow(model_all),9] <- rownames(as.data.frame(pred$matRRfit)[row_max,])
    high.RR <-  as.data.frame(t(as.data.frame(pred$matRRfit)[row_max,]))
    high.RR[,2] <- as.numeric(substr(rownames(high.RR),4,8))
    colnames(high.RR) <- c("RR","lag")
    high.RR.low <-  as.data.frame(t(as.data.frame(pred$matRRlow)[row_max,]))
    high.RR.low[,2] <- as.numeric(substr(rownames(high.RR.low),4,8))
    colnames(high.RR.low) <- c("RR.low","lag")
    high.RR.up <-  as.data.frame(t(as.data.frame(pred$matRRhigh)[row_max,]))
    high.RR.up[,2] <- as.numeric(substr(rownames(high.RR.up),4,8))
    colnames(high.RR.up) <- c("RR.up","lag")
    high <- merge(merge(high.RR,high.RR.up,by="lag"),high.RR.low,by="lag")
    rm(low.RR,low.RR.up,low.RR.low,high.RR,high.RR.up,high.RR.low,df_cb,df_time,row_min,row_max)
    write.csv(comRR,file = paste0(save_path,"DLNM_Plots/RawData/",respon,"_",var,"_P2_",Sys.Date(),".csv"),row.names = FALSE)
    write.csv(low,file = paste0(save_path,"DLNM_Plots/RawData/",respon,"_",var,"_P3_",Sys.Date(),".csv"),row.names = FALSE)
    write.csv(high,file = paste0(save_path,"DLNM_Plots/RawData/",respon,"_",var,"_p4_",Sys.Date(),".csv"),row.names = FALSE)
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
    rm(comRR, low, high)  }  }
write_excel_csv(picture_all, file = paste0(save_path,"DLNM_Plots/RawData/DLNM_P2-P4_KeyPoints_Summary",Sys.Date(),".csv"))
write_excel_csv(model_all, file = paste0(save_path,"DLNM_Plots/RawData/DLNM_Model_Summary",Sys.Date(),".csv"))
rm(M, lag, respon, var, vars_list,flu_types, map, pred, RRfit, breaks_lines, fqaic, cut_fc, model_all, picture_all)
