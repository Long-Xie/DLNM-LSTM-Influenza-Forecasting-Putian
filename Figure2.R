

# install.packages("ash")
# install.packages("maps")
# install.packages("proj4")

# Plot settings
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

# Paths relative to project root
raw_dir    <- "./data/weather"
data_dir   <- "./data"
output_dir <- "./output"
if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}

# Import data
rm(list=ls()); gc()
options(scipen = 100)
options(digits = 10)

library(openxlsx)
WB <- createWorkbook()
library(lubridate)

# Import weather data
wt <- data.frame()
weather_files <- list.files(path = raw_dir, pattern = "Putian")
for (i in weather_files) {
  wt <- rbind(wt, read.csv(file.path(raw_dir, i), header = TRUE))   
}; rm(i)

wt$datetime <- as.Date(gsub("-","/",wt$datetime))
wt$tempdiff <- wt$tempmax - wt$tempmin
wt$year_day <- yday(wt$datetime)
wt <- wt[,c(2:5,10:11,18:19,23,25,34,35)]

# Import flu surveillance data
data1 <- read.csv(file.path(data_dir, "2018-2023_Putian_influenza_surveillance-UTF8.csv"),
                  header = TRUE, fileEncoding = "UTF-8")
data2 <- read.csv(file.path(data_dir, "2023_Putian_influenza_surveillance-UTF8.csv"),
                  header = TRUE, fileEncoding = "UTF-8")
cases <- unique(rbind(data1, data2[,-1])); rm(data1); rm(data2)

cases$onset_date <- as.Date(cases$onset_date)
cases <- subset(cases, cases$onset_date >= as.Date("2018-01-01"))
cases <- cases[which(cases$pcr_result != ""),]
cases <- cases[which(cases$pcr_result != "mixed"),]

cases$flutype <- cases$pcr_result
cases$flutype[which(cases$flutype=="H1"|cases$flutype=="H3")] <- "typeA"
cases$flutype[which(cases$flutype=="Victoria"|cases$flutype=="Yamagata")] <- "typeB"

casesili <- as.data.frame.array(table(cases$onset_date))
casesili$Date <- as.Date(row.names(casesili))
casesAB <- as.data.frame.array(table(cases$onset_date, cases$flutype))
casesAB$Date <- as.Date(row.names(casesAB))
casesflu <- as.data.frame.array(table(cases$onset_date, cases$pcr_result))
casesflu$Date <- as.Date(row.names(casesflu))

cases_sum <- merge(casesili,merge(casesAB, casesflu, by = "Date"), by = "Date")
rm(casesili);rm(casesAB);rm(casesflu)
cases_sum$positive <- cases_sum$typeA+cases_sum$typeB
cases_sum <- cases_sum[,c("Date","table(cases$onset_date)","positive","typeA","typeB","H1","H3","Victoria","Yamagata","negative.x")]
colnames(cases_sum) <- c("Date","ili","positive","typeA","typeB","H1","H3","Victoria","Yamagata","negative")

wtcases <- merge(cases_sum, wt, by.x = "Date", by.y = "datetime", all.y = TRUE)
for (i in 2:10) { 
  wtcases[,i][is.na(wtcases[,i])] <-0
  wtcases[,i] <- as.numeric(wtcases[,i]) 
}; rm(i)
wtcases <- subset(wtcases, wtcases$Date >= as.Date("2018-01-01"))

cases[["age_year"]]<- cases$age
cases[["age_year"]][which(cases["age_unit"]=="month" & cases$age < 12)] <- 0 
cases[["age_year"]][which(cases["age_unit"]=="month" & cases$age > 11 & cases$age < 24)] <- 1  
cases[["age_year"]][which(cases["age_unit"]=="month" & cases$age > 23)] <- 2

addWorksheet(WB, sheetName = "Daily_Influenza_Meteorology")
writeData(WB, sheet = "Daily_Influenza_Meteorology", wtcases, rowNames = FALSE)

# Descriptive statistics table
cols <- c("Var", "Min","5%","1st Qu","Median","Mean","3rd Qu","95%","Max","Sd")
summary_dt <- as.data.frame(matrix(nrow = ncol(wtcases)-2, ncol = length(cols)))
colnames(summary_dt) <- cols; rm(cols)
summary_dt$Var <- colnames(wtcases[,c(2:(ncol(wtcases)-1))])

for (i in c(2:(ncol(wtcases)-1))) {
  summary_dt$Min[which(summary_dt$Var==colnames(wtcases)[i])] <- min(wtcases[,i],na.rm = TRUE)
  summary_dt$`5%`[which(summary_dt$Var==colnames(wtcases)[i])] <- round(quantile(wtcases[,i],0.05,na.rm = TRUE),3)
  summary_dt$`1st Qu`[which(summary_dt$Var==colnames(wtcases)[i])] <- round(quantile(wtcases[,i],0.25,na.rm = TRUE),3)
  summary_dt$Median[which(summary_dt$Var==colnames(wtcases)[i])] <- round(median(wtcases[,i],na.rm = TRUE),3)
  summary_dt$Mean[which(summary_dt$Var==colnames(wtcases)[i])] <- round(mean(wtcases[,i],na.rm = TRUE),3)
  summary_dt$`3rd Qu`[which(summary_dt$Var==colnames(wtcases)[i])] <- round(quantile(wtcases[,i],0.75,na.rm = TRUE),3)
  summary_dt$`95%`[which(summary_dt$Var==colnames(wtcases)[i])] <- round(quantile(wtcases[,i],0.95,na.rm = TRUE),3)
  summary_dt$Max[which(summary_dt$Var==colnames(wtcases)[i])] <- max(wtcases[,i],na.rm = TRUE)
  summary_dt$Sd[which(summary_dt$Var==colnames(wtcases)[i])] <- round(sd(wtcases[,i],na.rm = TRUE),3)
}; rm(i)

addWorksheet(WB, sheetName = "Descriptive_Statistics")  
writeData(WB, sheet = "Descriptive_Statistics", summary_dt, rowNames = FALSE)
rm(summary_dt)

# Save workbook
saveWorkbook(WB, file.path(output_dir, "summarized_data.xlsx"), overwrite = TRUE)







###############figure2
#####################Fig2A
# Data file: ./data/Fig2A.csv

# Load package
library(ggplot2)

# Read data with relative path
df <- read.csv("./data/Fig2A.csv")

df$months <- factor(df$months, 
                    levels = c("Jan","Feb","Mar","Apr","May","Jun",
                               "Jul","Aug","Sep","Oct","Nov","Dec"),
                    ordered = TRUE)



ggplot(df, aes(x = months, y = count, color = type, group = type)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  facet_wrap(~ year, nrow = 2, ncol = 3) +
  scale_color_manual(values = my_colors) +
  theme_bw(base_size = 15) +
  theme(
    text = element_text(family = "Times New Roman"),  
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  ) +
  labs(
    title = "",
    x = "Month",
    y = "Cases",
    color = "Type"
  )


#######
colors <- c( "#F1DFA4" , "#74B69F",  "#E2C8D8", "#A6CDE4" )
cases_sum$Date <- as.Date(cases_sum$Date, origin = "1899-12-30")
cases_sum$years <- year(cases_sum$Date)
cases_sum$months <- month(cases_sum$Date)
month_hp <- cases_sum %>%
  group_by(years, months) %>%
  dplyr::summarise(H1 = sum(H1),  H3 = sum(H3), Victoria = sum(Victoria),  Yamagata = sum(Yamagata))
month_hp <- gather(month_hp, key = "type", value = "count", -years, -months)
for (i in 1:nrow(month_hp)) {
  month_hp$years[i] <- paste0(month_hp$years[i], month_hp$type[i]) }; rm(i)
month_hp$months <- factor(month_hp$months, levels = c("12","11","10","9","8","7","6","5","4","3","2","1"),
                          labels = c("Dec","Nov","Oct","Sep","Aug","Jul","Jun","May","Apr","Mar","Feb","Jan"))
month_hp$type <- factor(month_hp$type, levels = c("H1", "H3", "Victoria","Yamagata"),
                        labels = c("A / H1N1","A / H3N2","B / Victoria","B / Yamagata"))
flu <- cases[which(cases$PCR_result != "negative"),]
flu$age_group <- cut(flu$age_years, breaks = seq(0,90,10), include.lowest = T,right = FALSE,ordered_result = TRUE)
age <- flu[,c("gender", "age_group", "PCR_result")]
age <- age %>%
  group_by(gender, age_group,  PCR_result) %>%
  dplyr::summarise(case_count = n())
age[which(age$gender=="male"),4] <- 0-age[which(age$gender=="male"),4]
age$PCR_result <- factor(age$PCR_result, levels = c("H1", "H3", "Victoria","Yamagata"),
                         labels = c("A / H1N1","A / H3N2","B / Victoria","B / Yamagata"))
group_plot <- ggplot(month_hp,aes(x=years,y=months))+
  geom_tile(aes(fill=type))+
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  geom_vline(xintercept=c(4.5, 8.5, 12.5, 16.5, 20.5),size=line_size, color = "black") +
  scale_fill_manual(values = colors, name="" )+
  theme(panel.border = element_rect(fill=NA, color="black", size=line_size,linetype="solid"),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "right",
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.key.size  = unit(0, 'cm'),
        legend.key.spacing.y = unit(0, 'cm'),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
hp <- ggplot(month_hp,aes(x=years,y=months,fill=count)) + geom_raster() + 
  scale_fill_gradient(low="white", high="#912C2C" ,limits = c(0,300), breaks = seq(0,300,100)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0), labels = c("","","2018","","","","2019","","","","2020","","","","2021","","","","2022","","","","2023","")) +
  geom_vline(xintercept=c(4.5, 8.5, 12.5, 16.5, 20.5), size=line_size, color = "black") +
  coord_fixed(ratio = 1) +
  labs(x="Years", y="Month", fill="Cases") +
  theme(panel.border = element_rect(fill=NA, color="black", size=line_size, linetype="solid"),
        axis.ticks.y = element_line(color="black", size=line_size),
        axis.ticks.x = element_blank(),
        text = element_text(size = text_size, family=family, colour = "black"),
        axis.text.y = element_text(size = text_size, family=family, colour = "black"),
        axis.text.x = element_text(hjust = 1,size = text_size, family=family, colour = "black"),
        legend.key.width = unit(0.3, 'cm'),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt") )
hp


###########################Fig2B

fig2b <- ggplot(age,aes(x=age_group,y=case_count,fill = PCR_result)) +
  geom_bar(stat="identity", width = 1.8, position = position_dodge(preserve = 'single')) + 
  coord_flip() + geom_hline(yintercept = 0, size = line_size, color="black") +
  scale_y_continuous(expand = c(0,0), breaks = seq(-900, 900, 100), labels = as.character(abs(seq(-900, 900, 100))), limits = c(-400, 400)) +
  scale_x_discrete(expand = c(0,0)) +
  labs(x="Age, years", y="Influenza cases") +
  scale_fill_manual(values = colors, name="Type" )+
  theme_bw() + theme_classic() +
  theme(legend.key.size  = unit(0.3, 'cm'),
        legend.key.spacing.y = unit(0.1, 'cm'),
        legend.direction = "vertical",
        legend.position = "right",
        legend.text = element_text(size = text_size, family=family, colour = "black"),
        legend.title = element_text(size = text_size, family=family, colour = "black"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        axis.text = element_text(size = text_size, family=family, colour = "black"),
        axis.title = element_text(size = text_size, family=family, colour = "black"),
        plot.margin = margin(t = 5, r = 10, b = 1, l = 1, unit = "pt")        ) +
  annotate('text',label = 'Female', x = 8.5, y = 200, size = text_size/.pt, family=family, colour = "black") +
  annotate('text',label = 'Male', x = 8.5, y = -200, size = text_size/.pt, family=family, colour = "black")
fig2b

ggsave(paste0(save_path,"Fig2B.tiff"), height = 3.5, width = 7,units = "in",bg = NULL, dpi=dpi)

group_plot +theme(legend.key.size  = unit(0.3, 'cm'),
                  legend.key.spacing.y = unit(0.1, 'cm'),
                  legend.direction = "vertical",
                  legend.position = "right",
                  legend.text = element_text(size = text_size, family=family, colour = "black"),
                  legend.title = element_text(size = text_size, family=family, colour = "black"),
                  legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) + 
  hp + 
  plot_layout(ncol = 1, heights = c(0.05, 0.95), guides = "collect") 

ggsave(paste0(save_path,"Fig2A.tiff"), height = 3.5, width = 7,units = "in",bg = NULL, dpi=dpi)

patchwork <- ((group_plot/hp) ) / fig2b
patchwork[[1]] <- patchwork[[1]] + plot_layout(tag_level = 'new')
patchwork[[2]] <- patchwork[[2]] + plot_layout(tag_level = 'new')

patchwork + 
  plot_layout(ncol = 1, heights = c(0.025,0.45,0.5), guides = 'collect') +
  plot_annotation(tag_levels = list(c('A', '',"B"), ''))&
  theme(panel.spacing.y = unit(0, "cm"),
        legend.position = "right",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.spacing.y = unit(c(-3,2), "cm"),
        plot.tag = element_text(size = text_size, colour = "black", family = family),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt"))

ggsave(paste0(save_path,"Fig2_upper.tiff"), height = 7, width = 7,units = "in",bg = NULL, dpi=dpi)

patchwork + 
  plot_layout(ncol = 1, heights = c(0.025,0.45,0.5), guides = 'collect') +
  plot_annotation(tag_levels = list(c('a', '',"b"), ''))&
  theme(panel.spacing.y = unit(0, "cm"),
        legend.position = "right",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.spacing.y = unit(c(-3,2), "cm"),
        plot.tag = element_text(size = text_size, colour = "black", family = family),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt"))

ggsave(paste0(save_path,"Fig2_lower.tiff"), height = 7, width = 7,units = "in",bg = NULL, dpi=dpi)

addWorksheet(WB, sheetName = "Fig2A")  
writeData(WB, sheet = "Fig2A", month_hp, rowNames = FALSE)
addWorksheet(WB, sheetName = "Fig2B")  
writeData(WB, sheet = "Fig2B", age, rowNames = FALSE)
saveWorkbook(WB, paste0(save_path, "processed_data.xlsx"), overwrite = TRUE, returnValue = FALSE)
rm(WB)
rm(fig2b,group_plot, hp, month_hp,patchwork,flu,age,colors)




