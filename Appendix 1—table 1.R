
# ======================= Load packages & fonts ==============================
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

font_add("arial", 'arial.ttf', bold = 'arialbd.ttf',
         italic = 'ariali.ttf', bolditalic = 'arialbi.ttf')
font_add("times", 'times.ttf', bold = 'timesbd.ttf',
         italic = 'timesi.ttf', bolditalic = 'timesbi.ttf')
showtext_auto(enable = TRUE)

line_size <- 1/.pt
text_pt   <- 8
dpi       <- 300
text_size <- dpi/72 * text_pt
family    <- "times"

# Pick a suitable tick interval from a predefined range
find_interval <- function(min, max,
                          n_range = c(0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,
                                      0.5,1,2,3,4,5,10,20,50,100,200,300,400,
                                      500,600)) {
  ideal_step <- (max - min) / 6.5
  best_step  <- min(n_range[which(n_range >= ideal_step)])
  return(best_step)
}

# ======================= Project paths (relative) ===========================

project_root <- "./PutianInfluenza"   
work_dir     <- file.path(project_root, "WorkData")      
raw_dir      <- file.path(project_root, "RawData")       
weather_dir  <- file.path(raw_dir, "Putian_Weather")     

setwd(work_dir)

# ======================= Import data ========================================
rm(list = ls()); gc()
options(scipen = 100)
options(digits = 10)

save_path <- paste0(getwd(), "/")                       # output directory = work folder
if (!dir.exists(save_path)) { dir.create(save_path, recursive = TRUE) }

library(openxlsx)
WB <- createWorkbook()
library(lubridate)

# --- Weather data: read all Putian weather CSV files and bind rows ---
wt <- data.frame()
for (i in list.files(path = weather_dir, pattern = "Putian")) {
  wt <- rbind(wt, read.csv(file.path(weather_dir, i), header = TRUE))
}; rm(i)

wt$datetime  <- as.Date(gsub("-", "/", wt$datetime))
wt$tempdiff  <- wt$tempmax - wt$tempmin
wt$year_day  <- yday(wt$datetime)
wt <- wt[, c(2:5, 10:11, 18:19, 23, 25, 34, 35)]

# --- Influenza surveillance data ---
data1 <- read.csv(file.path(work_dir, "Putian_Flu_Surveillance_2018-2023-UTF8.csv"),
                  header = TRUE, fileEncoding = "UTF-8")
data2 <- read.csv(file.path(work_dir, "Putian_Flu_Surveillance_2023-UTF8.csv"),
                  header = TRUE, fileEncoding = "UTF-8")
cases <- unique(rbind(data1, data2[, -1])); rm(data1); rm(data2)

cases$report_date <- as.Date(cases$report_date)
cases <- subset(cases, cases$report_date >= as.Date("2018-01-01"))
cases <- cases[which(cases$pcr_result != ""), ]
cases <- cases[which(cases$pcr_result != "Mixed"), ]

# Map raw PCR results into broad type A / type B groups
cases$flutype <- cases$pcr_result
cases$flutype[which(cases$flutype == "A_H1N1" | cases$flutype == "A_H3N2")] <- "typeA"
cases$flutype[which(cases$flutype == "Victoria" | cases$flutype == "Yamagata")] <- "typeB"

# --- Daily case counts (overall / by broad type / by subtype) ---
casesili <- as.data.frame.array(table(cases$report_date))
casesili$Date <- as.Date(row.names(casesili))

casesAB <- as.data.frame.array(table(cases$report_date, cases$flutype))
casesAB$Date <- as.Date(row.names(casesAB))

casesflu <- as.data.frame.array(table(cases$report_date, cases$pcr_result))
casesflu$Date <- as.Date(row.names(casesflu))

cases_sum <- merge(casesili, merge(casesAB, casesflu, by = "Date"), by = "Date")
rm(casesili); rm(casesAB); rm(casesflu)

cases_sum$positive <- cases_sum$typeA + cases_sum$typeB
cases_sum <- cases_sum[, c("Date", "table(cases$report_date)", "positive",
                           "typeA", "typeB", "A_H1N1", "A_H3N2",
                           "Victoria", "Yamagata", "negative.x")]
colnames(cases_sum) <- c("Date", "ili", "positive", "typeA", "typeB",
                         "H1", "H3", "Victoria", "Yamagata", "negative")

# --- Merge daily cases with weather (keep all weather dates) ---
wtcases <- merge(cases_sum, wt, by.x = "Date", by.y = "datetime", all.y = TRUE)
for (i in 2:10) {
  wtcases[, i][is.na(wtcases[, i])] <- 0
  wtcases[, i] <- as.numeric(wtcases[, i])
}; rm(i)

wtcases <- subset(wtcases, wtcases$Date >= as.Date("2018-01-01"))

# --- Convert age to full years ---
cases[["age_years"]] <- cases$age
cases[["age_years"]][which(cases$age_unit == "month" & cases$age < 12)]  <- 0
cases[["age_years"]][which(cases$age_unit == "month" & cases$age > 11 & cases$age < 24)] <- 1
cases[["age_years"]][which(cases$age_unit == "month" & cases$age > 23)] <- 2
i <- 1
rm(i)

# --- Write daily influenza + weather summary sheet ---
addWorksheet(WB, sheetName = "DailyFluWeatherSummary")
writeData(WB, sheet = "DailyFluWeatherSummary", wtcases, rowNames = FALSE)

# ======================= Descriptive statistics table =======================
cols <- c("Var", "Min", "5%", "1st Qu", "Median", "Mean", "3rd Qu", "95%", "Max", "Sd")
summary_dt <- as.data.frame(matrix(nrow = ncol(wtcases) - 2, ncol = length(cols)))
colnames(summary_dt) <- cols; rm(cols)
summary_dt$Var <- colnames(wtcases[, c(2:(ncol(wtcases) - 1))])

for (i in c(2:(ncol(wtcases) - 1))) {
  summary_dt$Min[which(summary_dt$Var == colnames(wtcases)[i])]   <- min(wtcases[, i], na.rm = TRUE)
  summary_dt$`5%`[which(summary_dt$Var == colnames(wtcases)[i])]  <- round(quantile(wtcases[, i], 0.05, na.rm = TRUE), 3)
  summary_dt$`1st Qu`[which(summary_dt$Var == colnames(wtcases)[i])] <- round(quantile(wtcases[, i], 0.25, na.rm = TRUE), 3)
  summary_dt$Median[which(summary_dt$Var == colnames(wtcases)[i])] <- round(median(wtcases[, i], na.rm = TRUE), 3)
  summary_dt$Mean[which(summary_dt$Var == colnames(wtcases)[i])]  <- round(mean(wtcases[, i], na.rm = TRUE), 3)
  summary_dt$`3rd Qu`[which(summary_dt$Var == colnames(wtcases)[i])] <- round(quantile(wtcases[, i], 0.75, na.rm = TRUE), 3)
  summary_dt$`95%`[which(summary_dt$Var == colnames(wtcases)[i])] <- round(quantile(wtcases[, i], 0.95, na.rm = TRUE), 3)
  summary_dt$Max[which(summary_dt$Var == colnames(wtcases)[i])]   <- max(wtcases[, i], na.rm = TRUE)
  summary_dt$Sd[which(summary_dt$Var == colnames(wtcases)[i])]    <- round(sd(wtcases[, i], na.rm = TRUE), 3)
}; rm(i)

addWorksheet(WB, sheetName = "SummaryStatistics")
writeData(WB, sheet = "SummaryStatistics", summary_dt, rowNames = FALSE)
rm(summary_dt)
