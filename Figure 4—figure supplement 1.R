
project_root <- "./PutianInfluenza"                           
save_path    <- file.path(project_root, "Plots")              
raw_dir      <- file.path(project_root, "RawData")            
weather_dir  <- file.path(raw_dir, "Putian_Weather")          

rm(list = ls()); gc()
options(scipen = 100)
options(digits = 10)

if (!dir.exists(save_path)) { dir.create(save_path, recursive = TRUE) }
library(openxlsx)
WB <- createWorkbook()
library(lubridate)

# --- Weather data: read all Putian weather CSV files and bind rows ---
wt <- data.frame()
for (i in list.files(path = weather_dir, pattern = "Putian")) {
  wt <- rbind(wt, read.csv(file.path(weather_dir, i), header = TRUE))
}; rm(i)

wt$datetime <- as.Date(gsub("-", "/", wt$datetime))
wt$tempdiff <- wt$tempmax - wt$tempmin
wt$year_day <- yday(wt$datetime)
wt <- wt[, c(2:5, 10:11, 18:19, 23, 25, 34, 35)]

# --- Influenza surveillance data ---
data1 <- read.csv(file.path(raw_dir, "Putian_Flu_Surveillance_2018-2023-UTF8.csv"),
                  header = TRUE, fileEncoding = "UTF-8")
data2 <- read.csv(file.path(raw_dir, "Putian_Flu_Surveillance_2023-UTF8.csv"),
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
cases[["age_years"]][which(cases$age_unit == "month" & cases$age < 12)] <- 0
cases[["age_years"]][which(cases$age_unit == "month" & cases$age > 11 & cases$age < 24)] <- 1
cases[["age_years"]][which(cases$age_unit == "month" & cases$age > 23)] <- 2
i <- 1
rm(i)

# --- Write daily influenza + weather summary sheet ---
addWorksheet(WB, sheetName = "DailyFluWeatherSummary")
writeData(WB, sheet = "DailyFluWeatherSummary", wtcases, rowNames = FALSE)

# ======================= Fig S. Correlation matrix ==========================
cor_fn <- function(data, mapping, method = "p", use = "pairwise", ...) {
  x    <- eval_data_col(data, mapping$x)
  y    <- eval_data_col(data, mapping$y)
  corr <- cor(x, y, method = method, use = use)
  colFn <- colorRampPalette(c("darkgreen", "white", "darkred"))
  fill  <- colFn(100)[findInterval(corr, seq(-1, 1, length = 100))]
  ggally_cor(data = data, mapping = mapping, title = "Cor", colour = "black",
             family = family, size = (dpi/72 * 6)/.pt, lineheight = 0.3) +
    theme(panel.background = element_rect(fill = fill, colour = NA),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black"))
}

point_fc <- function(data, mapping, ...) {
  ggally_points(data = data, mapping = mapping, color = "#377EB8",
                alpha = 0.2, size = 2/.pt) +
    theme_bw() +
    theme(panel.grid = element_blank())
}

den_fc <- function(data, mapping, ...) {
  ggally_densityDiag(data = data, mapping = mapping,
                     color = "#003366", size = line_size) +
    theme_bw() +
    theme(panel.grid = element_blank())
}

p_cor <- ggpairs(wtcases[, c(2, 4, 5, 13, 20, 14, 15, 16, 17, 18, 19)],
                 upper = list(continuous = cor_fn),
                 diag = list(continuous = den_fc),
                 lower = list(continuous = point_fc),
                 columnLabels = c("Cases\nof ILI",
                                  "influenza\nA cases", "influenza\nB cases",
                                  "Avg.Temp\n(°C)", "Temp\nRange\n(°C)",
                                  "Humidity\n(%)", "Precipi-\ntation\n(mm)",
                                  "Wind\nspeed\n(km/h)", "Wind\ndirection",
                                  "Solar\nradiation\n(W/m2)", "UV index"),
                 switch = "both",
                 progress = FALSE)

p_cor + theme(
  strip.text = element_text(size = (dpi/72 * 6), family = family,
                            colour = "black", hjust = 0.5, lineheight = 0.25),
  axis.text.x = element_text(size = (dpi/72 * 6), family = family,
                             colour = "black", angle = 90, vjust = 0.5),
  axis.text.y = element_text(size = (dpi/72 * 6), family = family,
                             colour = "black"),
  plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt"))

ggsave(paste0(save_path, "Correlation_Matrix.tiff"),
       height = 7, width = 7, units = "in", bg = NULL, dpi = dpi)
rm(cor_fn, point_fc, den_fc, p_cor)










