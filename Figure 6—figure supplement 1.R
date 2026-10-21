###################

library(readr)
library(openxlsx)
library(forecast)
library(ggplot2)
library(patchwork)
if(!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
setwd("./")
wtcases <- read.xlsx("df_merge3_joined.xlsx")
save_path <- "./ARIMA_results(multivariate)/20260701_update/"


colnames(wtcases)
str(wtcases[, c("temp", "tempdiff", "humidity", "precip",  "windspeed", "uvindex", "winddir", "solarradiation","DOW", "weekly_detection","facial_coverings","outsiders_proportion")])
wtcases$outsiders_proportion <- as.numeric(wtcases$outsiders_proportion)

wtcases$Date <- as.Date(wtcases$Date, origin = "1899-12-30")
wtcases$typeA_rate <- as.numeric(wtcases$typeA_rate)
wtcases$typeB_rate <- as.numeric(wtcases$typeB_rate)

arima_fc_fixed <- function(var) {
  data <- wtcases
  data$Date <- as.Date(data$Date)
  data$outsiders_proportion <- as.numeric(data$outsiders_proportion)
  data[[var]] <- as.numeric(data[[var]])
  

  xreg_vars <- c("temp", "tempdiff", "humidity", "precip", "windspeed",
                 "uvindex", "winddir", "solarradiation", 
                 "DOW", "weekly_detection", "facial_coverings", "outsiders_proportion")
  
  train_idx <- which(data$Date < "2023-01-01")
  test_idx <- which(data$Date >= "2023-01-01")
  

  train <- ts(data[[var]][train_idx], frequency = 7)
  test <- ts(data[[var]][test_idx], frequency = 7)
  
  X_train <- as.matrix(data[train_idx, xreg_vars])
  X_test <- as.matrix(data[test_idx, xreg_vars])
  

  fit_arima <- auto.arima(train, xreg = X_train, trace = TRUE)
  
  fc <- forecast(fit_arima, xreg = X_test, h = length(test))
  

  pred_vals <- pmax(0, pmin(1, as.numeric(fc$mean)))
  
  results <- data.frame(actual = as.numeric(test), predicted = pred_vals)
  
  actual <- results$actual
  predicted <- results$predicted
  valid_idx <- actual != 0
  
  MAE <- mean(abs(actual - predicted))
  RMSE <- sqrt(mean((actual - predicted)^2))
  MAPE <- mean(abs((actual[valid_idx] - predicted[valid_idx]) / actual[valid_idx]))
  SMAPE <- mean(2 * abs(actual[valid_idx] - predicted[valid_idx]) /
                  (abs(actual[valid_idx]) + abs(predicted[valid_idx])))
  
  write.csv(results, paste0(save_path, "ARIMAX_Forecasting_", var, "_fixed.csv"), row.names = FALSE)
  
  cat("=====", var, "(fixed) =====\n")
  cat("MAE:", MAE, "\nRMSE:", RMSE, "\nMAPE:", MAPE, "\nSMAPE:", SMAPE, "\n\n")
  
  return(results)
}


arima_a <- arima_fc_fixed("typeA_rate")
arima_b <- arima_fc_fixed("typeB_rate")


#####

train_idx <- which(wtcases$Date < "2023-01-01")
test_idx <- which(wtcases$Date >= "2023-01-01")

arima_a$Date <- wtcases$Date[test_idx]
arima_b$Date <- wtcases$Date[test_idx]





arima_a$Date <- wtcases$Date[test_idx]
arima_a$Observations <- wtcases$PositiveA[test_idx] # actual → Observations
arima_a$Predictions <- arima_a$Point.Forecast # predicted → Predictions

arima_b$Date <- wtcases$Date[test_idx]
arima_b$Observations <- wtcases$PositiveB[test_idx]
arima_b$Predictions <- arima_b$Point.Forecast

names(arima_a) <- c("Observations", "Predictions", "Date")
names(arima_b) <- c("Observations", "Predictions", "Date")





month_abbr <- c("Jan","Feb","Mar","Apr","May","Jun",
                "Jul","Aug","Sep","Oct","Nov","Dec")


plot_arima <- function(df, laby = "Positive rate") {
  ggplot(df, aes(x = Date)) +
    geom_line(aes(y = Observations, color = "Observations"), linewidth = 0.8) +
    geom_line(aes(y = Predictions, color = "Predictions"), linewidth = 0.8, linetype = "dashed") +
    scale_color_manual(
      values = c("Observations" = "#E15759", "Predictions" = "#377EB8"),
      name = NULL
    ) +
    scale_x_date(
      labels = function(x) month_abbr[as.numeric(format(x, "%m"))],
      date_breaks = "2 months"
    ) +
    labs(x = "Date", y = laby) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    )
}



p <- (plot_arima(arima_a, laby = "Positivity rate of influenza A virus") +
        plot_arima(arima_b, laby = "Positivity rate of influenza B virus")) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A") &
  theme(
    panel.spacing = unit(5, "mm"),
    legend.position = "top",
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(t = -5, b = 10),
    plot.margin = margin(t = 15, r = 10, b = 10, l = 10),
    plot.tag = element_text(size = 12, colour = "black")
  )

print(p)




