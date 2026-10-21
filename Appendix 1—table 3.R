
##ARIMA

library(readr)
library(openxlsx)
library(forecast)
library(ggplot2)
library(patchwork)

setwd("./1.Putian Influenza Project/operation_files")

wtcases <- read.xlsx("data_with_positive_rate_migrant_ratio/df_merge3_joined.xlsx")

save_path <- "./ARIMA_results_multivariate/20260701_update/"

if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

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

##########################XGB

library(xgboost)
library(readr)
library(openxlsx)
library(ggplot2)
library(patchwork)

setwd("./1.Putian_Influenza_Project/operation_files")

# read data
wtcases <- read.xlsx("data_positive_rate_migrant_ratio/df_merge3_joined.xlsx")


wtcases$Date <- as.Date(wtcases$Date, origin = "1899-12-30")
wtcases$outsiders_proportion <- as.numeric(wtcases$outsiders_proportion)
wtcases$typeA_rate <- as.numeric(wtcases$typeA_rate)
wtcases$typeB_rate <- as.numeric(wtcases$typeB_rate)


param_grid <- expand.grid(
  max_depth = c(3,4,5,6),
  eta       = c(0.01, 0.05, 0.1),
  lambda    = c(0.5, 1, 2)
)


xgb_fc_grid <- function(var){
  
  data <- wtcases
  data <- na.omit(data)
  
  feature_cols <- c("temp", "tempdiff", "humidity", "precip", "windspeed",
                    "uvindex", "winddir",
                    "DOW", "weekly_detection", "facial_coverings",
                    "outsiders_proportion")
  
  train_all_idx <- which(data$Date < "2023-01-01")
  test_idx      <- which(data$Date >= "2023-01-01")
  
  df_train_all <- data[train_all_idx, ]
  df_test      <- data[test_idx, ]
  
  X_train_all <- as.matrix(df_train_all[, feature_cols])
  y_train_all <- df_train_all[[var]]
  
  X_test <- as.matrix(df_test[, feature_cols])
  y_test <- df_test[[var]]
  
 
  ts_cv_mae <- function(params_row){
    md  <- params_row$max_depth
    et  <- params_row$eta
    lam <- params_row$lambda
    
    n <- nrow(df_train_all)
   
    init_train_size <- floor(n * 0.6)
    step <- 12
    cv_mae_vec <- c()
    
    for(i in seq(from = init_train_size, to = n - step, by = step)){
      tr_idx <- 1:i
      val_idx <- (i+1):(i+step)
      
      dt <- xgb.DMatrix(X_train_all[tr_idx, ], label = y_train_all[tr_idx])
      dv <- xgb.DMatrix(X_train_all[val_idx, ], label = y_train_all[val_idx])
      
      p <- xgb.train(
        list(objective = "reg:squarederror", eval_metric = "mae",
             max_depth = md, eta = et, lambda = lam,
             subsample = 0.8, colsample_bytree = 0.8, alpha = 0, seed = 42),
        dt, nrounds = 800, watchlist = list(val=dv),
        early_stopping_rounds = 30, print_every_n = 0, verbose = 0
      )
      pred_val <- predict(p, dv)
      mae_fold <- mean(abs(pred_val - y_train_all[val_idx]))
      cv_mae_vec <- c(cv_mae_vec, mae_fold)
    }
    return(mean(cv_mae_vec))
  }
  
  
  param_grid$cv_mae <- apply(param_grid, 1, function(r){
    ts_cv_mae(as.data.frame(t(r)))
  })
  
  
  best_row <- param_grid[which.min(param_grid$cv_mae), ]
  cat("\n==== ", var, " Optimal grid parameters ====\n")
  print(best_row)
  
  
  write.csv(param_grid, paste0(save_path,"GridSearch_",var,".csv"), row.names = FALSE)
  
  # ===
  dtrain_final <- xgb.DMatrix(X_train_all, label = y_train_all)
  dtest_final  <- xgb.DMatrix(X_test, label = y_test)
  
  best_params <- list(
    objective = "reg:squarederror",
    eval_metric = "mae",
    max_depth = best_row$max_depth,
    eta       = best_row$eta,
    lambda    = best_row$lambda,
    subsample = 0.8,
    colsample_bytree = 0.8,
    alpha = 0,
    seed = 42
  )
  
  fit_final <- xgb.train(
    best_params,
    dtrain_final,
    nrounds = 1000,
    watchlist = list(train = dtrain_final),
    print_every_n = 100
  )
  
  pred <- pmax(0, pmin(1, predict(fit_final, dtest_final)))
  
  results <- data.frame(
    actual = y_test,
    predicted = pred,
    Date = df_test$Date
  )
  
  valid_idx <- results$actual != 0
  MAE  <- mean(abs(results$actual - results$predicted))
  RMSE <- sqrt(mean((results$actual - results$predicted)^2))
  MAPE <- mean(abs((results$actual[valid_idx] - results$predicted[valid_idx]) / results$actual[valid_idx]))
  SMAPE<- mean(2 * abs(results$actual[valid_idx] - results$predicted[valid_idx]) /
                 (abs(results$actual[valid_idx]) + abs(results$predicted[valid_idx])))
  
  write.csv(results, paste0(save_path, "XGB_Forecasting_", var, ".csv"), row.names = FALSE)
  
  imp <- xgb.importance(feature_names = feature_cols, model = fit_final)
  xgb.plot.importance(imp)
  
  cat("\n===== ", var, "(XGBoost grid‑search final model) =====\n")
  cat("MAE:", MAE, "\nRMSE:", RMSE, "\nMAPE:", MAPE, "\nSMAPE:", SMAPE, "\n\n")
  
  return(list(results = results, best_param = best_row, importance = imp))
}


res_a <- xgb_fc_grid("typeA_rate")
res_b <- xgb_fc_grid("typeB_rate")

xgb_a <- res_a$results
xgb_b <- res_b$results


#########Bayesian-optimised LSTM
#The MAE, RMSE, MAPE and SMAPE indicators of the Bayesian-optimised LSTM need to be calculated using the program in Figure6.py.




