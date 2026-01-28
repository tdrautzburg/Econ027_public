#' Multivariate Lag-Augmented LP with Optional Contemporaneous Controls
#' 
#' @param extra_controls A 'ts' object or list of 'ts' objects to be added as 
#'                       CONTEMPORANEOUS controls (at time t).
f_run_lag_augmented_LP <- function(lhs_vars_list, shock_ts, 
                                 extra_controls = NULL, # <--- NEW ARGUMENT
                                 h_max = 48,
                                 n_lag_vars = 24, 
                                 n_lag_shock = 36,
                                 sample_start = c(1970, 1), 
                                 sample_end = c(1996, 12)) {
  
  # 1. Construct Base Lagged Controls (Same as before)
  control_list <- list()
  
  # A. Shock Lags
  for(k in 1:n_lag_shock) control_list[[paste0("LagS_", k)]] <- stats::lag(shock_ts, -k)
  
  # B. Variable Lags
  var_names <- names(lhs_vars_list)
  for (v_name in var_names) {
    series <- lhs_vars_list[[v_name]]
    for(k in 1:n_lag_vars) {
      control_list[[paste0("Lag_", v_name, "_", k)]] <- stats::lag(series, -k)
    }
  }
  
  # 2. Add Extra Contemporaneous Controls (NEW STEP)
  # We do NOT lag these. They enter at time t.
  if (!is.null(extra_controls)) {
    # If it's a list, add items individually; if single ts, add directly
    if (is.list(extra_controls)) {
      control_list <- c(control_list, extra_controls)
    } else {
      control_list[["ExtraCtrl"]] <- extra_controls
    }
  }
  
  # 3. Intersect and align
  controls_ts <- do.call(ts.intersect, control_list)
  
  # 4. Estimation Loop
  all_results <- list()
  
  for (target_name in var_names) {
    target_series <- lhs_vars_list[[target_name]]
    var_res <- data.frame(Variable = target_name, Horizon = 0:h_max, IRF = NA, SE = NA)
    
    for (h in 0:h_max) {
      # A. LHS at t+h
      Y_future <- stats::lag(target_series, h)
      
      # B. Merge LHS, Shock, and Controls
      # Check for Singularity: If h=0 and target_series is in extra_controls, 
      # this regression is invalid. We assume user handles this logic.
      reg_ts <- ts.intersect(LHS = Y_future, Shock = shock_ts, controls_ts)
      
      # C. Window
      reg_data <- window(reg_ts, start = sample_start, end = sample_end)
      df <- as.data.frame(reg_data)
      
      # D. Seasonality (Monthly Dummies)
      df$Season <- factor(cycle(reg_data))
      
      # E. Run Regression
      fit <- lm(LHS ~ ., data = df)
      
      coeftest_res <- coeftest(fit, vcov = vcovHC(fit, type = "HC1"))
      var_res$IRF[h + 1] <- coef(fit)["Shock"]
      var_res$SE[h + 1]  <- coeftest_res["Shock", "Std. Error"]
    }
    var_res$Lower <- var_res$IRF - 1.96 * var_res$SE
    var_res$Upper <- var_res$IRF + 1.96 * var_res$SE
    all_results[[target_name]] <- var_res
  }
  
  final_df <- do.call(rbind, all_results)
  rownames(final_df) <- NULL
  return(final_df)
}