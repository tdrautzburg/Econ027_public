# --- Construct the Data Frame Robustly ---
# We use ts.intersect to handle the shifting and alignment of 60 variables
build_lag_data <- function(Y_ts, S_ts, n_ar=24, n_s=36) {
  
  lag_list <- list(Y = Y_ts)
  
  # Add AR Lags (Delta y_{t-1} to Delta y_{t-24})
  for(i in 1:n_ar) {
    lag_list[[paste0("LagY_", i)]] <- stats::lag(Y_ts, -i)
  }
  
  # Add Shock Lags (S_{t-1} to S_{t-36})
  # Note: Image says sum j=1 to 36. 
  for(j in 1:n_s) {
    lag_list[[paste0("LagS_", j)]] <- stats::lag(S_ts, -j)
  }
  
  # Intersect to align dates
  combined_ts <- do.call(ts.intersect, lag_list)
  df <- as.data.frame(combined_ts)
  
  # Add Monthly Dummies (Seasonality)
  # We extract month from the row names or index
  # Note: We drop one dummy (January) automatically via lm's intercept, 
  # or we can create them explicitly. Let's use factor(cycle).
  # But we need the cycle corresponding to the *aligned* data.
  # Since 'combined_ts' is a ts object, cycle() works perfectly.
  df$Season <- factor(cycle(combined_ts))
  
  return(df)
}
