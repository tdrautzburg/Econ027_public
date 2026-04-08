library(fredr)
library(ggplot2)
library(tidyr)
library(dynlm)
library(zoo)
library(dplyr)
library(lubridate)
library(readr)

# 0. API Key (Uncomment and add yours)
fredr_set_key("bc0eee2e52d37c5e0acfe939cecacb94")

# Make sure the Results directory exists
if (!dir.exists("../Results")) {
  dir.create("../Results", recursive = TRUE)
}

# --- 1. Evaluation Parameters ---
# Stepping monthly from Jan 30, 2023 to April 3, 2026
vintages <- seq(as.Date("2023-01-30"), as.Date("2026-04-03"), by = "month")
vintages <- c(vintages, as.Date("2026-04-03"))
variables <- c("UNRATE", "PAYEMS")
h_max <- 6

# List to hold all generated forecasts
all_forecasts <- list()

# --- 2. Pseudo-Real-Time Forecasting Loop ---
# This is "pseudo" real-time, because we are doing it now after the fact.
for (var in variables) {
  for (v_date in vintages) {
    v_date <- as.Date(v_date)
    
    # Download data AS IT EXISTED on the vintage date
    raw_data <- tryCatch({
      fredr(series_id = var, 
            observation_start = as.Date("2000-01-01"),
            realtime_start = v_date, 
            realtime_end = v_date)
    }, error = function(e) NULL)
    
    # Skip if data pull failed or is empty
    if (is.null(raw_data) || nrow(raw_data) == 0) next
    
    # Convert to Time Series (Monthly)
    ts_data <- ts(raw_data$value, start = c(2000, 1), frequency = 12)
    if (var == "PAYEMS") {
      ts_data <- diff(ts_data) # Convert to absolute monthly change
    }
    
    # Address missing data due to the government shutdown (UNRATE Oct 2025)
    time_idx <- time(ts_data)
    if (var == "UNRATE") {
      # 2025 + 9/12 represents October
      target_idx <- which(abs(time_idx - (2025 + 9/12)) < 0.001)
      if (length(target_idx) > 0) ts_data[target_idx] <- 4.45
    }
    
    # Create estimation copy (Drop COVID anomaly 2020-2021)
    ts_estimation <- ts_data
    ts_estimation[time_idx >= 2020.0 & time_idx < 2022.0] <- NA
    
    # If there are too few observations, abandon this iteration of the loop.
    # FIX: Count non-NA values without using na.omit on a ts object
    if (sum(!is.na(ts_estimation)) < 24) next
    
    # --- Fit Models ---
    # dynlm automatically drops rows containing NAs in the merged lags
    model_ar1 <- dynlm(ts_estimation ~ L(ts_estimation, 1))
    model_ar12 <- dynlm(ts_estimation ~ L(ts_estimation, 1:12))
    
    # --- Forecast Init ---
    fc_ar1 <- numeric(h_max)
    fc_ar12 <- numeric(h_max)
    fc_direct_ar1 <- numeric(h_max)
    fc_direct_ar12 <- numeric(h_max)
    
    last_val <- tail(ts_data, 1)
    
    # === PART A: AR(1) Iterated ===
    c_1 <- coef(model_ar1)[1]
    rho_1 <- coef(model_ar1)[2]
    curr <- last_val
    for(i in 1:h_max){
      curr <- c_1 + rho_1 * curr
      fc_ar1[i] <- curr
    }
    
    # === PART B: AR(12) Iterated (Companion Form) ===
    coeffs <- coef(model_ar12)
    F_matrix <- matrix(0, nrow = 12, ncol = 12)
    F_matrix[1, ] <- coeffs[2:13]
    diag(F_matrix[-1, -ncol(F_matrix)]) <- 1
    
    C_vec <- matrix(0, nrow = 12, ncol = 1)
    C_vec[1, 1] <- coeffs[1]
    
    Z_current <- matrix(rev(tail(ts_data, 12)), ncol = 1)
    for(i in 1:h_max){
      Z_next <- C_vec + (F_matrix %*% Z_current)
      fc_ar12[i] <- Z_next[1,1]
      Z_current <- Z_next
    }
    
    # === PART C: Direct Forecasts ===
    for(h in 1:h_max) {
      # Direct AR(1)
      model_h_ar1 <- dynlm(L(ts_estimation, -h) ~ L(ts_estimation, 0))
      coeffs_1 <- coef(model_h_ar1)
      fc_direct_ar1[h] <- coeffs_1[1] + coeffs_1[2] * last_val
      
      # Direct AR(12)
      model_h_ar12 <- dynlm(L(ts_estimation, -h) ~ L(ts_estimation, 0:11))
      coeffs_12 <- coef(model_h_ar12)
      predictors <- rev(tail(ts_data, 12)) 
      fc_direct_ar12[h] <- coeffs_12[1] + sum(coeffs_12[2:13] * predictors)
    }
    
    # --- Store Results for this Vintage ---
    last_ym <- as.yearmon(time_idx[length(time_idx)])
    # target date = forecast date + horizon
    fc_dates <- as.Date(last_ym + (1:h_max)/12)
    
    temp_df <- data.frame(
      Variable = var,
      Vintage_Date = v_date,
      Target_Date = fc_dates,
      Horizon = 1:h_max,
      Iter_AR1 = fc_ar1,
      Iter_AR12 = fc_ar12,
      Direct_AR1 = fc_direct_ar1,
      Direct_AR12 = fc_direct_ar12
    )
    all_forecasts[[length(all_forecasts) + 1]] <- temp_df
  }
}

# --- 3. Consolidate and Fetch "First-Release Truth" ---
results_df <- bind_rows(all_forecasts) %>%
  pivot_longer(cols = c(Iter_AR1, Iter_AR12, Direct_AR1, Direct_AR12), 
               names_to = "Model", values_to = "Forecast")

truth_list <- list()
for (var in variables) {
  
  # Pull the entire revision history. 
  # Start observation date slightly earlier (Oct 2023) to ensure we have a lag 
  # to calculate the difference for January 2024.
  # Pull the entire revision history. 
  t_data_all <- fredr(series_id = var, 
                      observation_start = as.Date("2023-10-01"),
                      realtime_start = as.Date("2000-01-01")) 
  
  # Employment growth is hard: We need the employment level for t and t-1 
  # from the same vintage to compute the vintage employment growth
  if (var == "PAYEMS") {
    # 1. Group by the release date (vintage) and calculate the MoM change 
    #    using the data *exactly as it was reported together on that day*.
    first_releases <- t_data_all %>%
      arrange(realtime_start, date) %>%
      group_by(realtime_start) %>%
      mutate(value = value - lag(value, order_by = date)) %>%
      ungroup() %>%
      drop_na(value) %>% # Drop the first observation of each vintage which has no lag
      # 2. Now find the very first time a growth number was published for each target month
      group_by(date) %>%
      arrange(realtime_start) %>%
      slice(1) %>%
      ungroup()
    
  } else {
    # For UNRATE, just take the first published value for each observation month
    first_releases <- t_data_all %>%
      group_by(date) %>%
      arrange(realtime_start) %>%
      slice(1) %>%
      ungroup()
  }
  
  truth_list[[var]] <- first_releases %>% 
    select(Target_Date = date, Truth = value) %>% 
    mutate(Variable = var)
}
truth_df <- bind_rows(truth_list)

# Join First-Release Truth to Forecasts
final_df <- results_df %>% 
  left_join(truth_df, by = c("Variable", "Target_Date")) %>%
  drop_na(Truth)

# --- 4. Compute RMSE (Since Jan 2024) ---
rmse_df <- final_df %>%
  filter(Target_Date >= as.Date("2024-01-01")) %>% # <-- CHANGE THIS LINE
  group_by(Variable, Model, Horizon) %>%
  summarize(RMSE = sqrt(mean((Forecast - Truth)^2)), .groups = "drop") %>%
  arrange(Variable, Horizon, RMSE)

print("RMSE Since Jan 2024:")
print(rmse_df)

# --- 5. Export Spreadsheets ---
write_csv(final_df, "../Results/RealTime_Forecasts_Evaluation.csv")
write_csv(rmse_df, "../Results/RealTime_RMSE_Jan2026.csv")
print("Spreadsheets saved to ../Results/")

# --- 6. Spaghetti Plots ---
for (v in variables) {
  plot_data <- final_df %>% filter(Variable == v)
  t_plot_data <- truth_df %>% filter(Variable == v, Target_Date >= as.Date("2023-12-01"))
  
  p <- ggplot() +
    geom_line(data = plot_data, 
              aes(x = Target_Date, y = Forecast, group = interaction(Vintage_Date, Model), color = Model), 
              alpha = 0.4, size = 0.6) +
    geom_line(data = t_plot_data, aes(x = Target_Date, y = Truth), color = "black", size = 1.2) +
    facet_wrap(~Model) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = paste(v, "- Real-Time Forecast Spaghetti Plot"),
         subtitle = "Each colored line is a multi-step forecast generated on a different vintage date.",
         x = "Target Date", y = "Value")
  
  print(p)
  ggsave(paste0("../Results/Spaghetti_", v, "_first.png"), plot = p, width = 10, height = 6)
}
print("Plots saved to ../Results/")

# --- 7. RMSE Bar Plots ---
# Create a bar plot for each variable showing RMSE by Horizon and Model
for (v in variables) {
  
  # Filter RMSE data for the current variable
  rmse_plot_data <- rmse_df %>% filter(Variable == v)
  
  # Only plot if there is RMSE data to show (e.g., if target dates actually reached 2026)
  if (nrow(rmse_plot_data) > 0) {
    p_rmse <- ggplot(rmse_plot_data, aes(x = factor(Horizon), y = RMSE, fill = Model)) +
      # geom_col is equivalent to geom_bar(stat = "identity")
      # position_dodge places the 4 estimator bars side-by-side per horizon
      geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
      scale_fill_brewer(palette = "Set1") +
      theme_minimal() +
      theme(legend.position = "bottom",
            plot.title = element_text(face = "bold")) +
      labs(title = paste(v, "- Forecast RMSE by Horizon and Model"),
           subtitle = "Evaluation for target dates since January 2024",
           x = "Forecast Horizon (Months ahead)",
           y = "Root Mean Squared Error (RMSE)",
           fill = "Estimator")
    
    print(p_rmse)
    ggsave(paste0("../Results/RMSE_BarPlot_", v, "_first.png"), plot = p_rmse, width = 8, height = 5)
  } else {
    print(paste("No RMSE data available to plot for", v, "- check target dates."))
  }
}
print("RMSE bar plots saved to ../Results/")

# --- 8. Compute and Plot Mean Bias (Mean Error) ---

# Calculate Mean Bias (Forecast - Truth)
bias_df <- final_df %>%
  filter(Target_Date >= as.Date("2024-01-01")) %>% # <-- CHANGE THIS LINE TOO
  group_by(Variable, Model, Horizon) %>%
  summarize(Mean_Bias = mean(Forecast - Truth), .groups = "drop") %>%
  arrange(Variable, Horizon)

print("Mean Bias Since Jan 2026:")
print(bias_df)

# Export the Bias table
write_csv(bias_df, "../Results/RealTime_Bias_Jan2024.csv")

# Create a bar plot for each variable showing Bias by Horizon and Model
for (v in variables) {
  
  # Filter Bias data for the current variable
  bias_plot_data <- bias_df %>% filter(Variable == v)
  
  if (nrow(bias_plot_data) > 0) {
    p_bias <- ggplot(bias_plot_data, aes(x = factor(Horizon), y = Mean_Bias, fill = Model)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
      # Add a dashed line at 0 to easily see over vs under prediction
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
      scale_fill_brewer(palette = "Set1") +
      theme_minimal() +
      theme(legend.position = "bottom",
            plot.title = element_text(face = "bold")) +
      labs(title = paste(v, "- Forecast Mean Bias by Horizon and Model"),
           subtitle = "Mean Error (Forecast - Truth) since January 2024",
           x = "Forecast Horizon (Months ahead)",
           y = "Mean Bias (Positive = Over-predict, Negative = Under-predict)",
           fill = "Estimator")
    
    print(p_bias)
    ggsave(paste0("../Results/Bias_BarPlot_", v, "_first.png"), plot = p_bias, width = 8, height = 5)
  } else {
    print(paste("No Bias data available to plot for", v, "- check target dates."))
  }
}
print("Bias spreadsheets and bar plots saved to ../Results/")