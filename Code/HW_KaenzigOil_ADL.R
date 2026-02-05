library(fredr)
require("purrr")
require("dplyr")
require("lubridate")
require("janitor")
require("ggplot2") # Added explicit ggplot2 loading

#################
# 0
# Housekeeping
#################

library(here)

# Clear the Environment/console/plots
rm(list = ls())
cat("\014")
if(!is.null(dev.list())) dev.off()

# Set project directory
i_am("HW_KaenzigOil_ADL.R")

# --- Source Functions ---
source(here("functions", "f_general_wild_bootstrap.R"))
source(here("functions", "f_IRFs_from_general_ADL.R"))
source(here("functions", "f_build_lag_data.R"))
source(here("functions", "f_run_lag_augmented_LP.R"))

################################################################
#
# Config
#
################################################################

#################
# 1
# Load data
#################

final_data <- read.csv(here("../data", "data_oil_Macro.csv"))
# IP transform
final_data$log_ip <- 100*log(final_data$ip)
# WTI transform
final_data$log_p <- 100*log(final_data$oilprice/final_data$coreCPI)

S_raw <- final_data$oil_supply_news_shock

# 2. Convert to Time Series (Start Jan 1966)
start_date <- c(1975, 1)

IP_ts_level  <- ts(final_data$log_ip, start = start_date, frequency = 12)
P_ts_level   <- ts(final_data$log_p, start = start_date, frequency = 12)
S_ts         <- ts(S_raw, start = start_date, frequency = 12)

dIP_ts  <- diff(IP_ts_level)
dP_ts   <- diff(P_ts_level)

################################################################
#
# PRE-CALCULATION: DETERMINE SCALING FACTOR
#
################################################################
cat("Calculating Scaling Factor based on Real Oil Price Response...\n")

# We run a reference ADL for Real Oil Price (dP_ts) to find the impact of S_t
# Using 12 lags as the standard specification
ref_n_y <- 12
ref_n_s <- 12

# A. Lags of Y (Oil Price Diff)
ref_y_lags <- list(Y = dP_ts) 
for(i in 1:ref_n_y) {
  ref_y_lags[[paste0("LagY_", i)]] <- stats::lag(dP_ts, -i)
}

# B. Lags of S
ref_s_lags_list <- list() 
for(j in 0:ref_n_s) {
  ref_s_lags_list[[paste0("LagS_", j)]] <- stats::lag(S_ts, -j)
}

# C. Matrix & Window
ref_matrix <- do.call(ts.intersect, c(ref_y_lags, ref_s_lags_list))
ref_windowed <- window(ref_matrix, start = c(1970, 1), end = c(1996, 12))
ref_df <- as.data.frame(ref_windowed)

# D. Estimate
ref_model <- lm(Y ~ ., data = ref_df)

# E. Get Impact Coefficient (LagS_0)
impact_coef <- coef(ref_model)["LagS_0"]

# F. Calculate Scaling Factor
# We want the impact to be 5.3%. Currently it is 'impact_coef'.
# Scaled * impact_coef = 5.3  =>  Scaled = 5.3 / impact_coef
scale_factor <- 5.3 / impact_coef

cat(paste0("Impact Coefficient: ", round(impact_coef, 4), "\n"))
cat(paste0("Scaling Factor calculated: ", round(scale_factor, 4), "\n"))

################################################################
#
# MAIN LOOP
#
################################################################

adl_results_all <- data.frame()

vars_config <- list(
  list(name = "Industrial_Production_lvl",   data = IP_ts_level, n_s = 12,  n_y = 12),  
  list(name = "Industrial_Production_lvl",   data = IP_ts_level, n_s = 12,  n_y = 24),  
  list(name = "Industrial_Production_diff",   data = dIP_ts, n_s = 12,  n_y = 24),
  list(name = "Industrial_Production_diff",  data = dIP_ts,  n_s = 2,  n_y = 2),
  list(name = "Industrial_Production_diff",  data = dIP_ts,  n_s = 2,  n_y = 1),
  list(name = "Industrial_Production_diff",   data = dIP_ts, n_s = 1,  n_y = 1),
  list(name = "Real_Oil_Price_diff", data=dP_ts, n_s = 12,  n_y = 24),
  list(name = "Real_Oil_Price_diff", data=dP_ts, n_s = 1,  n_y = 1)
)

cat("\nRunning ADL Benchmark Estimates Loop...\n")

for (cfg in vars_config) {
  
  # --- 1. Construct Lag Matrix Manually ---
  y_lags <- list(Y = cfg$data) 
  for(i in 1:cfg$n_y) {
    y_lags[[paste0("LagY_", i)]] <- stats::lag(cfg$data, -i)
  }
  
  s_lags_list <- list() 
  current_s_lags <- 0:cfg$n_s
  for(j in current_s_lags) {
    s_lags_list[[paste0("LagS_", j)]] <- stats::lag(S_ts, -j)
  }
  
  ts_matrix <- do.call(ts.intersect, c(y_lags, s_lags_list))
  ts_windowed <- window(ts_matrix, start = c(1970, 1), end = c(1996, 12))
  df_adl <- as.data.frame(ts_windowed)
  
  # --- 3. Estimate Model ---
  model_adl <- lm(Y ~ ., data = df_adl)
  
  # --- 4. Calculate IRF (Point Estimate) ---
  coef_main <- coef(model_adl)
  irf_point <- calc_general_irf(coef_main, 
                                s_lags = current_s_lags, 
                                horizon = 48)
  
  # --- 5. Bootstrap Standard Errors ---
  boot_results <- run_general_bootstrap(
    model_obj = model_adl, 
    ar_lags    = 1:24,
    s_lags     = current_s_lags, 
    B          = 200, 
    h          = 48
  )
  irf_se <- boot_results$se_cumul
  
  # --- 6. Store RAW Results ---
  temp_df <- data.frame(
    Variable = cfg$name,
    Method   = "ADL (Romer & Romer)",
    Horizon  = 0:48,
    IRF      = irf_point,
    Lower    = irf_point - 1.96 * irf_se,
    Upper    = irf_point + 1.96 * irf_se
  )
  
  adl_results_all <- rbind(adl_results_all, temp_df)
  
  # --- 7. Plotting Function ---
  # Defined as a helper closure to avoid repeating code for Scaled vs Raw
  create_plot <- function(data_df, title_suffix) {
    ggplot(data_df, aes(x = Horizon, y = IRF)) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
      geom_line(aes(y = Lower), linetype = "dashed", color = "black") +
      geom_line(aes(y = Upper), linetype = "dashed", color = "black") +
      geom_line(color = "black", size = 1.2) +
      scale_x_continuous(breaks = seq(0, 48, by = 6), expand = c(0, 0)) +
      theme_bw() + 
      theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(color = "black", size = 11),
        axis.title = element_text(size = 12)
      ) +
      labs(
        title = paste0(gsub("_", " ", cfg$name), title_suffix), 
        y = "Percent",
        x = "Months after Shock"
      )
  }
  
  # --- 7A. Save RAW Plot ---
  p_single <- create_plot(temp_df, " (Raw)")
  file_name <- paste0("Kaenzig_ADL_Plot_", cfg$name, "_xlag", cfg$n_s, "_ylag", cfg$n_y, ".png")
  ggsave(filename = here("../Results", "Homework", file_name), plot = p_single, width = 8, height = 5)
  
  # --- 7B. Save SCALED Plot ---
  # Apply Scaling Factor to Estimates and SEs
  temp_df_scaled <- temp_df
  temp_df_scaled$IRF   <- temp_df$IRF * scale_factor
  temp_df_scaled$Lower <- temp_df$Lower * scale_factor
  temp_df_scaled$Upper <- temp_df$Upper * scale_factor
  
  p_scaled <- create_plot(temp_df_scaled, " (Scaled to 5.3% Oil Price)")
  file_name_scaled <- paste0("Kaenzig_ADL_Plot_", cfg$name, "_xlag", cfg$n_s, "_ylag", cfg$n_y, "_scaled.png")
  ggsave(filename = here("../Results", "Homework", file_name_scaled), plot = p_scaled, width = 8, height = 5)
  
  cat("Saved plots (Raw & Scaled) for:", cfg$name, "(", cfg$n_s, "lags)\n")
}