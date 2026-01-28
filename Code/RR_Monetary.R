if (!require("pacman")) install.packages("pacman")
pacman::p_load(AER, dynlm, orcutt, nlme, stargazer, zoo, quantmod, 
               sandwich, lmtest, readxl, xtable, dplyr, lubridate, 
               janitor, here, ggplot2, tidyr)

# Clear Environment
rm(list = ls())
cat("\014")
if(!is.null(dev.list())) dev.off()

# Set project directory (Adjust as needed)
# setwd("YOUR_PATH_HERE") 
i_am("RR_Monetary.R") 

# --- Source Functions ---
# Assuming these files contain the functions defined in our previous steps
source(here("functions", "f_general_wild_bootstrap.R"))
source(here("functions", "f_IRFs_from_general_ADL.R"))
source(here("functions", "f_build_lag_data.R"))
source(here("functions", "f_run_lag_augmented_LP.R")) # Contains 'run_lag_augmented_LP'

################################################################
#
# Config
#
################################################################

ZeroImpactRestriction <- FALSE

################################################################
# 1. Load and Prep Data
################################################################

# Read Data
raw_data <- read_excel(here("../", "data/RomerandRomerDataAppendix.xls"), sheet = "DATA BY MONTH")
raw_data$DateObj <- as.yearmon(raw_data$DATE, format = "%b-%y")

# 1. Extract Raw Series
IP_raw  <- 100*raw_data$LNIPNSA    # Log IP
PPI_raw <- 100*raw_data$LNPPINSA   # Log PPI
DFF_raw <- as.numeric(raw_data$DFF) # Federal Funds Rate (Level)
S_raw   <- raw_data$RESID      # Romer Shock

# 2. Convert to Time Series (Start Jan 1966)
start_date <- c(1966, 1)

IP_ts_level  <- ts(IP_raw, start = start_date, frequency = 12)
PPI_ts_level <- ts(PPI_raw, start = start_date, frequency = 12)
FF_ts_level  <- ts(DFF_raw, start = start_date, frequency = 12) # Level of FFR
S_ts         <- ts(S_raw, start = start_date, frequency = 12)

# 3. Create Differences for ADL Inputs
# The ADL function cumulates the input. 
# To get Level Response of IP/PPI -> Input Delta Log (Growth)
# To get Level Response of FFR    -> Input Delta FFR (Change)
dIP_ts  <- diff(IP_ts_level)
dPPI_ts <- diff(PPI_ts_level)
dFFR_ts <- diff(FF_ts_level) # Change in Funds Rate

################################################################
# 2. Loop for Romer & Romer ADL Estimates (Benchmark)
################################################################

# Ensure the restriction toggle is defined (Default to TRUE if missing)
if (!exists("ZeroImpactRestriction")) ZeroImpactRestriction <- TRUE

# Store results for the final combined comparison later
adl_results_all <- data.frame()

# Define the 3 variables to loop over
# Using Clean Names for nicer plots/filenames
vars_config <- list(
  list(name = "Industrial_Production",  data = dIP_ts,  n_s = 36),
  list(name = "Producer_Price_Index",   data = dPPI_ts, n_s = 48), 
  list(name = "Federal_Funds_Rate",     data = dFFR_ts, n_s = 36)
)

cat("Running ADL Benchmark Estimates (ZeroImpactRestriction =", ZeroImpactRestriction, ")...\n")

for (cfg in vars_config) {
  
  # --- 1. Construct Lag Matrix Manually ---
  
  # A. Create Lags of Dependent Variable (Y) - Always 1 to 24
  y_lags <- list(Y = cfg$data) 
  for(i in 1:24) {
    y_lags[[paste0("LagY_", i)]] <- stats::lag(cfg$data, -i)
  }
  
  # B. Create Lags of Shock (S) - Controlled by ZeroImpactRestriction
  s_lags_list <- list() 
  
  if (ZeroImpactRestriction) {
    # Romer & Romer Replication: Start at Lag 1 (No immediate impact)
    current_s_lags <- 1:cfg$n_s
    for(j in current_s_lags) {
      s_lags_list[[paste0("LagS_", j)]] <- stats::lag(S_ts, -j)
    }
  } else {
    # Unrestricted: Start at Lag 0 (Immediate impact allowed)
    current_s_lags <- 0:cfg$n_s
    for(j in current_s_lags) {
      s_lags_list[[paste0("LagS_", j)]] <- stats::lag(S_ts, -j)
    }
  }
  
  # C. Combine into a Time Series Matrix
  ts_matrix <- do.call(ts.intersect, c(y_lags, s_lags_list))
  
  # --- 2. Window the Data ---
  ts_windowed <- window(ts_matrix, start = c(1970, 1), end = c(1996, 12))
  df_adl <- as.data.frame(ts_windowed)
  
  # --- 3. Estimate Model ---
  model_adl <- lm(Y ~ ., data = df_adl)
  
  # --- 4. Calculate IRF (Point Estimate) ---
  coef_main <- coef(model_adl)
  
  # Pass the correct lag vector (current_s_lags) determined above
  irf_point <- calc_general_irf(coef_main, 
                                s_lags = current_s_lags, 
                                horizon = 48)
  
  # --- 5. Bootstrap Standard Errors ---
  # Pass the same lag vector to the bootstrap engine
  boot_results <- run_general_bootstrap(
    model_obj = model_adl, 
    ar_lags   = 1:24,
    s_lags    = current_s_lags, # Dynamic: 0:36 or 1:36 based on config & boolean
    B         = 200, 
    h         = 48
  )
  irf_se <- boot_results$se_cumul
  
  # --- 6. Store Results ---
  temp_df <- data.frame(
    Variable = cfg$name,
    Method   = "ADL (Romer & Romer)",
    Horizon  = 0:48,
    IRF      = irf_point,
    Lower    = irf_point - 1.96 * irf_se,
    Upper    = irf_point + 1.96 * irf_se
  )
  
  adl_results_all <- rbind(adl_results_all, temp_df)
  
  # --- 7. PLOT AND SAVE INDIVIDUALLY ---
  
  p_single <- ggplot(temp_df, aes(x = Horizon, y = IRF)) +
    # Zero line
    geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
    
    # Confidence Intervals (Dashed Lines)
    geom_line(aes(y = Lower), linetype = "dashed", color = "black") +
    geom_line(aes(y = Upper), linetype = "dashed", color = "black") +
    
    # Point Estimate (Solid Line)
    geom_line(color = "black", size = 1.2) +
    
    # Axis Scales (Match Romer Figure 1)
    scale_x_continuous(breaks = seq(0, 48, by = 6), expand = c(0, 0)) +
    
    # Styling
    theme_bw() + 
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_text(color = "black", size = 11),
      axis.title = element_text(size = 12)
    ) +
    
    # Labels
    labs(
      title = gsub("_", " ", cfg$name), 
      y = "Percent",
      x = "Months after Shock"
    )
  
  # Save Plot (PDF)
  if (ZeroImpactRestriction) {
    file_name <- paste0("RR_ADL_Plot_", cfg$name, "_benchmark.png")
  } else {
    file_name <- paste0("RR_ADL_Plot_", cfg$name, "_unrestricted.png")
  }
  ggsave(filename = here("../Papers","GraphsTables", file_name), plot = p_single, width = 8, height = 5)
  
  cat("Saved plot:", file_name, "\n")
}

################################################################
# 3. Run Lag-Augmented Local Projections (Robustness)
################################################################

cat("Running Lag-Augmented Local Projections...\n")

# A. Prepare List of Target Variables (in LEVELS)
# Note: For LP, we simply project the future Level on the shock.
lp_vars_list <- list(
  "Log_IP"  = IP_ts_level,
  "Log_PPI" = PPI_ts_level,
  "FFR"     = FF_ts_level
)

# B. Define Contemporaneous Controls (Optional)
# WARNING: Including Curr_IP when Target is Log_IP causes singularity at h=0.
# The code below runs, but results at h=0 might be NA or 0 depending on lm() handling.
if (ZeroImpactRestriction) {
  current_controls <- list(
    "Curr_IP"  = IP_ts_level, 
    "Curr_PPI" = PPI_ts_level
  )
} else {
  current_controls <- NULL  
}

# C. Run LP
# Note: Ensure the function 'run_lag_augmented_LP' is loaded from your source file
lp_out_df <- f_run_lag_augmented_LP(
  lhs_vars_list  = lp_vars_list,
  shock_ts       = S_ts,
  extra_controls = current_controls, 
  h_max          = 48,
  n_lag_vars     = 24, # RR use 24 lags of outcome
  n_lag_shock    = 36  # RR use 36 lags of shock
)

# D. Format LP Results to match ADL structure
lp_out_df$Method <- "Lag-Augmented LP"
# Select only relevant columns
lp_clean <- lp_out_df %>% select(Variable, Method, Horizon, IRF, Lower, Upper)

################################################################
# 4. Combine and Plot (Corrected)
################################################################

# 1. Combine ADL and LP results
# Ensure we are using the 'adl_results_all' from the loop above
# and 'lp_clean' from the LP section (if you ran that part).
if (!exists("lp_clean")) {
  # If you haven't run the LP part yet, just plot ADL
  final_plot_data <- adl_results_all
} else {
  final_plot_data <- rbind(adl_results_all, lp_clean)
}

# 2. CLEANING: Remove any NA variables (Fixes the 4th panel issue)
final_plot_data <- final_plot_data[!is.na(final_plot_data$Variable), ]

# 3. Rename Variables for nice labels
# This ensures the panels have professional titles
final_plot_data$Variable <- factor(final_plot_data$Variable, 
                                   levels = c("Industrial_Production", "Log_IP", 
                                              "Producer_Price_Index",  "Log_PPI", 
                                              "Federal_Funds_Rate",    "FFR"),
                                   labels = c("Industrial Production", "Industrial Production",
                                              "Producer Price Index",  "Producer Price Index",
                                              "Federal Funds Rate",    "Federal Funds Rate"))
################################################################
# 4. Combine and Plot (Fixed)
################################################################

# 1. Combine ADL and LP results
if (!exists("lp_clean")) {
  final_plot_data <- adl_results_all
} else {
  final_plot_data <- rbind(adl_results_all, lp_clean)
}

# 2. Robust Renaming (Handles both "Log_IP" and "Industrial_Production")
# We use a safe replacement strategy that won't create NAs
library(dplyr)

final_plot_data <- final_plot_data %>%
  mutate(Variable = case_when(
    # Industrial Production aliases
    Variable %in% c("Log_IP", "Industrial_Production", "IP") ~ "Industrial Production",
    
    # PPI aliases
    Variable %in% c("Log_PPI", "Producer_Price_Index", "PPI") ~ "Producer Price Index",
    
    # FFR aliases
    Variable %in% c("FFR", "Federal_Funds_Rate", "DFF") ~ "Federal Funds Rate",
    
    # Fallback (keep original if no match)
    TRUE ~ as.character(Variable)
  ))

# 3. Factor Order (Optional: Forces specific order of panels)
final_plot_data$Variable <- factor(final_plot_data$Variable, 
                                   levels = c("Industrial Production", "Producer Price Index", "Federal Funds Rate"))

# 4. Create the Combined Plot
combined_plot_1 <- ggplot(final_plot_data, aes(x = Horizon, y = IRF, color = Method)) +
  # Zero line
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  
  # Confidence Intervals (Ribbons with low opacity)
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Method), alpha = 0.1, color = NA) +
  
  # IRF Lines
  geom_line(size = 1) +
  
  # Facet by Variable
  facet_wrap(~Variable, scales = "free_y", ncol = 1) +
  
  # Colors: Black for RR, Red for LP
  scale_color_manual(values = c("ADL (Romer & Romer)" = "black", "Lag-Augmented LP" = "red")) +
  scale_fill_manual(values = c("ADL (Romer & Romer)" = "gray20", "Lag-Augmented LP" = "red")) +
  
  # Labels
  labs(title = "Monetary Policy Shock: Replication vs. Robustness",
       subtitle = "Comparison of Romer & Romer (2004) ADL vs. Lag-Augmented Local Projection",
       x = "Months after Shock",
       y = "Response (Level)") +
  
  # Theme Adjustments
  theme_bw() +
  theme(
    legend.position = "bottom",
    
    # Remove the gray box from the strip labels
    strip.background = element_rect(fill = "white", color = NA), 
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    
    # Clean up grid lines
    panel.grid.minor = element_blank()
  )

# 5. Print and Save
print(combined_plot_1)

# Save the plot
if (ZeroImpactRestriction) {
  ggsave(here("../Papers","GraphsTables", "RR_ADL_LP_benchmark.png"), combined_plot_1, width = 10, height = 8)
} else {
    ggsave(here("../Papers","GraphsTables", "RR_ADL_LP_unrestricted.png"), combined_plot_1, width = 10, height = 8)
}