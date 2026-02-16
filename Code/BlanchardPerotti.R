# -----------------------------------------------------------------------------
# 1. Reset and Setup
# -----------------------------------------------------------------------------
rm(list = ls())     # Clear workspace
graphics.off()      # Close all plot windows
cat("\014")         # Clear console

library(fredr)
library(sandwich)
library(lmtest)
library(dynlm) 
library(zoo) 

# --- CONFIGURATION SWITCHES ---
PRE_COVID_ONLY <- TRUE 
INCLUDE_CUBIC_TREND <- FALSE  # Set to TRUE to include t, t^2, t^3

# Placeholder for API Key
# Get your FRED API key here
# https://fred.stlouisfed.org/docs/api/api_key.html
fredr_set_key("your API key here")

# Ensure Results directory exists
dir.create("../Results", showWarnings = FALSE, recursive = TRUE)

# Track saved files for final output
saved_files <- c()

# -----------------------------------------------------------------------------
# 2. Data Download and Transformation
# -----------------------------------------------------------------------------

# (1) Read time series
raw_gdp <- fredr(series_id = "GDPC1")
raw_gov <- fredr(series_id = "GCEC1")

# (2) Log transform
df_gdp <- data.frame(date = raw_gdp$date, gdp = log(raw_gdp$value))
df_gov <- data.frame(date = raw_gov$date, g = log(raw_gov$value))

# (3) Merge and Clean
df_merged <- merge(df_gdp, df_gov, by = "date")
df_merged <- df_merged[order(df_merged$date), ]
df_merged <- na.omit(df_merged)

# (4) Create Time Series Object
ts_mat <- as.matrix(df_merged[, c("g", "gdp")])

start_year <- as.numeric(format(df_merged$date[1], "%Y"))
start_qtr  <- as.numeric(substr(quarters(df_merged$date[1]), 2, 2))

data_ts <- ts(ts_mat, start = c(start_year, start_qtr), frequency = 4)

# -----------------------------------------------------------------------------
# 3. Shock Identification (Model & Time Series Plot)
# -----------------------------------------------------------------------------

sample_suffix <- if (PRE_COVID_ONLY) "pre_covid" else "full_sample"

# Define cutoff for shock estimation
cutoff_time_shock <- if (PRE_COVID_ONLY) 2020.0 else end(data_ts)[1] + 1
data_shock <- window(data_ts, end = cutoff_time_shock)

# --- Define Trend String ---
# We use trend(g) which dynlm interprets as the time index
trend_str <- ""
if (INCLUDE_CUBIC_TREND) {
  trend_str <- "+ trend(g) + I(trend(g)^2) + I(trend(g)^3)"
}

# Estimate innovation model
# Note: We must include the trend here to ensure residuals are orthogonal to it
shock_fmla_str <- paste("g ~ L(g, 1:8) + L(gdp, 1:8)", trend_str)
shock_model <- dynlm(as.formula(shock_fmla_str), data = data_shock)

dg <- sd(residuals(shock_model))
shock_resids <- residuals(shock_model)

cat(paste("Standard Deviation of Shock (dg):", round(dg, 5), "\n"))
if (INCLUDE_CUBIC_TREND) cat("Note: Cubic deterministic trend included in estimation.\n")

# --- PLOT 1: Time Series (Raw g vs Residuals) ---
ts_filename <- paste0("../Results/TimeSeries_Shock_Residuals_", sample_suffix, ".png")

par(mfrow = c(2, 1), mar = c(4, 5, 3, 2) + 0.1)

# Top Panel
plot(data_shock[, "g"], col = "black", lwd = 2,
     main = "Log Government Spending (g)", 
     ylab = "Log Level", xlab = "")
grid()

# Bottom Panel
plot(shock_resids, col = "blue", lwd = 1, type = "h",
     main = "Structural Shock (Residuals)", 
     ylab = "Deviation", xlab = "Time")
abline(h = 0, col = "black")
grid()

dev.copy(png, filename = ts_filename, width = 800, height = 800)
dev.off()
saved_files <- c(saved_files, ts_filename)

par(mfrow = c(1, 1)) # Reset layout

# -----------------------------------------------------------------------------
# 4. Estimation Loop (IRFs)
# -----------------------------------------------------------------------------

variables_to_plot <- c("g", "gdp") 
h_max <- 20
shock_var <- "g" 
raw_betas_store <- list()

for (var_name in variables_to_plot) {
  
  cat(paste("Estimating IRF for LHS variable:", var_name, "\n"))
  
  betas <- numeric(h_max + 1)
  ses   <- numeric(h_max + 1)
  horizons <- 0:h_max
  
  for (h in horizons) {
    
    if (h == 0 && var_name == shock_var) {
      betas[h + 1] <- 1
      ses[h + 1]   <- 0
      next
    }
    
    # Formula: Uses negative lags to "time-travel" -- pushes the LHS into the future
    lhs_str <- paste0("L(", var_name, ", -", h, ")")
    
    # Append trend string to RHS
    rhs_str <- paste(shock_var, "+ L(g, 1:8) + L(gdp, 1:8)", trend_str)
    
    fmla_str <- paste(lhs_str, "~", rhs_str)
    
    # Windowed Data
    if (PRE_COVID_ONLY) {
      cutoff_time <- 2020.0 - (h / 4)
      data_subset <- window(data_ts, end = cutoff_time)
      model <- dynlm(as.formula(fmla_str), data = data_subset)
    } else {
      model <- dynlm(as.formula(fmla_str), data = data_ts)
    }
    
    # Standard Errors
    nw_vcov <- NeweyWest(model, lag = h + 1, prewhite = FALSE, adjust = TRUE)
    ctest   <- coeftest(model, vcov. = nw_vcov)
    
    betas[h + 1] <- ctest[shock_var, "Estimate"]
    ses[h + 1]   <- ctest[shock_var, "Std. Error"]
  }
  
  raw_betas_store[[var_name]] <- betas
  
  # ===========================================================================
  # A. UNSCALED IRF (Raw Coefficients)
  # ===========================================================================
  
  lower_raw <- betas - 1.65 * ses
  upper_raw <- betas + 1.65 * ses
  
  plot_title_raw <- paste0("Response of ", toupper(var_name), " (Unscaled)")
  filename_raw <- paste0("../Results/IRF_", var_name, "_unscaled_", sample_suffix, ".png")
  
  # Render
  par(mar = c(5, 5, 4, 2) + 0.1)
  plot(horizons, betas, type = "n",
       ylim = c(min(lower_raw), max(upper_raw)),
       main = plot_title_raw,
       xlab = "Horizon (h)",
       ylab = "Log Deviation (Raw)",
       cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.2)
  
  grid(lwd = 1.5)
  abline(h = 0, col = "black", lwd = 2)
  polygon(c(horizons, rev(horizons)), c(upper_raw, rev(lower_raw)),
          col = rgb(1, 0, 0, 0.2), border = NA) # Red shading for unscaled
  lines(horizons, betas, col = "red", lwd = 3)
  
  # Save
  dev.copy(png, filename = filename_raw, width = 800, height = 600)
  dev.off()
  saved_files <- c(saved_files, filename_raw)
  
  # ===========================================================================
  # B. SCALED IRF (1 SD Shock)
  # ===========================================================================
  
  betas_scaled <- betas * dg
  ses_scaled   <- ses * dg
  
  lower <- betas_scaled - 1.65 * ses_scaled
  upper <- betas_scaled + 1.65 * ses_scaled
  
  plot_title_scaled <- paste0("Response of ", toupper(var_name), " (1 SD Shock)")
  filename_scaled <- paste0("../Results/IRF_", var_name, "_scaled_", sample_suffix, ".png")
  
  # Render
  plot(horizons, betas_scaled, type = "n",
       ylim = c(min(lower), max(upper)),
       main = plot_title_scaled,
       xlab = "Horizon (h)",
       ylab = "Log Deviation (Scaled)",
       cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.2)
  
  grid(lwd = 1.5)
  abline(h = 0, col = "black", lwd = 2)
  polygon(c(horizons, rev(horizons)), c(upper, rev(lower)),
          col = rgb(0, 0, 1, 0.2), border = NA) # Blue shading for scaled
  lines(horizons, betas_scaled, col = "blue", lwd = 3)
  
  # Save
  dev.copy(png, filename = filename_scaled, width = 800, height = 600)
  dev.off()
  saved_files <- c(saved_files, filename_scaled)
}

# -----------------------------------------------------------------------------
# 5. Cumulative Multiplier Calculation and Plot
# -----------------------------------------------------------------------------

cat("Computing Cumulative Multipliers...\n")

b_gdp <- raw_betas_store[["gdp"]]
b_g   <- raw_betas_store[["g"]]

cum_gdp <- cumsum(b_gdp)
cum_g   <- cumsum(b_g)

multiplier <- (cum_gdp / cum_g) * 5

# --- PLOT 4: Multiplier ---
mult_filename <- paste0("../Results/Multiplier_", sample_suffix, ".png")

par(mar = c(5, 5, 4, 2) + 0.1)
plot(horizons, multiplier, type = "o", pch = 19, col = "darkred", lwd = 3,
     main = "Cumulative Spending Multiplier (x5)",
     xlab = "Horizon (h)",
     ylab = "Cumulative Multiplier",
     cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.2)

grid(lwd = 1.5)
abline(h = 0, col = "black", lwd = 1)
abline(h = 0, col = "gray", lty = 2) 

dev.copy(png, filename = mult_filename, width = 800, height = 600)
dev.off()
saved_files <- c(saved_files, mult_filename)

# -----------------------------------------------------------------------------
# 6. Final Output
# -----------------------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Processing Complete. The following plots have been saved:\n")
cat("------------------------------------------------------------\n")
for (f in saved_files) {
  cat(paste("  -", f, "\n"))
}
cat("------------------------------------------------------------\n")