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
INCLUDE_CUBIC_TREND <- FALSE  
USE_UNEMP_INDICATOR <- TRUE   # Set to FALSE to use NBER USRECQ instead
UNEMP_THRESHOLD <- 6.5        # Quarterly average threshold for high unemployment
NLAGS <- 4                    # Switch for number of RHS lags

fredr_set_key("bc0eee2e52d37c5e0acfe939cecacb94")

dir.create("../Results", showWarnings = FALSE, recursive = TRUE)
dir.create("../data", showWarnings = FALSE, recursive = TRUE) 
saved_files <- c()

# -----------------------------------------------------------------------------
# 2. Data Download and Transformation
# -----------------------------------------------------------------------------
raw_gdp <- fredr(series_id = "GDPC1")
raw_gov <- fredr(series_id = "GCEC1")

# Fetch state indicator based on switch
if (USE_UNEMP_INDICATOR) {
  # Download UNRATE, automatically convert to quarterly average via FRED API
  raw_state <- fredr(series_id = "UNRATE", frequency = "q", aggregation_method = "avg")
  df_state <- data.frame(date = raw_state$date, R = ifelse(raw_state$value > UNEMP_THRESHOLD, 1, 0))
  state_suffix <- "HighUnemp"
  state_1_label <- paste0("High Unemp (> ", UNEMP_THRESHOLD, "%)")
  state_0_label <- paste0("Low Unemp (<= ", UNEMP_THRESHOLD, "%)")
} else {
  raw_state <- fredr(series_id = "USRECQ")
  df_state <- data.frame(date = raw_state$date, R = raw_state$value)
  state_suffix <- "Recession"
  state_1_label <- "Recession"
  state_0_label <- "Expansion"
}

df_gdp <- data.frame(date = raw_gdp$date, gdp = log(raw_gdp$value))
df_gov <- data.frame(date = raw_gov$date, g = log(raw_gov$value))

# Merge all three series
df_merged <- merge(df_gdp, df_gov, by = "date")
df_merged <- merge(df_merged, df_state, by = "date")
df_merged <- df_merged[order(df_merged$date), ]
df_merged <- na.omit(df_merged)

ts_mat <- as.matrix(df_merged[, c("g", "gdp", "R")])
start_date <- c(as.numeric(format(df_merged$date[1], "%Y")), 
                as.numeric(substr(quarters(df_merged$date[1]), 2, 2)))
data_ts <- ts(ts_mat, start = start_date, frequency = 4)

# -----------------------------------------------------------------------------
# 3. Shock Identification & CSV Save (State-Dependent)
# -----------------------------------------------------------------------------
sample_suffix <- if (PRE_COVID_ONLY) "pre_covid" else "full_sample"
cutoff_time_shock <- if (PRE_COVID_ONLY) 2020.0 else end(data_ts)[1] + 1
data_shock <- window(data_ts, end = cutoff_time_shock)

trend_str <- if (INCLUDE_CUBIC_TREND) "+ trend(g) + I(trend(g)^2) + I(trend(g)^3)" else ""

# State-dependent shock regression: interacting R with dynamic lags
shock_fmla <- as.formula(paste0("g ~ R * (L(g, 1:", NLAGS, ") + L(gdp, 1:", NLAGS, ")) ", trend_str))
shock_model <- dynlm(shock_fmla, data = data_shock)
shock_resids <- residuals(shock_model) 
dg <- sd(shock_resids)

csv_name <- paste0("../data/shock_resids_state_dep_", state_suffix, "_", NLAGS, "lags.csv")
shock_export <- data.frame(
  date = as.yearqtr(time(shock_resids)), 
  shock = as.numeric(shock_resids)
)
write.csv(shock_export, csv_name, row.names = FALSE)
cat(paste("Shock residuals saved to", csv_name, "\n"))

# Plot Raw g vs Shock
ts_filename <- paste0("../Results/TimeSeries_Shock_StateDep_", state_suffix, "_", sample_suffix, "_", NLAGS, "lags.png")
par(mfrow = c(2, 1), mar = c(4, 5, 3, 2))
plot(data_shock[, "g"], col = "black", lwd = 2, main = "Log Spending (g)", ylab = "Log Level", xlab = "")
grid()
plot(shock_resids, col = "darkmagenta", lwd = 1, type = "h", 
     main = paste("State-Dependent Shock (", state_suffix, ", ", NLAGS, " Lags)", sep=""), ylab = "Dev", xlab = "Time")
abline(h = 0)
grid()
dev.copy(png, filename = ts_filename, width = 800, height = 800); dev.off()
saved_files <- c(saved_files, ts_filename)

# -----------------------------------------------------------------------------
# 4. Estimation Loop (State-Dependent IRFs)
# -----------------------------------------------------------------------------
par(mfrow = c(1, 1)) 
h_max <- 20
shock_var <- "g" 

raw_betas_0_store <- list() # Store baseline (R=0)
raw_betas_1_store <- list() # Store high state (R=1)

for (var_name in c("g", "gdp")) {
  b_0 <- numeric(h_max + 1); se_0 <- numeric(h_max + 1)
  b_1 <- numeric(h_max + 1); se_1 <- numeric(h_max + 1)
  horizons <- 0:h_max
  
  for (h in horizons) {
    if (h == 0 && var_name == shock_var) { 
      b_0[h+1] <- 1; se_0[h+1] <- 0
      b_1[h+1] <- 1; se_1[h+1] <- 0
      next 
    }
    
    # R interacts with the shock AND dynamic number of lags
    fmla <- as.formula(paste0("L(", var_name, ", -", h, ") ~ R * (", shock_var, " + L(g, 1:", NLAGS, ") + L(gdp, 1:", NLAGS, ")) ", trend_str))
    d_sub <- if (PRE_COVID_ONLY) window(data_ts, end = 2020.0 - (h/4)) else data_ts
    model <- dynlm(fmla, data = d_sub)
    
    vmat <- NeweyWest(model, lag = h + 1, prewhite = FALSE, adjust = TRUE)
    cfs <- coef(model)
    
    name_g <- shock_var
    name_Rg <- paste0("R:", shock_var)
    if (!name_Rg %in% names(cfs)) name_Rg <- paste0(shock_var, ":R")
    
    # --- State 0 (Low Unemp / Expansion) ---
    b_0[h+1] <- cfs[name_g]
    se_0[h+1] <- sqrt(vmat[name_g, name_g])
    
    # --- State 1 (High Unemp / Recession) ---
    b_1[h+1] <- cfs[name_g] + cfs[name_Rg]
    se_1[h+1] <- sqrt(vmat[name_g, name_g] + vmat[name_Rg, name_Rg] + 2 * vmat[name_g, name_Rg])
  }
  
  raw_betas_0_store[[var_name]] <- b_0
  raw_betas_1_store[[var_name]] <- b_1
  
  # --- PLOT SCALED IRFs (Both States) ---
  fn_s <- paste0("../Results/IRF_", var_name, "_StateDep_", state_suffix, "_", sample_suffix, "_", NLAGS, "lags.png")
  
  bs_0 <- b_0 * dg; ss_0 <- se_0 * dg
  bs_1 <- b_1 * dg; ss_1 <- se_1 * dg
  
  y_min <- min(c(bs_0 - 1.65*ss_0, bs_1 - 1.65*ss_1), na.rm = TRUE)
  y_max <- max(c(bs_0 + 1.65*ss_0, bs_1 + 1.65*ss_1), na.rm = TRUE)
  
  plot(horizons, bs_0, type = "n", ylim = c(y_min, y_max), 
       main = paste("Response of", toupper(var_name), "to 1 SD Shock (", NLAGS, " Lags)"), 
       xlab = "Horizon", ylab = "Log Dev")
  grid(); abline(h = 0, lwd = 2)
  
  # State 1 Plot (Red)
  polygon(c(horizons, rev(horizons)), c(bs_1+1.65*ss_1, rev(bs_1-1.65*ss_1)), col = rgb(1,0,0,0.15), border = NA)
  lines(horizons, bs_1, col = "red", lwd = 3, lty = 2)
  
  # State 0 Plot (Blue)
  polygon(c(horizons, rev(horizons)), c(bs_0+1.65*ss_0, rev(bs_0-1.65*ss_0)), col = rgb(0,0,1,0.15), border = NA)
  lines(horizons, bs_0, col = "blue", lwd = 3, lty = 1)
  
  legend("topright", legend = c(state_0_label, state_1_label), 
         col = c("blue", "red"), lwd = 3, lty = c(1, 2), bty = "n")
  
  dev.copy(png, filename = fn_s, width = 800, height = 600); dev.off()
  saved_files <- c(saved_files, fn_s)
}

# -----------------------------------------------------------------------------
# 5. Multipliers & Final List
# -----------------------------------------------------------------------------
multiplier_0 <- (cumsum(raw_betas_0_store[["gdp"]]) / cumsum(raw_betas_0_store[["g"]])) * 5
multiplier_1 <- (cumsum(raw_betas_1_store[["gdp"]]) / cumsum(raw_betas_1_store[["g"]])) * 5

fn_m <- paste0("../Results/Multiplier_StateDep_", state_suffix, "_", sample_suffix, "_", NLAGS, "lags.png")
y_lim_m <- range(c(multiplier_0, multiplier_1), na.rm = TRUE)

plot(0:h_max, multiplier_0, type = "o", pch = 19, col = "blue", lwd = 3, 
     ylim = y_lim_m, main = "Cumulative Multiplier (x5) by State", xlab="Horizon", ylab="Multiplier")
lines(0:h_max, multiplier_1, type = "o", pch = 17, col = "red", lwd = 3, lty = 2)
grid(); abline(h = 0, lwd=2)
legend("topleft", legend = c(state_0_label, state_1_label), 
       col = c("blue", "red"), lwd = 3, lty = c(1, 2), pch=c(19, 17), bty = "n")

dev.copy(png, filename = fn_m, width = 800, height = 600); dev.off()
saved_files <- c(saved_files, fn_m)

cat("\nDone! Files saved:\n"); print(saved_files)