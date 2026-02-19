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
INCLUDE_CUBIC_TREND <- TRUE  
fredr_set_key("bc0eee2e52d37c5e0acfe939cecacb94")

dir.create("../Results", showWarnings = FALSE, recursive = TRUE)
dir.create("../data", showWarnings = FALSE, recursive = TRUE) # Ensure data dir exists
saved_files <- c()

# -----------------------------------------------------------------------------
# 2. Data Download and Transformation
# -----------------------------------------------------------------------------
raw_gdp <- fredr(series_id = "GDPC1")
raw_gov <- fredr(series_id = "GCEC1")

df_gdp <- data.frame(date = raw_gdp$date, gdp = log(raw_gdp$value))
df_gov <- data.frame(date = raw_gov$date, g = log(raw_gov$value))

df_merged <- merge(df_gdp, df_gov, by = "date")
df_merged <- df_merged[order(df_merged$date), ]
df_merged <- na.omit(df_merged)

ts_mat <- as.matrix(df_merged[, c("g", "gdp")])
start_date <- c(as.numeric(format(df_merged$date[1], "%Y")), 
                as.numeric(substr(quarters(df_merged$date[1]), 2, 2)))
data_ts <- ts(ts_mat, start = start_date, frequency = 4)

# -----------------------------------------------------------------------------
# 3. Shock Identification & CSV Save
# -----------------------------------------------------------------------------
sample_suffix <- if (PRE_COVID_ONLY) "pre_covid" else "full_sample"
cutoff_time_shock <- if (PRE_COVID_ONLY) 2020.0 else end(data_ts)[1] + 1
data_shock <- window(data_ts, end = cutoff_time_shock)

trend_str <- if (INCLUDE_CUBIC_TREND) "+ trend(g) + I(trend(g)^2) + I(trend(g)^3)" else ""

shock_fmla <- as.formula(paste("g ~ L(g, 1:8) + L(gdp, 1:8)", trend_str))
shock_model <- dynlm(shock_fmla, data = data_shock)
shock_resids <- residuals(shock_model) 
dg <- sd(shock_resids)

# --- NEW: Save Shock Residuals to CSV ---
shock_export <- data.frame(
  date = as.yearqtr(time(shock_resids)), 
  shock = as.numeric(shock_resids)
)
write.csv(shock_export, "../data/shock_resids.csv", row.names = FALSE)
cat("Shock residuals saved to ../data/shock_resids.csv\n")

# Plot Raw g vs Shock
ts_filename <- paste0("../Results/TimeSeries_Shock_", sample_suffix, ".png")
par(mfrow = c(2, 1), mar = c(4, 5, 3, 2))
plot(data_shock[, "g"], col = "black", lwd = 2, main = "Log Spending (g)", ylab = "Log Level", xlab = "")
grid()
plot(shock_resids, col = "blue", lwd = 1, type = "h", main = "Shock (Residuals)", ylab = "Dev", xlab = "Time")
abline(h = 0)
grid()
dev.copy(png, filename = ts_filename, width = 800, height = 800); dev.off()
saved_files <- c(saved_files, ts_filename)

# -----------------------------------------------------------------------------
# 4. Identification Visualization (FWL 2x2 Panels)
# -----------------------------------------------------------------------------
viz_horizons <- c(0, 2, 4, 20)
for (var_name in c("g", "gdp")) {
  # We increase the resolution (res) to ensure the font looks crisp
  # We also adjust mar (margins) slightly to prevent labels from being cut off
  png(filename = paste0("../Results/Identification_", var_name, "_", sample_suffix, ".png"), 
      width = 1000, height = 1000, res = 120) 
  
  # cex = 1.2 scales all text up by 20% by default
  par(mfrow = c(2, 2), mar = c(5, 5, 4, 2), cex = 1.1)
  
  for (h in viz_horizons) {
    lhs_fmla <- as.formula(paste0("L(", var_name, ", -", h, ") ~ L(g, 1:8) + L(gdp, 1:8)", trend_str))
    d_sub <- if (PRE_COVID_ONLY) window(data_ts, end = 2020.0 - (h/4)) else data_ts
    lhs_model <- dynlm(lhs_fmla, data = d_sub)
    
    aligned <- merge.zoo(residuals(lhs_model), shock_resids, all = FALSE)
    
    # cex.main: Title size
    # cex.lab:  X and Y axis label size
    # cex.axis: Tick mark number size
    plot(as.numeric(aligned[,2]), as.numeric(aligned[,1]), 
         pch = 16, col = rgb(0.2, 0.2, 0.2, 0.4),
         main = paste(toupper(var_name), "| h =", h), 
         xlab = "Shock", 
         ylab = "LHS Resid",
         cex.main = 1.6, 
         cex.lab = 1.4, 
         cex.axis = 1.2)
    
    fit <- lm(aligned[,1] ~ aligned[,2])
    abline(fit, col = "red", lwd = 3) # Thicker line for Beamer
    
    # Legend text boosted with cex = 1.5
    legend("topleft", 
           legend = paste("Slope =", round(coef(fit)[2], 3)), 
           bty = "n", 
           text.col = "red", 
           cex = 1.5, 
           text.font = 2) # Bold text
    grid(lwd = 1.5)
  }
  
  dev.off()
  viz_fn <- paste0("../Results/Identification_", var_name, "_", sample_suffix, ".png")
  saved_files <- c(saved_files, viz_fn)
}
# -----------------------------------------------------------------------------
# 5. Estimation Loop (Unscaled and Scaled IRFs)
# -----------------------------------------------------------------------------
par(mfrow = c(1, 1)) 
h_max <- 20
shock_var <- "g" 
raw_betas_store <- list()

for (var_name in c("g", "gdp")) {
  betas <- numeric(h_max + 1); ses <- numeric(h_max + 1); horizons <- 0:h_max
  
  for (h in horizons) {
    if (h == 0 && var_name == shock_var) { betas[h+1] <- 1; ses[h+1] <- 0; next }
    fmla <- as.formula(paste0("L(", var_name, ", -", h, ") ~ ", shock_var, " + L(g, 1:8) + L(gdp, 1:8)", trend_str))
    d_sub <- if (PRE_COVID_ONLY) window(data_ts, end = 2020.0 - (h/4)) else data_ts
    model <- dynlm(fmla, data = d_sub)
    ctest <- coeftest(model, vcov. = NeweyWest(model, lag = h + 1, prewhite = FALSE, adjust = TRUE))
    betas[h+1] <- ctest[shock_var, "Estimate"]; ses[h+1] <- ctest[shock_var, "Std. Error"]
  }
  raw_betas_store[[var_name]] <- betas
  
  # --- A. UNSCALED (RED) ---
  fn_u <- paste0("../Results/IRF_", var_name, "_unscaled_", sample_suffix, ".png")
  plot(horizons, betas, type = "n", ylim = c(min(betas-1.65*ses), max(betas+1.65*ses)), 
       main = paste("Response of", toupper(var_name), "(Unscaled)"), xlab = "Horizon", ylab = "Raw Coeff")
  grid(); abline(h = 0, lwd = 2)
  polygon(c(horizons, rev(horizons)), c(betas+1.65*ses, rev(betas-1.65*ses)), col = rgb(1,0,0,0.1), border = NA)
  lines(horizons, betas, col = "red", lwd = 3)
  dev.copy(png, filename = fn_u, width = 800, height = 600); dev.off()
  saved_files <- c(saved_files, fn_u)
  
  # --- B. SCALED (BLUE) ---
  fn_s <- paste0("../Results/IRF_", var_name, "_scaled_", sample_suffix, ".png")
  b_s <- betas * dg; s_s <- ses * dg
  plot(horizons, b_s, type = "n", ylim = c(min(b_s-1.65*s_s), max(b_s+1.65*s_s)), 
       main = paste("Response of", toupper(var_name), "(1 SD Shock)"), xlab = "Horizon", ylab = "Log Dev")
  grid(); abline(h = 0, lwd = 2)
  polygon(c(horizons, rev(horizons)), c(b_s+1.65*s_s, rev(b_s-1.65*s_s)), col = rgb(0,0,1,0.1), border = NA)
  lines(horizons, b_s, col = "blue", lwd = 3)
  dev.copy(png, filename = fn_s, width = 800, height = 600); dev.off()
  saved_files <- c(saved_files, fn_s)
}

# -----------------------------------------------------------------------------
# 6. Multiplier & Final List
# -----------------------------------------------------------------------------
multiplier <- (cumsum(raw_betas_store[["gdp"]]) / cumsum(raw_betas_store[["g"]])) * 5
fn_m <- paste0("../Results/Multiplier_", sample_suffix, ".png")
plot(0:h_max, multiplier, type = "o", pch = 19, col = "darkred", lwd = 3, main = "Cumulative Multiplier (x5)")
grid(); abline(h = 0)
dev.copy(png, filename = fn_m, width = 800, height = 600); dev.off()
saved_files <- c(saved_files, fn_m)

cat("\nDone! Files saved:\n"); print(saved_files)