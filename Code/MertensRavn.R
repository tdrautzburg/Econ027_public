# -----------------------------------------------------------------------------
# 1. Reset and Setup
# -----------------------------------------------------------------------------
rm(list = ls())     # Clear workspace
graphics.off()      # Close all plot windows
cat("\014")         # Clear console

# Load necessary libraries
library(readxl)
library(zoo)
library(dynlm)
library(sandwich)
library(lmtest)

# Define paths
dir.create("../Results", showWarnings = FALSE, recursive = TRUE)
dir.create("../data", showWarnings = FALSE, recursive = TRUE)

MR_SHOCKS_ONLY <- TRUE 

# -----------------------------------------------------------------------------
# 2. Data Download and Transformation
# -----------------------------------------------------------------------------
# Read the Excel file
# Using .xlsx based on standard format; ensure file is in ../data/
file_path <- "../data/jme2014_data.xls" 
df <- read_excel(file_path, sheet = "Sheet1", skip = 1)

# Rename columns to make formulas cleaner
# Expected columns: Date, Tax Revenues, Govt Spending, Output, Tax Narrative
colnames(df)[1] <- "Date"
colnames(df)[2] <- "tax"
colnames(df)[3] <- "gov"
colnames(df)[4] <- "output"
# The shock column has spaces in the header
colnames(df)[5] <- "shock_col" 

# Ensure numeric and transform
# Log levels * 100 implies coefficients are percentage changes
df$tax <- as.numeric(df$tax) * 100
df$gov <- as.numeric(df$gov) * 100
df$output <- as.numeric(df$output) * 100
# "All Romer Shocks" column
if (MR_SHOCKS_ONLY) {
  df$shock <- -as.numeric(df$shock_col)
  suffix <- "_MR"
} else {
  df$shock <- -as.numeric(df$`All Romer Shocks`)
  suffix <- "_RR"
}

# Create Time Series Object for dynlm
# Assuming the first date is 1950.00 (1950 Q1)
start_year <- floor(df$Date[1])
start_qtr <- (df$Date[1] %% 1) * 4 + 1
data_ts <- ts(df[, c("tax", "gov", "output", "shock")], 
              start = c(start_year, start_qtr), 
              frequency = 4)

# -----------------------------------------------------------------------------
# 3. Estimation Loop (Local Projections)
# -----------------------------------------------------------------------------
outcomes <- c("tax", "gov", "output") 
h_max <- 20
saved_files <- c()

# Store betas for multiplier calculation later
store_betas <- list()

par(mfrow = c(1, 1))

for (var_name in outcomes) {
  
  betas <- numeric(h_max + 1)
  ses   <- numeric(h_max + 1)
  horizons <- 0:h_max
  
  cat(paste0("Estimating IRF for target: ", var_name, "...\n"))
  
  for (h in horizons) {
    
    # Construct Formula
    # Includes cubic trend of the shock/time
    fmla_str <- paste0("L(", var_name, ", -", h, ") ~ ",
                       "shock + trend(shock) + I(trend(shock)^2) + I(trend(shock)^3) +",
                       "L(tax, 1:4) + ",
                       "L(gov, 1:4) + ",
                       "L(output, 1:4) + ",
                       "L(shock, 1:4)")
    
    model <- dynlm(as.formula(fmla_str), data = data_ts)
    
    # Inference using Newey-West SE
    nw_vcov <- NeweyWest(model, lag = h + 1, prewhite = FALSE, adjust = TRUE)
    ctest <- coeftest(model, vcov. = nw_vcov)
    
    betas[h+1] <- ctest["shock", "Estimate"]
    ses[h+1]   <- ctest["shock", "Std. Error"]
  }
  
  # Store coefficients
  store_betas[[var_name]] <- betas
  
  # --- Plotting Individual IRFs ---
  upper <- betas + 1.645 * ses
  lower <- betas - 1.645 * ses
  plot_filename <- paste0("../Results/LP_IRF_", var_name, suffix, ".png")
  
  y_min <- min(lower)
  y_max <- max(upper)
  
  plot(horizons, betas, type = "n", ylim = c(y_min, y_max),
       main = paste0("Response of ", toupper(var_name), " to Tax Narrative Shock"),
       xlab = "Horizon (Quarters)", ylab = "Response [%]")
  grid()
  abline(h = 0, lwd = 1)
  polygon(c(horizons, rev(horizons)), c(upper, rev(lower)), col = rgb(0, 0, 1, 0.1), border = NA)
  lines(horizons, betas, col = "blue", lwd = 2)
  
  dev.copy(png, filename = plot_filename, width = 800, height = 600); dev.off()
  saved_files <- c(saved_files, plot_filename)
}

# -----------------------------------------------------------------------------
# 4. PDV Multiplier Calculation
# -----------------------------------------------------------------------------
cat("\nCalculating PDV Multiplier...\n")

# Extract stored IRFs
betas_tax <- store_betas[["tax"]]
betas_output <- store_betas[["output"]]

# Discount factor: 1 / (1 + 0.01)^h
discount_factors <- 1 / ((1.01) ^ (0:h_max))

# Compute PDVs (Present Discounted Values)
# cumsum gives the running sum of discounted effects up to horizon h
pdv_tax_irf <- cumsum(betas_tax * discount_factors)
pdv_output_irf <- cumsum(betas_output * discount_factors)

# Compute Multiplier
# Ratio of PDVs scaled by 1 / (Tax/Output Ratio) to convert elasticities to dollar units.
# Tax/Output Ratio = 0.175
ratio_tax_to_output <- 0.175
multiplier <- -(pdv_output_irf / pdv_tax_irf) / ratio_tax_to_output

# -----------------------------------------------------------------------------
# 5. Plot Multiplier
# -----------------------------------------------------------------------------
mult_filename <- paste0("../Results/LP_Multiplier_", var_name, suffix, ".png")

# Determine Y-axis limits (Starts at 0, extends to Max value + 10% buffer)
y_max_limit <- max(multiplier, na.rm = TRUE) * 1.1
# Ensure the top of the graph is at least slightly above 0 if the multiplier is very small
if (y_max_limit <= 0) y_max_limit <- 0.1 

plot(0:h_max, multiplier, type = "l", col = "darkred", lwd = 3,
     ylim = c(0, y_max_limit), # Vertical axis starts at 0
     main = "PDV Tax Multiplier",
     xlab = "Horizon (Quarters)",
     ylab = "Multiplier ($)")
grid()
abline(h = 0, lwd = 1)

# Save
dev.copy(png, filename = mult_filename, width = 800, height = 600); dev.off()
saved_files <- c(saved_files, mult_filename)

cat("\nAnalysis Complete. Saved plots:\n")
print(saved_files)