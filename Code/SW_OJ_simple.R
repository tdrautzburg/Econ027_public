# Main source: https://www.econometrics-with-r.org/15.1-the-orange-juice-data.html
# TD/Gemini-written extensions and illustrations

#################
# 0. Package Management
#################
if (!require("pacman")) install.packages("pacman")

# Define all required packages
libs <- c("AER", "dynlm", "orcutt", "nlme", "stargazer", "zoo", "quantmod", 
          "xtable", "dplyr", "lubridate", "janitor", "here", "sandwich", 
          "lmtest", "Matrix", "MatrixModels", "binsreg", "car")

# Install only what is missing, using binary for speed and to avoid 'staged install' locks
missing <- libs[!(libs %in% installed.packages()[, "Package"])]
if(length(missing)) install.packages(missing, type = "binary")

# Load all libraries silently
pacman::p_load(char = libs)

#################
# 1. Housekeeping
#################

# Clear environment and console
rm(list = ls())
cat("\014")
if(!is.null(dev.list())) dev.off()

# Set project directory - ensure this file exists in your root
library(here)
i_am("SW_OJ_simple.R")

# Source functions (helper scripts)
source(here("functions", "f_wild_bootstrap.R"))
source(here("functions", "f_IRFs_from_ADL.R"))

################
# 2. Describe Data
################

data("FrozenJuice")

# Compute price index and returns
FOJCPI <- FrozenJuice[, "price"]/FrozenJuice[, "ppi"]
FOJC_pctc <- 100 * diff(log(FOJCPI))
FDD <- FrozenJuice[, "fdd"]

# Descriptive Plots
# Plot 1: Orange juice price index
plot(as.zoo(FOJCPI),
     col = "steelblue", 
     lwd = 3, 
     xlab = "Date",
     ylab = "Price index", 
     main = "Frozen Concentrated Orange Juice",
     cex.main = 2, cex.lab = 1.5, cex.axis = 1.3)

# Plot 2: Combined Metrics (Vertical Stack)
par(mfrow = c(2, 1), mar = c(3, 4, 2, 1), oma = c(1, 0, 0, 0))

plot(as.zoo(FOJC_pctc),
     col = "steelblue", lwd = 1.5,
     xlab = "", ylab = "Price % Change",
     main = "Price Changes vs. Freezing Degree Days",
     cex.main = 1.2, cex.lab = 1.0, cex.axis = 0.9)

plot(as.zoo(FDD),
     col = "steelblue", lwd = 1.5,
     xlab = "Date", ylab = "FDD",
     main = "", cex.lab = 1.0, cex.axis = 0.9)


# Simple Regression
orange_SR <- dynlm(FOJC_pctc ~ FDD)
coeftest(orange_SR, vcov. = NeweyWest)

#################################
# 3. Estimate Dynamic Causal Effects
#################################

# run dynamic regression
orange_DLM <- dynlm(FOJC_pctc ~ FDD + L(FDD, 1:6))
# compute point estimates of cumulative multipliers
cum_mult <- cumsum(orange_DLM$coefficients[-1])
# rename columns
names(cum_mult) <- paste(0:6, "period CDM", sep = "-")
# display
print(cum_mult)

# Cumulative multipliers via modified regression
cum_mult_reg <- dynlm(FOJC_pctc ~ d(FDD) + d(L(FDD,1:5)) + L(FDD,6))
coeftest(cum_mult_reg, vcov. = NeweyWest)
 