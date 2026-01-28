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
i_am("SW_OJ.R")

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

# Outsourced Descriptive Plots
source(here("scripts", "SW_OJ_descriptive_plots.R"))

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

# Define HAC Parametersfor sub-scripts
N_obs <- nobs(orange_DLM)
nw_lags <- ceiling(0.75 * N_obs^(1/3))

# Outsourced Scatter Plots
source(here("scripts", "SW_OJ_scatter_plots.R"))

# Outsourced SE calculations for cumulative effects
source(here("scripts", "SW_OJ_cumulative_SE.R"))

##############################
# 4. Replicate HAC (Manual)
##############################

mf <- model.frame(orange_SR)
Y <- model.response(mf)
X_mat <- model.matrix(orange_SR, data = mf) 
X <- X_mat[, 2]
res <- orange_SR$residuals
N <- length(res)

acf_c <- function(x, j) {
  return(t(x[-c(1:j)]) %*% na.omit(Lag(x, j)) / t(x) %*% x)
}

v <- as.vector((X - mean(X)) * res)
var_beta_hat <- 1/N * (1/(N-2) * sum((X - mean(X))^2 * res^2) ) / (1/N * sum((X - mean(X))^2))^2
m <- floor(0.75 * N^(1/3)) 
f_hat_T <- 1 + 2 * sum((m - 1:(m-1))/m * sapply(1:(m - 1), function(i) acf_c(x = v, j = i))) 
manual_se <- sqrt(var_beta_hat * f_hat_T)

cat("Manual Calculation: ", manual_se, "\n")

################################
# 5. Replicate Table 16.1
################################

model_base <- dynlm(FOJC_pctc ~ L(FDD, 0:18))
model_seas <- dynlm(FOJC_pctc ~ L(FDD, 0:18) + season(FDD))

df_adl_ts <- ts.intersect(
  Y    = FOJC_pctc,
  LagY = stats::lag(FOJC_pctc, -1),
  L0   = FDD,
  L1   = stats::lag(FDD, -1),
  L2   = stats::lag(FDD, -2),
  L3   = stats::lag(FDD, -3),
  L4   = stats::lag(FDD, -4),
  L5   = stats::lag(FDD, -5),
  L6   = stats::lag(FDD, -6)
)

df_adl <- as.data.frame(df_adl_ts)
names(df_adl) <- c("Y", "LagY", "L0", "L1", "L2", "L3", "L4", "L5", "L6")
model_adl <- lm(Y ~ ., data = df_adl)

N <- nobs(model_base)
nw_lags     <- ceiling(0.75 * N^(1/3))
vcov_base   <- NeweyWest(model_base, lag = nw_lags, adjust = TRUE)
vcov_double <- NeweyWest(model_base, lag = 2*nw_lags, adjust = TRUE)
vcov_seas   <- NeweyWest(model_seas, lag = nw_lags, adjust = TRUE)
vcov_white  <- vcovHC(model_adl, type = "HC1")

# Extract Table Statistics
target_lags <- c(0:6, 12, 18)
col_lag <- c(); col_dyn <- c(); col_mod <- c()
col_cum_dl <- c(); col_cum_db <- c(); col_cum_ss <- c()
col_adl_w <- c(); col_adl_bs <- c(); col_adl_cm <- c()

adl_coefs <- coef(model_adl)
adl_cumul <- calc_adl_irf(adl_coefs)
boot_out <- run_wild_bootstrap(model_adl, B = 1000)

for (h in target_lags) {
  idx_beta <- 2 + h
  idx_cum  <- 2:(2 + h)
  
  val_dyn      <- coef(model_base)[idx_beta]
  se_dyn       <- sqrt(diag(vcov_base)[idx_beta])
  val_cum_base <- sum(coef(model_base)[idx_cum])
  se_cum_base  <- sqrt(sum(vcov_base[idx_cum, idx_cum]))
  se_cum_dbl   <- sqrt(sum(vcov_double[idx_cum, idx_cum]))
  val_cum_seas <- sum(coef(model_seas)[idx_cum])
  se_cum_seas  <- sqrt(sum(vcov_seas[idx_cum, idx_cum]))
  
  if (h <= 6) {
    term_name <- paste0("L", h)
    val_adl_beta <- adl_coefs[term_name]
    se_adl_white <- sqrt(diag(vcov_white)[term_name])
    se_adl_boot  <- boot_out$se_coef[term_name]
    str_adl_beta <- sprintf("%.2f", val_adl_beta)
    str_se_white <- sprintf("(%.2f)", se_adl_white)
    str_se_boot  <- sprintf("(%.2f)", se_adl_boot)
  } else {
    str_adl_beta <- "-"; str_se_white <- ""; str_se_boot  <- ""
  }
  
  val_adl_cum <- adl_cumul[h + 1]
  se_adl_cum  <- boot_out$se_cumul[h + 1]
  
  col_lag <- c(col_lag, as.character(h), "")
  col_dyn    <- c(col_dyn, sprintf("%.2f", val_dyn), sprintf("(%.2f)", se_dyn))
  col_mod    <- c(col_mod, sprintf("%.2f", val_cum_base), sprintf("(%.2f)", se_cum_base))
  col_cum_dl <- c(col_cum_dl, sprintf("%.2f", val_cum_base), sprintf("(%.2f)", se_cum_base))
  col_cum_db <- c(col_cum_db, sprintf("%.2f", val_cum_base), sprintf("(%.2f)", se_cum_dbl))
  col_cum_ss <- c(col_cum_ss, sprintf("%.2f", val_cum_seas), sprintf("(%.2f)", se_cum_seas))
  col_adl_w  <- c(col_adl_w, str_adl_beta, str_se_white)
  col_adl_bs <- c(col_adl_bs, str_adl_beta, str_se_boot)
  col_adl_cm <- c(col_adl_cm, sprintf("%.2f", val_adl_cum), sprintf("(%.2f)", se_adl_cum))
}

final_tab <- data.frame(
  Lag = col_lag, Dyn_DL = col_dyn, ModReg_Cum = col_mod,
  Cum_DL = col_cum_dl, Cum_DoubleNW = col_cum_db, Cum_Seas = col_cum_ss,
  ADL_Coeff_White = col_adl_w, ADL_Coeff_BS = col_adl_bs, ADL_Cumul_BS = col_adl_cm,
  stringsAsFactors = FALSE
)

# Export results to LaTeX
names(final_tab) <- gsub("_", " ", names(final_tab))
xtab <- xtable(final_tab)
print(xtab, include.rownames = FALSE, file = here("../Results", "OJ_dynamics.tex"))

source(here("scripts", "SW_OJ_data_save.R"))