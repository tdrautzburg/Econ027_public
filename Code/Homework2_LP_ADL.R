# ---------------------------------------------------------
# 1. Setup and Dependencies
# ---------------------------------------------------------
library(car)      # For deltaMethod
library(sandwich) # For vcovHC
library(lmtest)   # For coeftest

set.seed(27) # Critical for reproducibility

ADL_CORRECT <- FALSE

# ---------------------------------------------------------
# 2. Simulation (Data Generation)
# ---------------------------------------------------------
# True Parameters
rho1 <- 0.9
rho2 <- -0.5
theta <- 1
sigma_e <- 2 

# Simulation
T_total <- 2200
burn_in <- 200
T_final <- T_total - burn_in
T_short <- 2000

y <- numeric(T_total)
x <- rnorm(T_total, mean = 0, sd = 1)
epsilon <- rnorm(T_total, mean = 0, sd = sigma_e)

# Initialize
y[1] <- 0; y[2] <- 0

# Generate AR(2) process
for (t in 3:T_total) {
  y[t] <- rho1 * y[t-1] + rho2 * y[t-2] + theta * x[t] + epsilon[t]
}

# Create Data Frame
df_all <- data.frame(time = 1:T_total, y = y, x = x)
df_clean <- df_all[(burn_in + 1):T_total, ] # Drop burn-in
rownames(df_clean) <- NULL

# Create Lags (needed for LP and ADL)
# We use a helper to shift vectors
lag_var <- function(v, k) c(rep(NA, k), head(v, -k))

df_clean$y_lag1 <- lag_var(df_clean$y, 1)
df_clean$y_lag2 <- lag_var(df_clean$y, 2)
df_clean$y_lag3 <- lag_var(df_clean$y, 3) # Extra lag for LP control
df_clean$x_lag1 <- lag_var(df_clean$x, 1)
df_clean$x_lag2 <- lag_var(df_clean$x, 2)
df_clean$x_lag3 <- lag_var(df_clean$x, 3)

# Short Sample (First obs)
df_short <- df_clean[1:T_short, ]

# ---------------------------------------------------------
# 3. Calculate TRUE IRF Values (Benchmark)
# ---------------------------------------------------------
# Formulas from Part 1 of your homework
true_irf_0 <- theta
true_irf_4 <- (rho1^4 + 3 * rho1^2 * rho2 + rho2^2) * theta

# ---------------------------------------------------------
# 4. Estimate ADL Model (h=0 and h=4)
# ---------------------------------------------------------

if (ADL_CORRECT) {
  adl_model <- lm(y ~ y_lag1 + y_lag2 + x, data = df_short)
  vcov_adl <- vcovHC(adl_model, type = "HC1")
  
  # Define formulas for Delta Method
  form_h0 <- "x"
  form_h4 <- "x * (y_lag1^4 + 3*y_lag1^2*y_lag2 + y_lag2^2)"
  
  # Compute estimates
  adl_est_0 <- deltaMethod(adl_model, form_h0, vcov. = vcov_adl)
  adl_est_4 <- deltaMethod(adl_model, form_h4, vcov. = vcov_adl)
}  else {
  adl_model <- lm(y ~ y_lag1  + x, data = df_short)
  vcov_adl <- vcovHC(adl_model, type = "HC1")
  
  # Define formulas for Delta Method
  form_h0 <- "x"
  form_h4 <- "x * y_lag1^4 "
  
  # Compute estimates
  adl_est_0 <- deltaMethod(adl_model, form_h0, vcov. = vcov_adl)
  adl_est_4 <- deltaMethod(adl_model, form_h4, vcov. = vcov_adl)
}
# ---------------------------------------------------------
# 5. Estimate LP Model (h=0 and h=4)
# ---------------------------------------------------------
# Corrected LP Function (Handles h=0 and NAs robustly)
estimate_lp <- function(data, h) {
  n <- nrow(data)
  # Shift y forward by h (and pad end with NAs)
  y_target <- c(data$y[(h+1):n], rep(NA, h))
  
  # Run regression with controls
  # We use the existing lags in 'data' to ensure alignment
  lp_model <- lm(y_target ~ x + x_lag1 + x_lag2 + x_lag3 + 
                   y_lag1 + y_lag2 + y_lag3, data = data)
  
  vcov_lp <- vcovHC(lp_model, type = "HC1")
  res <- coeftest(lp_model, vcov. = vcov_lp)["x", , drop=FALSE]
  
  # Return formatted dataframe (Est, SE, 2.5%, 97.5%)
  ci <- confint(lp_model, level = 0.95, vcov. = vcov_lp)["x", ] # using 95% for consistency
  # Note: If your prompt asks for 90%, change level to 0.90
  
  return(c(Estimate = res[1,1], SE = res[1,2], 
           Lower = res[1,1] - 1.96*res[1,2], 
           Upper = res[1,1] + 1.96*res[1,2]))
}

lp_res_0 <- estimate_lp(df_short, 0)
lp_res_4 <- estimate_lp(df_short, 4)

# ---------------------------------------------------------
# 6. Final Comparison Output
# ---------------------------------------------------------
cat("\n=======================================================\n")
cat("               IRF COMPARISON RESULTS                  \n")
cat("=======================================================\n")

# --- HORIZON 0 ---
cat("\n--- HORIZON h = 0 ---\n")
cat(sprintf("TRUE Value:     %.4f\n", true_irf_0))
cat(sprintf("ADL Estimate:   %.4f  [CI: %.4f, %.4f] (Width: %.4f)\n", 
            adl_est_0$Estimate, adl_est_0$`2.5 %`, adl_est_0$`97.5 %`, 
            adl_est_0$`97.5 %` - adl_est_0$`2.5 %`))
cat(sprintf("LP Estimate:    %.4f  [CI: %.4f, %.4f] (Width: %.4f)\n", 
            lp_res_0['Estimate'], lp_res_0['Lower'], lp_res_0['Upper'], 
            lp_res_0['Upper'] - lp_res_0['Lower']))

# Check Coverage
adl_cover_0 <- (true_irf_0 >= adl_est_0$`2.5 %` & true_irf_0 <= adl_est_0$`97.5 %`)
lp_cover_0 <- (true_irf_0 >= lp_res_0['Lower'] & true_irf_0 <= lp_res_0['Upper'])
cat(sprintf("True value in ADL CI? %s\n", adl_cover_0))
cat(sprintf("True value in LP CI?  %s\n", lp_cover_0))

# --- HORIZON 4 ---
cat("\n--- HORIZON h = 4 ---\n")
cat(sprintf("TRUE Value:     %.4f\n", true_irf_4))
cat(sprintf("ADL Estimate:   %.4f  [CI: %.4f, %.4f] (Width: %.4f)\n", 
            adl_est_4$Estimate, adl_est_4$`2.5 %`, adl_est_4$`97.5 %`, 
            adl_est_4$`97.5 %` - adl_est_4$`2.5 %`))
cat(sprintf("LP Estimate:    %.4f  [CI: %.4f, %.4f] (Width: %.4f)\n", 
            lp_res_4['Estimate'], lp_res_4['Lower'], lp_res_4['Upper'], 
            lp_res_4['Upper'] - lp_res_4['Lower']))

# Check Coverage
adl_cover_4 <- (true_irf_4 >= adl_est_4$`2.5 %` & true_irf_4 <= adl_est_4$`97.5 %`)
lp_cover_4 <- (true_irf_4 >= lp_res_4['Lower'] & true_irf_4 <= lp_res_4['Upper'])
cat(sprintf("True value in ADL CI? %s\n", adl_cover_4))
cat(sprintf("True value in LP CI?  %s\n", lp_cover_4))
cat("\n=======================================================\n")