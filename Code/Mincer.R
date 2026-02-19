# 1. Setup Libraries
if (!require("wooldridge")) install.packages("wooldridge")
if (!require("sandwich")) install.packages("sandwich")
if (!require("lmtest")) install.packages("lmtest")
if (!require("car")) install.packages("car") # Needed for linearHypothesis

library(wooldridge)
library(sandwich)
library(lmtest)
library(car)


# Clear the Environment/console/plots
rm(list = ls())
cat("\014")
if(!is.null(dev.list())) dev.off()

data("wage1")

# ---------------------------------------------------------
# PART 1: The Models (Simple vs. Full)
# ---------------------------------------------------------

# Model A: Simple Regression (No Controls)
# This suffers from Omitted Variable Bias (OVB)
simple_model <- lm(lwage ~ educ, data = wage1)

# Model B: Full Mincer Regression
# Controlling for Experience (Quad) and Tenure
full_model <- lm(lwage ~ educ + exper +  tenure + female, data = wage1)

# ---------------------------------------------------------
# PART 2: Robust Standard Errors & Comparisons
# ---------------------------------------------------------

# Calculate Robust SEs for the Full Model
robust_se <- vcovHC(full_model, type = "HC1")
simple_se <- summary(simple_model)$coefficients["educ", "Std. Error"]
full_se_default <- summary(full_model)$coefficients["educ", "Std. Error"]
full_se_robust  <- sqrt(diag(robust_se))["educ"]

cat("\n--- Joint Significance Test (Wald Test) ---\n")
# Null Hypothesis: exper = 0 AND female  = 0 AND tenure = 0
# We use the robust covariance matrix (vcov = robust_vcov)
joint_test <- linearHypothesis(full_model, 
                               c("exper = 0", "female = 0", "tenure = 0"), 
                               vcov = robust_se)
print(joint_test)


# ---------------------------------------------------------
# PART 3: Frisch-Waugh-Lovell (FWL) Theorem
# ---------------------------------------------------------

# Step A: Partial out controls from Y
y_resid <- residuals(lm(lwage ~ exper + tenure + female, data = wage1))

# Step B: Partial out controls from X (Education)
x_resid <- residuals(lm(educ ~ exper  + tenure + female, data = wage1))

# Step C: Regress Y-residuals on X-residuals
fwl_model <- lm(y_resid ~ x_resid)
fwl_se <- summary(fwl_model)$coefficients["x_resid", "Std. Error"]

y_raw <- wage1$lwage #Otherwise, would have to move x_resid to the wage1 data frame
# AlternativeC: Regress raw Y on X-residuals
fwl_model_alt <- lm(y_raw ~ x_resid)
fwl_se_alt <- summary(fwl_model_alt)$coefficients["x_resid", "Std. Error"]

# ---------------------------------------------------------
# PART 4: The Comparison Output
# ---------------------------------------------------------
cat("\n===============================================\n")
cat("      COMPARISON OF ESTIMATES & ERRORS\n")
cat("===============================================\n")

cat("\n1. COEFFICIENTS (Omitted Variable Bias Check)\n")
cat(sprintf("  Simple Model (educ only):  %.5f\n", coef(simple_model)["educ"]))
cat(sprintf("  Full Model (w/ controls):  %.5f\n", coef(full_model)["educ"]))
cat("  -> Note: The simple model underestimates the return to education.\n")
cat(sprintf("  Residual Model:  %.5f\n", coef(fwl_model)["x_resid"]))
cat(sprintf("  Residual Model (raw LHS):  %.5f\n", coef(fwl_model_alt)["x_resid"]))
cat("  -> Note: Same results as full model.\n")

cat("\n2. STANDARD ERRORS (educ coefficient)\n")
cat(sprintf("  Full Model (Default SE):   %.5f\n", full_se_default))
cat(sprintf("  Full Model (Robust SE):    %.5f\n", full_se_robust))
cat(sprintf("  FWL Manual Step (SE):      %.5f\n", fwl_se))
cat("  -> Note: Manual FWL SE is slightly wrong (too small) because it\n")
cat("     uses incorrect Degrees of Freedom (n-2 instead of n-k-1).\n")

# ---------------------------------------------------------
# PART 5: Visualizing the Residualized Data
# ---------------------------------------------------------

plot(x_resid, y_resid,
     main = "FWL: Return to Education (Partialled Out)",
     xlab = "Education (Residuals after controls)",
     ylab = "Log Wage (Residuals after controls)",
     col = rgb(0, 0, 1, 0.5), pch = 16, 
     las = 1)

abline(fwl_model, col = "red", lwd = 2)
grid()
legend("topleft", legend = "Partial Regression Line",
       col = "red", lwd = 2, bty = "n")
