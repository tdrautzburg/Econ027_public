# Install packages if you haven't already: install.packages(c("lmtest", "sandwich"))
library(lmtest)
library(sandwich)

# Load the dataset
url <- "https://raw.githubusercontent.com/tdrautzburg/Econ027_public/main/Data/merged_olympics_macro.csv"
df <- read.csv(url)

# ==============================================================================
# PART 1: THE TRAP OF SPURIOUS REGRESSIONS (LEVELS)
# ==============================================================================

# 1. FWL Theorem Residuals
# Partial out the Host_nation control from all three main variables
m_fwl_cpi    <- lm(log_CPI ~ Host_nation, data = df, na.action = na.exclude)
m_fwl_medals <- lm(log_Medals ~ Host_nation, data = df, na.action = na.exclude)
m_fwl_sent   <- lm(UMCSENT_Q3 ~ Host_nation, data = df, na.action = na.exclude)

# Extract residuals (na.exclude ensures vectors match the dataframe length)
df$res_cpi    <- resid(m_fwl_cpi)
df$res_medals <- resid(m_fwl_medals)
df$res_sent   <- resid(m_fwl_sent)

# Plotting residuals against time with two y-axes (using base R for simplicity)
par(mar = c(5, 4, 4, 4) + 0.3)
plot(df$Year, df$res_cpi, type = "l", col = "blue", lwd = 2,
     ylab = "Residual log CPI", xlab = "Year", main = "FWL Residuals: CPI & Medals")
par(new = TRUE)
plot(df$Year, df$res_medals, type = "l", col = "red", lwd = 2, axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(df$res_medals, na.rm = TRUE)))
mtext("Residual log Medals", side = 4, line = 3)
legend("topleft", legend = c("Resid CPI", "Resid Medals"), col = c("blue", "red"), lty = 1, lwd = 2)

# Repeat for Sentiment
par(mar = c(5, 4, 4, 4) + 0.3)
plot(df$Year, df$res_sent, type = "l", col = "darkgreen", lwd = 2,
     ylab = "Residual Sentiment", xlab = "Year", main = "FWL Residuals: Sentiment & Medals")
par(new = TRUE)
plot(df$Year, df$res_medals, type = "l", col = "red", lwd = 2, axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(df$res_medals, na.rm = TRUE)))
mtext("Residual log Medals", side = 4, line = 3)
legend("topleft", legend = c("Resid Sentiment", "Resid Medals"), col = c("darkgreen", "red"), lty = 1, lwd = 2)

# 3. Scatter Plots of Levels Residuals
par(mfrow = c(1, 2))
plot(df$res_medals, df$res_cpi, main = "FWL: CPI vs Medals", 
     xlab = "Residual Medals", ylab = "Residual CPI", pch = 19, col = "blue")
plot(df$res_medals, df$res_sent, main = "FWL: Sent vs Medals", 
     xlab = "Residual Medals", ylab = "Residual Sentiment", pch = 19, col = "darkgreen")
par(mfrow = c(1, 1))

# 4. Estimate Full OLS Models with Newey-West HAC Standard Errors
m_cpi_full  <- lm(log_CPI ~ log_Medals + Host_nation, data = df, na.action = na.exclude)
m_sent_full <- lm(UMCSENT_Q3 ~ log_Medals + Host_nation, data = df, na.action = na.exclude)

cat("\n--- HAC Estimates for Levels CPI Model ---\n")
print(coeftest(m_cpi_full, vcov = NeweyWest(m_cpi_full)))

cat("\n--- HAC Estimates for Levels Sentiment Model ---\n")
print(coeftest(m_sent_full, vcov = NeweyWest(m_sent_full)))

# 6. Plot the full model residuals over time
df$res_full_cpi  <- resid(m_cpi_full)
df$res_full_sent <- resid(m_sent_full)

par(mfrow = c(1, 2))
plot(df$Year, df$res_full_cpi, type = "l", main = "CPI Model Residuals", xlab = "Year", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)

plot(df$Year, df$res_full_sent, type = "l", main = "Sentiment Model Residuals", xlab = "Year", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)
par(mfrow = c(1, 1))

# ==============================================================================
# PART 2: STATIONARY TRANSFORM (DIFFERENCES)
# ==============================================================================

# 7. Differenced FWL Residuals & Scatter Plots
m_diff_cpi    <- lm(diff_log_CPI ~ Host_nation, data = df, na.action = na.exclude)
m_diff_medals <- lm(diff_log_Medals ~ Host_nation, data = df, na.action = na.exclude)
m_diff_sent   <- lm(diff_sentiment ~ Host_nation, data = df, na.action = na.exclude)

df$res_diff_cpi    <- resid(m_diff_cpi)
df$res_diff_medals <- resid(m_diff_medals)
df$res_diff_sent   <- resid(m_diff_sent)

par(mfrow = c(1, 2))
plot(df$res_diff_medals, df$res_diff_cpi, main = "Diff FWL: CPI vs Medals", 
     xlab = "Resid Diff Medals", ylab = "Resid Diff CPI", pch = 19, col = "blue")
plot(df$res_diff_medals, df$res_diff_sent, main = "Diff FWL: Sent vs Medals", 
     xlab = "Resid Diff Medals", ylab = "Resid Diff Sentiment", pch = 19, col = "darkgreen")
par(mfrow = c(1, 1))

# 8. Estimate Updated Models in Differences
m_diff_cpi_full  <- lm(diff_log_CPI ~ diff_log_Medals + Host_nation, data = df, na.action = na.exclude)
m_diff_sent_full <- lm(diff_sentiment ~ diff_log_Medals + Host_nation, data = df, na.action = na.exclude)

cat("\n--- HAC Estimates for Differenced CPI Model ---\n")
print(coeftest(m_diff_cpi_full, vcov = NeweyWest(m_diff_cpi_full)))

cat("\n--- HAC Estimates for Differenced Sentiment Model ---\n")
print(coeftest(m_diff_sent_full, vcov = NeweyWest(m_diff_sent_full)))

# 10. Extract and Plot Differenced Residuals
df$res_diff_full_cpi  <- resid(m_diff_cpi_full)
df$res_diff_full_sent <- resid(m_diff_sent_full)

par(mfrow = c(1, 2))
plot(df$Year, df$res_diff_full_cpi, type = "l", main = "Diff CPI Residuals", xlab = "Year", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)

plot(df$Year, df$res_diff_full_sent, type = "l", main = "Diff Sentiment Residuals", xlab = "Year", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)
par(mfrow = c(1, 1))