library(sandwich)
library(lmtest)

#################
# 0
# Housekeeping
#################

library(here)

# Clear the Environment/console/plots
rm(list = ls())
cat("\014")
if(!is.null(dev.list())) dev.off()

# Set project directory
i_am("HAC_visualized.R")

#################
# 1
# Config 
#################

set.seed(027)
runs <- 2000
T <- 250
p <- 0.5

# Calculate fixed lag length
fixed_lag <- ceiling(.75 * T^(1/3))

#################
# 2
# Analysis 
#################

p_vals_ols <- numeric(runs)
p_vals_hac <- numeric(runs)

for(i in 1:runs) {
  # Generate "Sticky" Data
  x <- numeric(T)
  x[1] <- rnorm(1)
  for(t in 2:T) {
    x[t] <- if(runif(1) < p) x[t-1] else rnorm(1)
  }
  
  model <- lm(x ~ 1)
  
  # 1. Standard OLS p-value
  p_vals_ols[i] <- summary(model)$coefficients[1, 4]
  
  # 2. HAC p-value (Force fixed lag to match Stata's 'newey, lag(4)')
  # We set prewhite = FALSE to match Stata's default behavior
  p_vals_hac[i] <- coeftest(model, vcov = NeweyWest(model, 
                                                    lag = fixed_lag, 
                                                    prewhite = FALSE))[1, 4]
}

# Results
cat("Lag length used:", fixed_lag, "\n")
cat("Type I Error Rate (Standard OLS):", mean(p_vals_ols < 0.05), "\n")
cat("Type I Error Rate (HAC):", mean(p_vals_hac < 0.05), "\n")

# Visualization
# Clear any existing devices to prevent "corrupted" locks
while (!is.null(dev.list())) dev.off()
# Display the plot in R first
par(mfrow=c(1,2))
hist(p_vals_ols, main="OLS p-values", xlab="p-value", breaks=20, col="salmon")
hist(p_vals_hac, main="HAC p-values", xlab="p-value", breaks=20, col="lightblue")
# Save the CURRENT plot window to a file
file_name <- here("../Results", paste0("HAC_hist_p", p*100, "_T", T, ".png"))
dev.copy(png, file_name, width = 800, height = 500)
dev.off()

# Display the plot in R first
par(mfrow=c(1,1))
plot(x[1:min(100, T)], 
     type = "s",       # "s" (stairs) to better show the "stickiness"
     col = "darkblue", 
     lwd = 2,
     main = paste("The 'Sticky' Data (p =", p, ")"),
     xlab = "Time (t)", 
     ylab = "Value of x",
     panel.first = grid())

# Add a horizontal line at the true mean (0)
abline(h = 0, col = "red", lty = 2)

# 3. Save the CURRENT plot window to a file
file_name <- here("../Results", paste0("HAC_series_p", p*100, ".png"))

dev.copy(png, file_name, width = 800, height = 500)
dev.off() # Closes the file device
