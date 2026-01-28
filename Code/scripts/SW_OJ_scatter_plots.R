
# Plot 3: Scatter Plot (Slide Optimized)
png(here("../Results", "plot_scatter_price_fdd.png"), width = 1000, height = 800, res = 150)
par(mar = c(5, 5, 4, 2))
plot_data_model <- model.frame(orange_SR)
plot(plot_data_model$FDD, plot_data_model$FOJC_pctc,
     pch = 19, col = adjustcolor("steelblue", alpha.f = 0.5),
     xlab = "Freezing Degree Days", ylab = "Price Change (%)",
     main = "Price Sensitivity to Freezing Weather",
     cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.1)
abline(orange_SR, col = "darkred", lwd = 3)
grid()
dev.off()

# --- PLOT 4: Binscatter Plot with Superimposed OLS ---
# 1. Clean data for binning
df_plot <- data.frame(
  y = as.numeric(FOJC_pctc), 
  x = as.numeric(FDD[-1]) 
) %>% na.omit()

# 2. Run binsreg with quantile binning ("qs")
bin_out <- binsreg(y = df_plot$y, x = df_plot$x, nbins = 20, binspos = "qs")
dots_data <- bin_out$data.plot[[1]]$data.dots

# 3. Create Slide Plot
png(here("../Results", "plot_binscatter_price_fdd.png"), width = 1000, height = 800, res = 150)
par(mar = c(5, 5, 4, 2))
plot(dots_data$x, dots_data$y,
     pch = 19, col = "steelblue", cex = 1.8,
     xlab = "Freezing Degree Days (Quantile Bins)", 
     ylab = "Avg Price Change (%)",
     main = "Binscatter: OLS vs. Conditional Means",
     cex.main = 1.3, cex.lab = 1.1, cex.axis = 1.0)

# Re-estimate lm on clean df to ensure perfect overlay
model_fit <- lm(y ~ x, data = df_plot)
abline(model_fit, col = "darkred", lwd = 3)

legend("topleft", legend = c("Binned Means (20)", "OLS Line"),
       col = c("steelblue", "darkred"), pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 3), bty = "n")
grid()
dev.off()

# --- STATA-TYPE BINSCATTER (Manual Quantile Binning) ---

# 1. Align the data into a clean data frame
# We use as.numeric and handle the diff() lag to match model dimensions
df_stata_bin <- data.frame(
  y = as.numeric(FOJC_pctc),
  x = as.numeric(FDD[-1])
) %>% na.omit()

# 2. Re-estimate OLS on the exact same clean data frame to ensure the line matches
model_stata <- lm(y ~ x, data = df_stata_bin)

# 3. Create 20 equal-sized bins based on quantiles (Stata default)
num_bins <- 20
df_binned_stata <- df_stata_bin %>%
  mutate(bin = ntile(x, num_bins)) %>%
  group_by(bin) %>%
  summarise(
    mean_x = mean(x),
    mean_y = mean(y)
  )

# 4. Save as plot_binscatter_price_fdd_stata.png
png(here("../Results", "plot_binscatter_price_fdd_stata.png"), 
    width = 1000, height = 800, res = 150)
par(mar = c(5, 5, 4, 2))

# Plot binned means
plot(df_binned_stata$mean_x, df_binned_stata$mean_y,
     pch = 19, col = "steelblue", cex = 1.8,
     xlab = "Freezing Degree Days (Quantile Bins)", 
     ylab = "Avg Price Change (%)",
     main = "Stata-style Binscatter (Quantile Bins)",
     cex.main = 1.3, cex.lab = 1.1, cex.axis = 1.0)

# Superimpose the OLS line (Linear property is preserved by quantile binning)
abline(model_stata, col = "darkred", lwd = 3)

# Add clear legend
legend("topleft", 
       legend = c("Binned Means (Quantiles)", "OLS Regression Line"),
       col = c("steelblue", "darkred"), 
       pch = c(19, NA), 
       lty = c(NA, 1), 
       lwd = c(NA, 3), 
       bty = "n")

grid()
dev.off()
