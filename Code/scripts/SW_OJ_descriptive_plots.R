# Plot 1: Orange juice price index
png(here("../Results", "plot_price_index.png"), width = 1000, height = 700)
plot(as.zoo(FOJCPI),
     col = "steelblue", 
     lwd = 3, 
     xlab = "Date",
     ylab = "Price index", 
     main = "Frozen Concentrated Orange Juice",
     cex.main = 2, cex.lab = 1.5, cex.axis = 1.3)
dev.off()

# Plot 2: Combined Metrics (Vertical Stack)
png(here("../Results", "plot_combined_metrics.png"), width = 1200, height = 900, res = 150)
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
dev.off()
