spurious_sim <- function(T = 100, Reps = 1000, RW = TRUE) {
  set.seed(27)
  
  # Directory and naming logic
  suffix <- if(RW) "RW" else "IID"
  if(!dir.exists("../Results")) dir.create("../Results")
  
  t_homo <- numeric(Reps)
  t_nw   <- numeric(Reps)
  y <- numeric(T); z <- numeric(T); res <- numeric(T)
  
  for (i in 1:Reps) {
    if (RW == TRUE) {
      y <- cumsum(rnorm(T)); z <- cumsum(rnorm(T))
    } else {
      y <- rnorm(T); z <- rnorm(T)
    }
    
    mod <- lm(y ~ z)
    res <- residuals(mod)
    
    # Statistics
    t_homo[i] <- summary(mod)$coefficients[2, "t value"]
    x_demean <- z - mean(z)
    w <- x_demean * res
    sigma2_nw <- (sum(w^2) + 2 * sum(w[-1] * w[-length(w)])) / sum(x_demean^2)^2
    t_nw[i] <- coef(mod)[2] / sqrt(sigma2_nw)
  }
  
  # --- 1. Plot Histogram (Display & Save) ---
  # Open a new window for the histogram
  dev.new() 
  hist(abs(t_homo), breaks = 30, col = "skyblue", border = "white",
       main = paste("Distribution of |t-stat| (T =", T, ", Type =", suffix, ")"),
       xlab = "|t-statistic|")
  abline(v = 1.645, col = "red", lwd = 2, lty = 2)
  
  # Save the current window to PNG
  dev.copy(png, filename = paste0("../Results/SpuriousHist_", suffix, ".png"), 
           width = 800, height = 600)
  dev.off() # Closes the file device but keeps the window open (usually)
  
  # --- 2. Plot Time Series (Display & Save) ---
  dev.new()
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 5))
  
  # Top Panel: Dual Y-Axis
  plot(y, type = "l", col = "blue", lwd = 2, ylab = "y (Blue)", xlab = "Time", 
       main = paste("Last Simulation Path (", suffix, ")"))
  par(new = TRUE)
  plot(z, type = "l", col = "darkgreen", lwd = 2, axes = FALSE, xlab = "", ylab = "")
  axis(side = 4, col = "darkgreen", col.axis = "darkgreen")
  mtext("z (Green)", side = 4, line = 3, col = "darkgreen")
  
  # Bottom Panel: Residuals
  plot(res, type = "l", col = "red", lwd = 2, ylab = "Residual", xlab = "Time",
       main = "Regression Residuals (e_t)")
  abline(h = 0, lty = 3)
  
  # Save the current window to PNG
  dev.copy(png, filename = paste0("../Results/SpuriousReg_", suffix, ".png"), 
           width = 800, height = 800)
  dev.off()
  
  # Console Output
  cat("\n--- Simulation Summary (", suffix, ") ---\n")
  cat("Fraction |t| > 1.645 (Homoskedastic):", mean(abs(t_homo) > 1.645), "\n")
  cat("Fraction |t| > 1.645 (Newey-West):   ", mean(abs(t_nw) > 1.645), "\n")
}

# Run the simulations
spurious_sim(T = 100, Reps = 1000, RW = TRUE)
spurious_sim(T = 100, Reps = 1000, RW = FALSE)