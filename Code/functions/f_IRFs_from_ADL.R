# --- Function: Calculate Cumulative IRF ---
calc_adl_irf <- function(coef_vec, n_lags = 18) {
  # Flexible grep to find Phi (lagged Y) and Betas (lagged X)
  phi   <- coef_vec[grep("LagY", names(coef_vec))]
  betas <- coef_vec[grep("L[0-9]+", names(coef_vec))] # Matches L0, L1...L6
  
  irf <- numeric(n_lags + 1)
  for (i in 1:(n_lags + 1)) {
    lag_h <- i - 1
    # We map lag_h to column names "L0", "L1", etc.
    beta_name <- paste0("L", lag_h)
    b_val <- if (beta_name %in% names(betas)) betas[beta_name] else 0
    prev_effect <- if (i == 1) 0 else irf[i-1]
    irf[i] <- b_val + (phi * prev_effect)
  }
  return(cumsum(irf))
}
