# --- Function: Calculate Cumulative IRF for General ADL(p, q) ---
# Supports arbitrary lags for Y (AR) and S (Shock), including lag 0.
calc_general_irf <- function(coef_vec, ar_lags = 1:24, s_lags = 1:36, horizon = 48) {
  
  # 1. Extract AR coefficients (phi_i for Delta Y_t-i)
  phi <- numeric(max(ar_lags))
  for (i in ar_lags) {
    name <- paste0("LagY_", i)
    if (name %in% names(coef_vec)) phi[i] <- coef_vec[name]
  }
  
  # 2. Extract Shock coefficients (beta_j for S_t-j)
  # beta[1] = Lag 0, beta[2] = Lag 1, etc.
  beta <- numeric(max(s_lags) + 1) 
  
  for (j in s_lags) {
    if (j == 0) {
      # Handle Lag 0: Could be named "S", "Shock", or "LagS_0"
      if ("S" %in% names(coef_vec)) {
        beta[j + 1] <- coef_vec["S"]
      } else if ("Shock" %in% names(coef_vec)) {
        beta[j + 1] <- coef_vec["Shock"]
      } else if ("LagS_0" %in% names(coef_vec)) {
        beta[j + 1] <- coef_vec["LagS_0"]
      }
    } else {
      # Handle Lags 1+: Standard format "LagS_j"
      name <- paste0("LagS_", j)
      if (name %in% names(coef_vec)) beta[j + 1] <- coef_vec[name]
    }
  }
  
  # 3. Recursive Filter (The IRF Calculation)
  irf_delta <- numeric(horizon + 1)
  
  for (t in 1:(horizon + 1)) {
    h <- t - 1 # Current horizon step (0, 1, 2...)
    
    # Direct impact of shock lags (Moving Average part)
    # If horizon h matches a shock lag j, add beta_j
    # Note: h corresponds to index h+1 in beta vector
    term_shock <- if ((h + 1) <= length(beta)) beta[h + 1] else 0
    
    # Auto-regressive feedback (AR part)
    term_ar <- 0
    for (k in seq_along(phi)) {
      if (phi[k] != 0 && (t - k) > 0) {
        term_ar <- term_ar + (phi[k] * irf_delta[t - k])
      }
    }
    
    irf_delta[t] <- term_shock + term_ar
  }
  
  # Return Cumulative Sum
  return(cumsum(irf_delta))
}