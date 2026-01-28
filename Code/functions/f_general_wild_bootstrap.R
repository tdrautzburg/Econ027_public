# --- Function: Wild Bootstrap Engine (Corrected) ---
run_general_bootstrap <- function(model_obj, 
                                  ar_lags = 1:24, 
                                  s_lags = 1:36, 
                                  B = 500, 
                                  seed = 12345, 
                                  h = 48) {
  
  set.seed(seed)
  y_hat <- fitted(model_obj)
  u_hat <- residuals(model_obj)
  n     <- length(u_hat)
  X_mat <- model.matrix(model_obj) 
  
  # Pre-calculate X'X inverse for speed
  XX_inv <- solve(crossprod(X_mat))
  
  boot_irfs  <- matrix(NA, nrow = B, ncol = h + 1)
  
  for (b in 1:B) {
    # Rademacher distribution
    v_t <- sample(c(-1, 1), n, replace = TRUE)
    y_star <- y_hat + (u_hat * v_t)
    
    # Fast OLS
    coef_star <- XX_inv %*% crossprod(X_mat, y_star)
    coef_vec  <- as.vector(coef_star)
    names(coef_vec) <- colnames(X_mat)
    
    # PASS THE LAGS HERE so the IRF is calculated correctly
    boot_irfs[b, ]  <- calc_general_irf(coef_vec, 
                                        ar_lags = ar_lags, 
                                        s_lags = s_lags, 
                                        horizon = h)
  }
  
  return(list(se_cumul = apply(boot_irfs, 2, sd)))
}