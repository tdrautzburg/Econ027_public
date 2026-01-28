# --- Function: Wild Bootstrap Engine ---
run_wild_bootstrap <- function(model_obj, B = 1000, seed = 12345) {
  set.seed(seed)
  y_hat <- fitted(model_obj)
  u_hat <- residuals(model_obj)
  n     <- length(u_hat)
  X_mat <- model.matrix(model_obj)
  
  boot_coefs <- matrix(NA, nrow = B, ncol = length(coef(model_obj)))
  boot_irfs  <- matrix(NA, nrow = B, ncol = 19)
  colnames(boot_coefs) <- names(coef(model_obj))
  
  for (b in 1:B) {
    v_t <- sample(c(-1, 1), n, replace = TRUE)
    y_star <- y_hat + (u_hat * v_t)
    coef_star <- solve(crossprod(X_mat), crossprod(X_mat, y_star))
    coef_vec  <- as.vector(coef_star)
    names(coef_vec) <- colnames(X_mat)
    
    boot_coefs[b, ] <- coef_vec
    boot_irfs[b, ]  <- calc_adl_irf(coef_vec)
  }
  
  return(list(se_coef = apply(boot_coefs, 2, sd), 
              se_cumul = apply(boot_irfs, 2, sd)))
}
