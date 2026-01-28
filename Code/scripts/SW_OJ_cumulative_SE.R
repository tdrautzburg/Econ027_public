# 1. Use the robust VCOV matrix from your model
V <- NeweyWest(orange_DLM, lag = nw_lags, adjust = TRUE)

# 2. Identify the indices of the coefficients
idx1 <- "FDD"
idx2 <- "L(FDD, 1:6)1"

# 3. Calculate Variance of the sum
var_sum <- V[idx1, idx1] + V[idx2, idx2] + 2 * V[idx1, idx2]

# 4. Get the Standard Error
se_sum <- sqrt(var_sum)
cat("SE for the sum of coefficients:", se_sum, "\n")

pacman::p_load(car)

# 2. Define the hypothesis as a string
hyp_matrix <- c(0, 1, 1, rep(0, length(coef(orange_DLM)) - 3))

# Perform the test using the matrix form
lh_test <- linearHypothesis(orange_DLM, 
                            hypothesis.matrix = hyp_matrix, 
                            vcov. = NeweyWest(orange_DLM, lag = nw_lags, adjust = TRUE))

print(lh_test)

# 1. Get the point estimate (the sum)
# We calculate this manually from the coefficients
sum_coef <- sum(coef(orange_DLM)[c("FDD", "L(FDD, 1:6)1")])

# 2. Extract the F-statistic from the linearHypothesis object
# The F-stat is usually in the 3rd or 4th column of the 2nd row
f_stat <- lh_test$F[2]

# 3. Back out the Standard Error
# Since sqrt(F) = t = estimate / SE
se_backed_out <- abs(sum_coef) / sqrt(f_stat)

cat("Point Estimate (Sum):", sum_coef, "\n")
cat("Extracted SE:", se_backed_out, "\n")

# Use deltaMethod to get the estimate and SE in one data frame
# Note: Using indices (e.g., V1, V2) is safer for complex dynlm names
res_delta <- deltaMethod(orange_DLM, "FDD + `L(FDD, 1:6)1` ", 
                         vcov. = NeweyWest(orange_DLM, lag = nw_lags, adjust = TRUE))

print(res_delta)
# Extract just the SE
se_final <- res_delta$SE
cat("Delta-Method SE:", se_final, "\n")
