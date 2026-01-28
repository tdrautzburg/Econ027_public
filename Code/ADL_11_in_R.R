# --- 1. Generate Data for ADL(1,1) ---
set.seed(027)
n <- 100
phi <- 0.7      # Persistence in y
beta0 <- 1.5    # Contemporaneous effect of x
beta1 <- -0.8   # Lagged effect of x

x <- rnorm(n)
u <- rnorm(n)
y <- numeric(n)

# Initial value for y
y[1] <- (beta0 * x[1] + u[1]) / (1 - phi) 

# Simulation loop
for(t in 2:n) {
  y[t] <- phi * y[t-1] + beta0 * x[t] + beta1 * x[t-1] + u[t]
}

# --- 2. Manual OLS using lm() ---
# We must align the lags manually. 
# Note: we lose the first observation due to the lag.
y_t    <- y[2:n]
y_lag1 <- y[1:(n-1)]
x_t    <- x[2:n]
x_lag1 <- x[1:(n-1)]

fit_lm <- lm(y_t ~ y_lag1 + x_t + x_lag1)
summary(fit_lm)

# --- 3. Using dynlm() ---
# Install if needed: install.packages("dynlm")
library(dynlm)

# dynlm requires a time-series (ts) object
data_ts <- ts(data.frame(y, x))

# L() is the lag operator; it handles alignment automatically
fit_dyn <- dynlm(y ~ L(y) + x + L(x), data = data_ts)
summary(fit_dyn)

# Compare coefficients (they will be identical)
coef(fit_lm)
coef(fit_dyn)