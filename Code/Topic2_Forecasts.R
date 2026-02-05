library(fredr)
library(ggplot2)
library(tidyr) # For pivoting data
library(dynlm)
library(zoo) # dynlm relies on zoo/ts infrastructure

# 0. API Key (Uncomment and add yours)
# https://fred.stlouisfed.org/docs/api/api_key.html
fredr_set_key("YOUR API KEY HERE")


DO_U <- FALSE

# --- 1. Download & Prep Data ---
if (DO_U) {
  raw_data <- fredr(series_id = "UNRATE", observation_start = as.Date("2000-01-01"))
} else {
  raw_data <- fredr(series_id = "PAYEMS", observation_start = as.Date("2000-01-01"))
}

# Convert to Time Series (Monthly)
ts_data <- ts(raw_data$value, start = c(2000, 1), frequency = 12)
if (DO_U==FALSE) {
  ts_data <- diff(ts_data)
}

# Address missing data due to the government shutdown
if (DO_U) {
  window(ts_data, start = c(2025, 10), end = c(2025, 10)) <- 4.45
}

# Create a copy for estimation that treats 2020-2021 as "Missing" (NA)
# dynlm/lm will automatically drop these NAs when fitting the model
ts_estimation <- ts_data
  window(ts_estimation, start = c(2020, 1), end = c(2021, 12)) <- NA

# --- 2. Fit Models (Simple OLS with lags) ---

# AR(1): y_t = c + rho*y_{t-1}
# L(ts_estimation, 1) is the 1st lag
model_ar1 <- dynlm(ts_estimation ~ L(ts_estimation, 1))

# AR(12): y_t = c + phi_1*y_{t-1} + ... + phi_12*y_{t-12}
model_ar12 <- dynlm(ts_estimation ~ L(ts_estimation, 1:12))

# --- 3. Forecast ---
# We need the *latest* available data from the real series (not the one with NAs)
# to start our forecast.

h <- 6 # Forecast horizon

# === PART A: AR(1) Simple Iteration ===
# Get coefficients
c_1   <- coef(model_ar1)[1] # Intercept
rho_1 <- coef(model_ar1)[2] # Slope

# Start with the very last observation
last_val <- tail(ts_data, 1)
fc_ar1   <- numeric(h)

curr <- last_val
for(i in 1:h){
  # y_{t+1} = c + rho * y_t
  next_val <- c_1 + rho_1 * curr
  fc_ar1[i] <- next_val
  curr <- next_val # Update for next step
}

# === PART B: AR(12) as VAR(1) in Companion Form ===
# Equation: Z_t = C + F * Z_{t-1}
# Where Z_t is the vector of the last 12 values stacked.

# 1. Construct the Companion Matrix (F)
coeffs <- coef(model_ar12)
intercept <- coeffs[1]
phis      <- coeffs[2:13] # The 12 AR coefficients

# Create 12x12 matrix of zeros
F_matrix <- matrix(0, nrow = 12, ncol = 12)

# Top row is the AR coefficients
F_matrix[1, ] <- phis
# Sub-diagonal is 1s (to shift lags down: y_{t-1} becomes y_{t-2})
if(nrow(F_matrix) > 1) {
  diag(F_matrix[-1, -ncol(F_matrix)]) <- 1
}

# 2. Construct the Constant Vector (C_vec)
# Intercept goes in top row, 0 everywhere else
C_vec <- matrix(0, nrow = 12, ncol = 1)
C_vec[1, 1] <- intercept

# 3. Construct Initial State Vector (Z_current)
# We need the last 12 observations (reversed order: y_T at top, y_{T-11} at bottom)
Z_current <- matrix(rev(tail(ts_data, 12)), ncol = 1)

# 4. Iterate Matrix Multiplication
fc_ar12 <- numeric(h)

for(i in 1:h){
  # VAR(1) step: Z_{t+1} = C + F * Z_t
  Z_next <- C_vec + (F_matrix %*% Z_current)
  
  # The top element is our forecast for y_{t+1}
  fc_ar12[i] <- Z_next[1,1]
  
  # Update state
  Z_current <- Z_next
}


# Part C. Direct Forecasting Loop ---
# We estimate a unique model for each horizon h = 1 to 6.
# For direct forecasting, the RHS is always current information (Lags 0 to p-1).
# The LHS is the future value (Lead h).

horizons <- 6
fc_direct_ar1  <- numeric(horizons)
fc_direct_ar12 <- numeric(horizons)

# We use the FULL current dataset to generate the forecasts
current_data <- ts_data 

for(h in 1:horizons) {
  
  # === Model A: Direct AR(1) ===
  # Equation: y_{t+h} = alpha + beta * y_t
  # In dynlm, L(x, -h) is the LEAD (future value) of x.
  # L(x, 0) is the current value.
  model_h_ar1 <- dynlm(L(ts_estimation, -h) ~ L(ts_estimation, 0))
  
  # Predict: Intercept + Slope * (Latest Actual Value)
  coeffs_1 <- coef(model_h_ar1)
  fc_direct_ar1[h] <- coeffs_1[1] + coeffs_1[2] * tail(current_data, 1)
  
  
  # === Model B: Direct AR(12) ===
  # Equation: y_{t+h} = alpha + beta_1*y_t + ... + beta_12*y_{t-11}
  # We regress the future value (LHS) on the current 12 lags (RHS).
  model_h_ar12 <- dynlm(L(ts_estimation, -h) ~ L(ts_estimation, 0:11))
  
  # Predict:
  coeffs_12 <- coef(model_h_ar12)
  intercept <- coeffs_12[1]
  betas     <- coeffs_12[2:13] # The 12 coefficients
  
  # Get the 12 most recent actual observations
  # Note: dynlm orders lags 0, 1, 2... so we need data ordered Newest -> Oldest
  predictors <- rev(tail(current_data, 12)) 
  
  # Calculate: Intercept + Sum(Coeffs * Predictors)
  fc_direct_ar12[h] <- intercept + sum(betas * predictors)
}

# --- Output Results ---
print("Direct Forecasts (AR1, h=1..6):")
print(fc_direct_ar1)

print("Direct Forecasts (AR12, h=1..6):")
print(fc_direct_ar12)

# --- Output Results ---
print("AR(1) Forecast:")
print(fc_ar1)

print("AR(12) Forecast (Companion Form):")
print(fc_ar12)

# --- 5. Plotting & Saving ---

# A. Prepare the Data Objects
# We slice the last 12 months of ACTUAL history for context
# (start time is calculated dynamically based on the end of the data)
end_time <- time(ts_data)[length(ts_data)]
start_time_plot <- end_time - (11/12) # Go back 11 months (total 12 obs)
recent_history <- window(ts_data, start = start_time_plot)

# Create TS objects for the forecasts so they plot on the correct future dates
# We determine the start date for the forecasts (next month)
last_date <- as.numeric(end(recent_history)) # e.g. c(2025, 10)
fc_start_year <- if(last_date[2] == 12) last_date[1] + 1 else last_date[1]
fc_start_month <- if(last_date[2] == 12) 1 else last_date[2] + 1

# Convert the 4 forecast vectors into time series objects
ts_fc_iter_ar1    <- ts(fc_ar1, start = c(fc_start_year, fc_start_month), frequency = 12)
ts_fc_iter_ar12   <- ts(fc_ar12, start = c(fc_start_year, fc_start_month), frequency = 12)
ts_fc_direct_ar1  <- ts(fc_direct_ar1, start = c(fc_start_year, fc_start_month), frequency = 12)
ts_fc_direct_ar12 <- ts(fc_direct_ar12, start = c(fc_start_year, fc_start_month), frequency = 12)

# 1. Combine all series into one object
# ts.union automatically handles the date alignment and fills gaps with NA
combined_ts <- ts.union(
  History = recent_history, 
  Iter_AR1 = ts_fc_iter_ar1, 
  Iter_AR12 = ts_fc_iter_ar12,
  Direct_AR1 = ts_fc_direct_ar1, 
  Direct_AR12 = ts_fc_direct_ar12
)

# 2. Convert to a Data Frame for ggplot
df <- data.frame(combined_ts)

# 3. Add a proper Date column
# We use 'as.yearmon' from zoo to convert the decimal time (2025.833) to a Date
df$Date <- as.Date(as.yearmon(time(combined_ts)))

# 4. Pivot to "Long" format (ggplot loves long data)
df_long <- pivot_longer(df, cols = -Date, names_to = "Model", values_to = "Rate")

# Remove rows with NA (e.g., where History doesn't overlap with Forecasts)
df_long <- na.omit(df_long) 

# 1. Create the Base Plot (Common to both)
p <- ggplot(df_long, aes(x = Date, y = Rate, color = Model, linetype = Model)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c(
    "History" = "black", 
    "Iter_AR1" = "blue", "Direct_AR1" = "blue",
    "Iter_AR12" = "red", "Direct_AR12" = "red"
  )) +
  scale_linetype_manual(values = c(
    "History" = "solid", 
    "Iter_AR1" = "solid", "Direct_AR1" = "dashed",
    "Iter_AR12" = "solid", "Direct_AR12" = "dashed"
  )) +
  theme_minimal() +
  theme(legend.position = "bottom")

# 2. Handle Options (Define text and filename based on DO_U)
if (DO_U) {
  my_title <- "Unemployment Forecast Comparison"
  my_ylab  <- "Unemployment Rate (%)"
  my_file  <- "../Results/Forecast_Comparison_ggplot_U.png"
} else {
  my_title <- "Job Growth Forecast Comparison"
  my_ylab  <- "Monthly Growth Rate (%)" # Corrected from "Unemployment Rate"
  my_file  <- "../Results/Forecast_Comparison_ggplot_EMP.png"
}

# 3. Add Labels to Plot
p <- p + labs(title = my_title,
              subtitle = "Iterated vs. Direct Method (AR1 vs AR12)",
              y = my_ylab, 
              x = "")

# 4. Save
ggsave(my_file, plot = p, width = 8, height = 5)
print(p)

#######################################################
# Addendum: COVID-19 plot
######################################################



# 1. Setup Labels based on the Variable (DO_U)
if (DO_U) {
  plot_title <- "Unemployment Rate: The COVID Anomaly"
  y_label    <- "Unemployment Rate (%)"
  file_name  <- "../Results/Time_Series_UNRATE_Covid.png"
} else {
  plot_title <- "Nonfarm Payrolls (Growth): The COVID Anomaly"
  y_label    <- "Monthly Growth Rate (%)"
  file_name  <- "../Results/Time_Series_PAYEMS_Covid.png"
}

# 2. Convert ts_data to a Data Frame for ggplot
df_plot <- data.frame(
  Date = as.Date(as.yearmon(time(ts_data))),
  Value = as.numeric(ts_data)
)

# 3. Create the Plot
p_ts <- ggplot(df_plot, aes(x = Date, y = Value)) +
  
  # A. Add the Shaded Region for COVID (2020-2021)
  #    We do this FIRST so the line is drawn on top of the shading
  annotate("rect", xmin = as.Date("2020-01-01"), xmax = as.Date("2021-12-01"),
           ymin = -Inf, ymax = Inf, 
           alpha = 0.2, fill = "red") +
  
  # B. Add the Data Line
  geom_line(color = "black", size = 0.8) +
  
  # C. Add Text Annotation
  annotate("text", x = as.Date("2021-01-01"), y = max(df_plot$Value, na.rm=T), 
           label = "COVID-19\nShock", vjust = 1, color = "red", fontface = "bold") +
  
  # D. Styling
  theme_minimal() +
  labs(title = plot_title,
       subtitle = "Shaded region indicates the excluded 2020-2021 period",
       y = y_label, x = "") +
  theme(plot.title = element_text(face = "bold"))

# 4. Display and Save
print(p_ts)
ggsave(file_name, plot = p_ts, width = 8, height = 5)

print(paste("Plot saved to:", file_name))