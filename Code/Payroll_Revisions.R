# Gemini prompt:
# I want to analyze the most recent employment report in the U.S. using R. #
# Specifically, I am interested in using R and FRED to visualize the revisions 
# in nonfarm employment.
# Write simple code that loads the relevant libraries for downloading PAYEMS 
# with two custom vintage dates. Then, plot the two vintage dates using #
# base R and output it to an PNG file.
#
# What's missing: specifics about the vintage dates and plots.
# There's also no quality check, but this is simple enough for me to follow along

# 1. Load necessary libraries
# install.packages("fredr") # Uncomment if not installed
library(fredr)

# 2. Setup FRED API Key
# Replace 'YOUR_API_KEY_HERE' with your actual 32-character FRED API key
# Get your FRED API key here
# https://fred.stlouisfed.org/docs/api/api_key.html
fredr_set_key("your API key here")


# 3. Define Parameters
# PAYEMS = Total Nonfarm | USPRIV = Total Private
series_id <- "PAYEMS" 

# Set two vintage dates to compare revisions
# Example: Comparing data as known in early Jan 2026 vs. mid-Feb 2026
# Wrap the strings in as.Date()
vintage_date_1 <- as.Date("2026-01-01")
vintage_date_2 <- as.Date("2026-02-15")

# 4. Download Data Vintages
# Fetch data as it existed on Vintage Date 1
df_v1 <- fredr(
  series_id = series_id,
  realtime_start = vintage_date_1,
  realtime_end   = vintage_date_1
)

# Fetch data as it existed on Vintage Date 2
df_v2 <- fredr(
  series_id = series_id,
  realtime_start = vintage_date_2,
  realtime_end   = vintage_date_2
)

# Convert to millions
df_v1$value <- df_v1$value/1000
df_v2$value <- df_v2$value/1000

DO_NORMALIZE <- FALSE

# Normalize 2024m1 to 0?
if (DO_NORMALIZE) {
  df_v1$value <- df_v1$value-df_v1$value[df_v1$date==as.Date("2024-01-01")]
  df_v2$value <- df_v2$value-df_v2$value[df_v2$date==as.Date("2024-01-01")]
  my_label <- "Millions of Persons ['24m1=0]"
} else {
  my_label <- "Millions of Persons"
}

# 5. Filter Data for Clarity
# Subset to the last 24 months so revisions are visible in the plot
start_filter <- as.Date("2024-01-01")
  #as.Date(Sys.Date()) - 730*1.5 # Approx 2 years
df_v1_subset <- subset(df_v1, date >= start_filter)
df_v2_subset <- subset(df_v2, date >= start_filter)

# 6. Plot and Save to PNG
# Plot the first vintage (Older data)
my_plot <- plot(df_v1_subset$date, df_v1_subset$value,
     type = "l", col = "blue", lwd = 3,
     main = paste("Employment Real Time Estimates:", series_id),
     xlab = "Date", ylab = my_label,
     las = 1) # las=1 makes y-axis labels horizontal

# Add the second vintage (Newer data)
lines(df_v2_subset$date, df_v2_subset$value, 
      col = "red", lwd = 3, lty = 2)

# Add Grid and Legend
grid()
legend("topleft", 
       legend = c(paste("Vintage:", vintage_date_1), 
                  paste("Vintage:", vintage_date_2)),
       col = c("blue", "red"), 
       lty = c(1, 2), lwd = 2)

# 2. Copy the active plot (screen) to a PNG file
dev.copy(png, filename = "../Results/payems_revisions.png", width = 800, height = 600)

# 3. Close the file connection
dev.off()

message("Plot saved as 'payems_revisions.png'")
 