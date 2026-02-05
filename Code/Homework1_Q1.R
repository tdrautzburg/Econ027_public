# List of required packages
pkgs <- c("fredr", "dplyr", "purrr", "knitr")

# Install missing ones
installed_pkgs <- rownames(installed.packages())
for (p in pkgs) {
  if (!(p %in% installed_pkgs)) install.packages(p)
}

# Load them
library(fredr)
library(dplyr)
library(purrr)
library(knitr)

#################
# 0
# Housekeeping
#################

library(here)

# Clear the Environment/console/plots
rm(list = ls())
cat("\014")
if(!is.null(dev.list())) dev.off()

# Set project directory
i_am("Homework1_Q1.R")

#####################
# 1. Download data
#####################

# 1. Set your API key
# https://fred.stlouisfed.org/docs/api/api_key.html
fredr_set_key("YOUR KEY HERE")


# Define the indicators and their labels
indicators <- c(
  "PAYEMS"    = "Employment Growth (thousands, avg 3-mo chg)",
  "UNRATE"    = "Unemployment Rate (%)",
  "A191RL1Q225SBEA" = "GDP Growth (Quarterly %)",
  "PCEC96"       = "Real PCE (Annualized MoM)",
  "PCEPILFE"  = "Core PCE Inflation (Annualized MoM)",
  "PCEPI"     = "Headline PCE Inflation (Annualized MoM)",
  "CSUSHPISA"   = "House Price Growth (Annualized MoM)",
  "FEDFUNDS"  = "Federal Funds Rate (%)",
  "MORTGAGE30US" = "30-Year Mortgage Rate (%)",
  "SP500"  = "SP 500 Index (MoM %)"
)

# Function to fetch and transform data
get_economic_data <- function(id, label) {
  df <- fredr(series_id = id, observation_start = as.Date("2022-01-01"))
  
  df <- df %>%
    arrange(date) %>%
    mutate(
      transformed_val = case_when(
        # 1) Quarterly change in PAYEMS
        series_id == "PAYEMS" ~ (value - lag(value, 3))/3,
        
        # 2) Annualized monthly inflation and house price growth
        series_id %in% c("PCEPILFE", "PCEPI", "CSUSHPISA", "SP500", "PCEC96") ~ ((value / lag(value, 1))^12 - 1) * 100,
        
        # 3) Raw stock index growth
        series_id %in% c( "SP500") ~ ((value / lag(value, 1)) - 1) * 100,

                # Others remain as raw values (rates or levels)
        TRUE ~ value
      ),
      Indicator = label
    ) %>%
    filter(!is.na(transformed_val))
  
  return(df)
}

#####################
# 2. Make dashoboard table
#####################
# Map over the indicators to create a combined dataset
all_data <- imap_dfr(indicators, ~ get_economic_data(.y, .x))

# 1. Final processing of the summary table
summary_table <- all_data %>%
  group_by(Indicator) %>%
  summarise(
    Latest_Date = last(date),
    Latest_Value = last(transformed_val),
    Value_12m_Ago = transformed_val[which.min(abs(date - (last(date) - 365)))]
  ) %>%
  mutate(Difference = Latest_Value - Value_12m_Ago) %>%
  
  # 2. Automatically use the order from your 'indicators' vector
  mutate(Indicator = factor(Indicator, levels = unname(indicators))) %>%
  arrange(Indicator)

# 3. Output
if (requireNamespace("knitr", quietly = TRUE)) {
  print(knitr::kable(summary_table, digits = 2))
} else {
  print(as.data.frame(summary_table))
}

# Generate TeX output
# Note: 'booktabs = TRUE' requires \usepackage{booktabs} in your LaTeX preamble
tex_output <- knitr::kable(
  summary_table, 
  format = "latex", 
  booktabs = TRUE, 
  digits = 2,
  caption = "Summary of Economic Indicators",
  col.names = c("Indicator", "Latest Date", "Current Value", "12m Ago", "Diff")
)

# Print to console so you can copy/paste
cat(tex_output)

writeLines(tex_output, "../Results/Homework/HW1_Q1_dashboard.tex")