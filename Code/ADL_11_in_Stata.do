// --- 1. Setup and Data Generation ---
clear
set obs 100
set seed 027

// Generate a time index (Required for time-series operators)
gen t = _n
tsset t

// Define parameters
scalar phi   = 0.7     // Persistence in y
scalar beta0 = 1.5     // Contemporaneous effect of x
scalar beta1 = -0.8    // Lagged effect of x

// Generate x and u
gen x = rnormal()
gen u = rnormal()

// Initialize y
gen y = .
replace y = (beta0 * x + u) / (1 - phi) in 1

// Simulation loop (Recursive)
forval i = 2/100 {
    replace y = phi*y[`i'-1] + beta0*x[`i'] + beta1*x[`i'-1] + u[`i'] in `i'
}

// --- 2. Manual OLS ---
// In Stata, "Manual" still uses the lag operators because 
// they are safer than manual indexing.
regress y L1.y x L1.x

// --- 3. Dynamic Estimation ---
// Stata's 'regress' command with factor variables IS the 
// equivalent of dynlm(). It handles the loss of the first 
// observation automatically.
regress y L.y x L.x

// --- 4. Post-Estimation Comparison ---
// View coefficients to verify they match your R results
matrix list e(b)