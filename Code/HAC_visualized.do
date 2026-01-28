clear all
set seed 027
local runs 2000
local N 1000
local p 0.2

* Calculate dynamic lag length based on S&W rule-of-thumb
local lags = ceiling(.75* (`N')^(1/3))
display "Using `lags' lags for HAC adjustment"

frame create results p_ols p_hac

forvalues i = 1/`runs' {
    quietly {
        clear
        set obs `N'
        gen t = _n
        gen x = rnormal()
        
        * Generate "Sticky" Data
        forval j = 2/`N' {
            replace x = x[`j'-1] in `j' if runiform() < `p'
        }
        
        tsset t
        
        * 1. Standard OLS
        reg x
        local p_ols = 2*ttail(e(df_r), abs(_b[_cons]/_se[_cons]))
        
        * 2. Newey-West (HAC) with Dynamic Lags
        newey x, lag(`lags')
        local p_hac = 2*ttail(e(df_r), abs(_b[_cons]/_se[_cons]))
        
        frame post results (`p_ols') (`p_hac')
    }
}

cwf results
* Switch to results frame
cwf results

* Calculate Rejection Rates
gen rej_ols = p_ols < .05
gen rej_hac = p_hac < .05

di "Type I Error Rate (Standard OLS):"
summarize rej_ols
di "Type I Error Rate (HAC):"
summarize rej_hac

* Plot 1: Standard OLS (Should show a huge spike near 0)
hist p_ols, width(0.05) start(0) freq ///
    title("OLS p-values") ///
    subtitle("Highly biased toward False Significance") ///
    xtitle("p-value") name(ols_hist, replace)

* Plot 2: HAC (Should look roughly flat/uniform)
hist p_hac, width(0.05) start(0) freq ///
    title("HAC p-values") ///
    subtitle("Correctly distributed under the Null") ///
    xtitle("p-value") name(hac_hac, replace)

* Combine them for the slide
graph combine ols_hist hac_hac, cols(2) xsize(10) ysize(5)