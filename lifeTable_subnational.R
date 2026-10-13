##################
### Life Table ###
##################

# For now it proves as a robust and faster version for life table,
# reducing for about 3/4 time.
# It use data.table for boosting process, suggest using parallel computation.

library(data.table)
library(dplyr)

leTable_dt <- function(data_dt, year_val, sex_val, area_val){
  if(!is.data.table(data_dt)){
    data_dt <- as.data.table(data_dt)
  }
  # Filter
  filtered_dt <- data_dt[year == year_val & sex == sex_val & area == area_val,
                         .(sex, year, area, age, dths, pop)]
  return(filtered_dt)
}

### Kannisto-Thatcher life table extension in high age ###
ktExtension <- function(lx70, return_e85 = FALSE){
  # make sure input of lx70 is a vector
  # For age groups >= 70, calculate hazard rate using the approximation
  # mu(x+1/2) ~ -log(1 - q) = -log(p) where p is the probability of
  # survival to the next age group and equals l(x+5) / lx
  
  mux <- (log(lx70) - log(c(lx70[-1], NA))) / 5
  mux <- case_when(is.na(mux) ~ 0, T ~ mux)
  
  valid_mux_indices <- !is.na(mux)
  mux_valid <- mux[valid_mux_indices]
  lx70_valid_for_mux <- lx70[valid_mux_indices]
  
  # # Check if the data is enable to do following steps
  if(length(mux_valid) < 3){
    warning("Insufficient valid transition rates (mux) calculated from lx70. Need at least 3 for intervals 70-74, 75-79, 80-84. Returning default ax values (2.5).")
    return(if (return_e85) NA_real_ else rep(2.5, 4))
  }
  # Calculate lx for 1-year age groups from 70 to 84 (or up to available data)
  # The original code assumes 4 groups (70-74, 75-79, 80-84, 85+)
  # Let's replicate the structure assuming 4 groups input
  # Rep lengths depend on having 3 intervals for mux (70-74 -> 75-79, etc.)
  
  # Generate AGP sequence
  agp_seq <- 70:84 # Original code goes up to 85 for calculation setup
  
  # Create the data frame/data table structure
  ktData_list <- list(
    AGP = agp_seq,
    lx = lx70[rep(1: 3, each = 5)], # Replicate lx for each 5-year group up to 84
    int = rep(0:4, 3),             # Interval within the 5-year group
    mux = rep(mux_valid[1:3], each = 5) # Replicate mux
  )
  
  # Add the starting point for age 85+
  ktData_list$AGP <- c(ktData_list$AGP, 85)
  ktData_list$lx <- c(ktData_list$lx, lx70[4]) # lx at start of 85+ group
  ktData_list$int <- c(ktData_list$int, 0)     # Interval is 0 for the start of 85+
  ktData_list$mux <- c(ktData_list$mux, NA)    # Mux is not applicable here
  
  ktData <- as.data.table(ktData_list) # Convert to data.table
  
  # Calculate single-year lx up to 84, and qx / logitqx
  ktData[, lx := ifelse(AGP == 85, lx, lx * exp(mux * -int))] # Calculate lx within 5-year groups
  ktData[, dx := lx - shift(lx, type = "lead")] # dx = lx - l(x+1)
  ktData[AGP == 85, dx := NA] # dx for the last age (85) needs the next lx (from model)
  ktData[, qx := dx / lx]
  ktData[, logitqx := log(qx / (1 - qx))]
  ktData <- ktData[AGP < 85] # Keep only data up to age 84 for regression
  
  # Run regression on logit of probability of dying by logit
  y <- ktData$logitqx
  # x <- ktData$AGP + 0.5 # Age at midpoint
  keep <- seq(1, nrow(ktData), by = 5)   # 70-74, 75-79, 80-84 各取第一列
  y <- ktData$logitqx[keep]
  x <- ktData$AGP[keep] + 2.5
  
  # Check for NA/Inf values before regression
  valid_reg <- !is.na(y) & is.finite(y) & !is.na(x) & is.finite(x)
  if(sum(valid_reg) < 2){
    warning("Not enough valid data points (need at least 2) for logit qx regression for ages 70-84 in ktExtension. Unable to extrapolate mortality rates beyond 85. Returning default ax values (2.5) for all groups >= 70.")
    return(if (return_e85) NA_real_ else rep(2.5, 4)) # Return ax for 70-74, 75-79, 80-84, 85+
  }
  
  # Used weights = x
  mod <- lm(y[valid_reg] ~ x[valid_reg], weights = x[valid_reg])
  logA <- coef(mod)[1]
  B <- coef(mod)[2]
  
  # Calculate qx for age x >= 85
  ages_extrapolate <- 85: 129
  logitqx85 <- logA + B * (ages_extrapolate + 0.5)
  qx85 <- exp(logitqx85) / (1 + exp(logitqx85))
  
  # Prepare data table for ages 85+
  ktData85 <- data.table(
    AGP = ages_extrapolate,
    lx = NA_real_, # Initialize lx column
    qx = qx85
  )
  ktData85[1, lx := lx70[4]] # Starting lx for age 85
  
  # Calculate lx and dx for age x >= 85 using data.table efficiently
  # Calculate px first
  ktData85[, px := 1 - qx]
  
  # Calculate lx by iteration
  for(i in 2: nrow(ktData85)){
    ktData85$lx[i] <- ktData85$lx[i - 1] * ktData85$px[i - 1]
  }
  # Calculate dx
  ktData85[, dx := lx * qx]

  # Person-years above 85 from the extrapolated single-year table (ax = 0.5)
  e85 <- sum(ktData85$lx - ktData85$dx / 2) / lx70[4]
  if(return_e85) return(e85)

  # Combine data < 85 and >= 85
  # We only need the calculated dx for the final ax calculation
  # The structure needs dx for ages 70-129
  dx_full <- c(ktData$dx, ktData85$dx) # Combine dx values
  agp_full <- 70: 129
  
  # Calculate ax for the original 5-year age groups (70-74, 75-79, 80-84, 85+)
  # Define the original age groups for aggregation
  agp_groups <- c(70, 75, 80, 85)
  # Assign each single year age (70-129) to its corresponding group start age
  group_map <- agp_groups[findInterval(agp_full, agp_groups)]
  
  # Create a temporary data table for aggregation
  temp_ax_calc <- data.table(agp_full, dx_full, group_map)
  
  # Midpoint calculation adjustment: yl = age - group_start + 0.5
  temp_ax_calc[, yl := agp_full - group_map + 0.5]
  
  # Calculate sum(dx * yl) and sum(dx) per group
  ax_calc <- temp_ax_calc[, .(sum_dx_yl = sum(dx_full * yl, na.rm = TRUE),
                              sum_dx = sum(dx_full, na.rm = TRUE)), by = group_map]
  
  # Calculate ax = sum(dx * yl) / sum(dx)
  ax_calc[, ax := ifelse(sum_dx == 0, 2.5, sum_dx_yl / sum_dx)] # Use 2.5 if sum(dx) is 0 (or handle as needed)
  
  # Return the ax values for the 70-74, 75-79, 80-84, 85+ groups
  # Ensure the order matches the expected output (corresponding to lx70 input)
  # The groups correspond to indices 16, 17, 18, 19 in the main leCalc data
  # Need to return ax for 70-74, 75-79, 80-84, 85+ (4 values)
  final_ax <- ax_calc[order(group_map)]$ax
  
  # The result should be 4 values
  if(length(final_ax) != 4){
    warning("ktExtension did not produce exactly 4 ax values. Check input lx70 length and logic.")
    return(if (return_e85) NA_real_ else rep(2.5, 4))
  }
  
  return(final_ax)
}


### Young age ax calculation ###
# This function is simple, no significant optimization needed
youngExtension_JPN <- function(mx5, sex){
  # mx5[1] = m(0–4) or proxy of infant mortality
  m0 <- mx5[1]
  
  if(sex == 1){  # Male
    a0   <- ifelse(m0 >= 0.107, 0.330, 0.045 + 2.684 * m0)
    a1_4 <- ifelse(m0 >= 0.107, 1.352, 1.651 - 2.816 * m0)
  }else if(sex == 2){  # Female
    a0   <- ifelse(m0 >= 0.107, 0.350, 0.053 + 2.800 * m0)
    a1_4 <- ifelse(m0 >= 0.107, 1.361, 1.522 - 1.518 * m0)
  }else{
    stop("sex must be 1 (male) or 2 (female)")
  }

  a0_4 <- (1 * a0 + 4 * a1_4) / 5
  return(a0_4)
}


### LE calculation using data.table ###
leCalc_dt_70q0 <- function(data_dt){
  
  setorder(data_dt, age)
  
  
  
  sex_val <- data_dt$sex[1]
  n_rows <- nrow(data_dt)
  
  # Define nx values directly
  nx_vals <- c(rep(5, n_rows - 1), NA)
  
  data_dt[, id := .I] # Row number
  data_dt[, ax := 2.5]
  data_dt[, nx := nx_vals]
  
  # Calculate mx, handle POP < 1 or mx == 0
  data_dt[, mx := ifelse(pop < 1, 1, dths / pop)]
  data_dt[mx == 0, mx := 1e-5] # Handle zero mortality case
  
  # Calculate initial qx, px, lx, dx
  data_dt[, qx := ifelse(age == 85, 1, (nx * mx) / (1 + (nx - ax) * mx))]
  data_dt[qx > 1, qx := 0.99999]
  
  px_shifted <- shift(1 - data_dt$qx, n=1, fill=1)
  data_dt[, lx := 1e5 * cumprod(px_shifted)]
  data_dt[, dx := lx * qx]
  
  lx70_vals <- data_dt[age >= 70, lx]
  
  mx5_vals <- data_dt[age == 0, mx] 
  if(length(mx5_vals) >= 1){
    data_dt[age == 0, ax := youngExtension_JPN(mx5_vals, sex_val)]
  }
  
  lx70_vals <- data_dt[age >= 70, lx]
  if(length(lx70_vals) >= 4){
    data_dt[age %in% c(70, 75, 80, 85), ax := ktExtension(lx70_vals)]
  }
  
  
  
  ### Iteration using data.table ###
  iter <- 4 # Number of iterations

  for(i in seq(iter)){
    dx_old <- data_dt$dx
    ax_new <- data_dt$ax

    for(k in 2:(n_rows - 1)){ 

      if(data_dt$age[k] < 70){
        dx_km1 <- dx_old[k - 1]
        dx_k   <- dx_old[k]
        dx_kp1 <- dx_old[k + 1]
        
        if(dx_k != 0){

          ax_new[k] <- (-5 / 24 * dx_km1 + 2.5 * dx_k + 5 / 24 * dx_kp1) / dx_k
          ax_new[k] <- pmax(pmin(ax_new[k], 5), 0)
        }
      }
    }
    data_dt[, ax := ax_new]
    
    data_dt[, qx := ifelse(age == 85, 1, (nx * mx) / (1 + (nx - ax) * mx))]
    data_dt[qx > 1, qx := 0.99999]
    data_dt[, px := 1 - qx]
    data_dt[, lx := 1e5 * cumprod(shift(px, n=1, fill=1))]
    data_dt[, dx := lx * qx]
  }
  

  # Close the 85+ interval with Kannisto extrapolation (fallback: lx / mx)
  e85 <- ktExtension(data_dt[age >= 70, lx], return_e85 = TRUE)

  data_dt[, lx_lead := shift(lx, type = "lead")]
  data_dt[, Lx := nx * lx_lead + ax * dx]
  data_dt[id == n_rows, Lx := if (is.na(e85)) lx / mx else lx * e85] # 封閉年齡組 85+
  
  data_dt[, Tx := rev(cumsum(rev(Lx)))]
  data_dt[, ex := ifelse(lx == 0, 0, Tx / lx)]
  # Calculate Standard Error and Confidence Intervals
  # spx = variance of qx
  # W_spx = contribution to variance of ex[0]
  # STx = variance of Tx
  # SeSE = standard error of ex
  
  data_dt[, spx := fifelse(id == n_rows,
                           4 / (dths * (mx ^ 2)), # Approximation for last interval variance
                           fifelse(dths == 0, 0, (qx ^ 2) * (1 - qx) / dths))]
  
  data_dt[, spx := fifelse(id == n_rows,
                           fifelse(dths == 0 | mx == 0, 0, 4 / (dths * (mx^2))), # Handle COUNT=0 or mx=0
                           fifelse(dths == 0, 0, (qx^2) * (1 - qx) / dths))]
  
  # Need ex_lead (ex[x+1]) for W_spx calculation
  data_dt[, ex_lead := shift(ex, type = "lead", fill = 0)]
  
  data_dt[, W_spx := fifelse(id < n_rows,
                             spx * (lx^2) * (((1 - ax) * nx + ex_lead)^2),
                             ((lx / 2)^2) * spx)]
  
  data_dt[, STx := rev(cumsum(rev(W_spx)))] # Variance of Tx (cumulative from bottom)
  data_dt[, SeSE := fifelse(lx == 0, 0, sqrt(STx / (lx^2)))] # Standard Error of ex
  
  # Confidence intervals
  z <- 1.96 # For 95% CI (Original code had z=0.95, which is likely wrong, z=1.96 for 95%)
  # Let's use 1.96 assuming 95% CI is desired.
  
  z <- qnorm(1 - (1 - 0.95) / 2) # Standard way to get z for 95% CI
  
  data_dt[, lowercl := ex - z * SeSE]
  data_dt[, uppercl := ex + z * SeSE]
  
  # Ensure lower CI is not negative
  data_dt[lowercl < 0, lowercl := 0]
  
  # Prepare output table
  # Add age group labels
  agp_labels <- c("0-4", paste0(seq(5, 80, by = 5), "-", seq(9, 84, by = 5)), "85+")
  data_dt[, agp := agp_labels]
  # Select and order columns for the final output
  output_dt <- data_dt[, .(agp, age, deaths=dths, pops=pop, nx, ax, mx, qx, dx, lx, Lx, Tx, ex, lowercl, uppercl)]
  
  return(output_dt)
}


### LE compute by loop wrapper (using data.table functions) ###
leCompute_dt_70q0 <- function(data_dt, year_val, sex_val, area_val){
  # Ensure input is data.table
  if(!is.data.table(data_dt)){
    data_dt <- as.data.table(data_dt)
  }
  lt_dt <- leTable_dt(data_dt, year_val, sex_val, area_val)
  
  # Check if filtering resulted in data
  if(nrow(lt_dt) == 0){
    warning(paste("No data found for year =", year_val, "sex =", sex_val, "area =", area_val))
    return(NULL) # Return NULL or an empty data.table structure
  }
  # Check if lt_dt has the expected 19 rows
  if(nrow(lt_dt) != 18){
    warning(paste("Data for year =", year_val, "sex =", sex_val, "area =", area_val,
                  "does not have the expected 18 age groups. Skipping calculation."))
    return(NULL)
  }
  
  result_dt <- leCalc_dt_70q0(lt_dt)
  result_dt[, year := year_val]
  result_dt[, sex := sex_val]
  result_dt[, area := area_val]
  
  output_dt <- result_dt[, .(year, sex, area, agp, age, ax, dx, mx, lx, Lx, Tx, ex, lowercl, uppercl)]
  
  return(output_dt)
}

