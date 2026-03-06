# ==============================================================================
# MODEL 2: INDUSTRIAL DECOUPLING (High-Precision ARDL Bounds Testing)
# Theory: Pesaran, Shin, and Smith (2001) | Cobb-Douglas Spec (Sarkar & Mathavan)
# Variables: Industrial GDP, Coal %, Capital (GFCF), Labour Force (LFPR)
# ==============================================================================

# 0. CLEAN ENVIRONMENT & LIBRARIES
# ------------------------------------------------------------------------------
rm(list = ls())
library(dplyr)
library(WDI)          # For Live Data Pulls
library(urca)         # For Zivot-Andrews & ADF tests
library(ARDL)         # For Bounds testing and ECM
library(lmtest)       # For post-estimation diagnostics
library(tseries)      # For Jarque-Bera normality test
library(strucchange)  # For CUSUM stability plots

cat("\n=================================================================")
cat("\n PART 1: DATA ACQUISITION & PREPARATION")
cat("\n=================================================================\n")

# A. Pull Capital (GFCF) and Labour (LFPR) Live from WDI
live_macro <- WDI(country = "IN", 
                  indicator = c("NE.GDI.FTOT.KD", "SL.TLF.CACT.ZS"), 
                  start = 1990, end = 2024) %>%
  dplyr::rename(gfcf_val = NE.GDI.FTOT.KD, lab_val = SL.TLF.CACT.ZS) %>%
  dplyr::select(year, gfcf_val, lab_val)

# B. Load QoG Data for GDP and Specific Coal Percentage
raw_data <- read.csv("/Users/anshuman/Downloads/qog_std_ts_jan26 (1).csv")

# C. Merge, Calculate Absolute Industrial GDP, and Transform
df_model <- raw_data %>%
  dplyr::filter(ccodealp == "IND") %>%
  dplyr::arrange(year) %>%
  dplyr::select(year, 
                gdp_total_pc  = wdi_gdpcapcon2015, 
                ind_share     = wdi_gdpind,         
                mix_coal      = wdi_elprodcoal) %>% 
  dplyr::left_join(live_macro, by = "year") %>%
  dplyr::filter(year >= 1990) %>%
  stats::na.omit() %>%
  dplyr::mutate(
    # Absolute Industrial GDP (Removes Service/Agri Sector Dilution)
    ind_gdp_abs = gdp_total_pc * (ind_share / 100),
    
    # Log transformations for Elasticity Interpretation
    log_ind_gdp = log(ind_gdp_abs),
    log_gfcf    = log(gfcf_val),
    log_coal    = log(mix_coal),
    log_lab     = log(lab_val)
  )

cat("Data Prep Complete. Sample:", min(df_model$year), "-", max(df_model$year), "| N =", nrow(df_model), "\n")

cat("\n=================================================================")
cat("\n PART 2: ENDOGENOUS STRUCTURAL BREAK DETECTION (ZIVOT-ANDREWS)")
cat("\n=================================================================\n")

# Automatically detect the true structural break year in Industrial GDP
za_test <- urca::ur.za(df_model$log_ind_gdp, model = "both")
break_index <- za_test@bpoint
break_year <- df_model$year[break_index]

cat("\n--- Zivot-Andrews Test Detected Structural Break in Year:", break_year, "---\n")

# Create the dynamic dummy variable based on the scientific break year
df_model <- df_model %>%
  dplyr::mutate(d_break = ifelse(year >= break_year, 1, 0))

cat("\n=================================================================")
cat("\n PART 3: ARDL MODEL ESTIMATION & BOUNDS TEST")
cat("\n=================================================================\n")

# 1. OPTIMAL ARDL LAG SELECTION (AIC)
cat("\n===== 1. LAG SELECTION =====\n")
# Spec includes Labour and the dynamic Break Dummy
ardl_models <- ARDL::auto_ardl(log_ind_gdp ~ log_gfcf + log_coal + log_lab + d_break, 
                               data = df_model, max_order = 2)

best_ardl <- ardl_models$best_model
cat("Optimal ARDL Order Found:", ardl_models$best_order, "\n")

# 2. THE BOUNDS TEST (Testing for Cointegration / Lock-in)
cat("\n===== 2. PESARAN BOUNDS TEST FOR COINTEGRATION =====\n")
# Case 3: Unrestricted intercept, no trend
bounds_res <- ARDL::bounds_f_test(best_ardl, case = 3)
print(bounds_res)

# 3. UNRESTRICTED ECM (Short-Run Dynamics)
cat("\n===== 3. SHORT-RUN DYNAMICS & ERROR CORRECTION TERM (ECT) =====\n")
uecm_model <- ARDL::uecm(best_ardl)
summary_uecm <- summary(uecm_model)

cat("\n--- Error Correction Term (Speed of Adjustment) ---\n")
print(summary_uecm$coefficients["L(log_ind_gdp, 1)", c("Estimate", "Std. Error", "t value", "Pr(>|t|)")])

cat("\n--- Full Unrestricted ECM Equation ---\n")
print(summary_uecm$coefficients)

# 4. LONG-RUN ELASTICITIES
cat("\n===== 4. LONG-RUN COEFFICIENTS (BETA) =====\n")
lr_multipliers <- ARDL::multipliers(best_ardl)
print(lr_multipliers)

cat("\n=================================================================")
cat("\n PART 4: POST-ESTIMATION DIAGNOSTICS & STABILITY (THE TRINITY)")
cat("\n=================================================================\n")

# A. Serial Correlation (Breusch-Godfrey Test)
cat("\n--- 1. Serial Correlation (Breusch-Godfrey) ---\n")
print(lmtest::bgtest(best_ardl, order = 1))

# B. Heteroskedasticity (Breusch-Pagan Test)
cat("\n--- 2. Heteroskedasticity (Breusch-Pagan) ---\n")
print(lmtest::bptest(best_ardl))

# C. Normality of Residuals (Jarque-Bera Test)
cat("\n--- 3. Residual Normality (Jarque-Bera) ---\n")
print(tseries::jarque.bera.test(residuals(best_ardl)))

# D. Functional Form / Specification (Ramsey RESET Test)
cat("\n--- 4. Model Specification (Ramsey RESET) ---\n")
print(lmtest::resettest(best_ardl, power = 2, type = "fitted"))

# ==============================================================================
# PART 5: PUBLICATION-QUALITY VISUALIZATIONS DASHBOARD
# ==============================================================================
cat("\n===== GENERATING DIAGNOSTIC PLOTS DASHBOARD =====\n")

# Set up a 2x2 plotting grid
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0), bg = "white")

# ------------------------------------------------------------------------------
# Plot 1: Actual vs. Fitted Values
# ------------------------------------------------------------------------------
fitted_vals <- as.numeric(fitted(best_ardl))
actual_vals <- tail(df_model$log_ind_gdp, length(fitted_vals))
time_index  <- seq_along(actual_vals)

plot(time_index, actual_vals, type = "l", col = "black", lwd = 2, 
     ylab = "Log Industrial GDP", xlab = "Time (Index)",
     main = "Graph A: Actual vs. Fitted Values")
lines(time_index, fitted_vals, col = "blue", lwd = 2, lty = 2)
legend("topleft", legend = c("Actual", "Fitted"), col = c("black", "blue"), 
       lty = c(1, 2), lwd = 2, bty = "n")

# ------------------------------------------------------------------------------
# Plot 2: Residuals Plot with 95% Confidence Bands
# ------------------------------------------------------------------------------
res <- as.numeric(residuals(best_ardl))
sigma <- sd(res)
max_bound <- max(abs(res), 2.5 * sigma) * 1.1 

plot(time_index, res, type = "h", col = "darkred", lwd = 2,
     ylab = "Residuals", xlab = "Time (Index)",
     main = "Graph B: Model Residuals (±2σ Bands)",
     ylim = c(-max_bound, max_bound))
abline(h = 0, col = "black", lwd = 2)
abline(h = 2 * sigma, col = "blue", lty = 2, lwd = 1.5)
abline(h = -2 * sigma, col = "blue", lty = 2, lwd = 1.5)



# ------------------------------------------------------------------------------
# Plot 3 & 4: THE CUSUM FIX (Bypassing ARDL formula parser errors)
# ------------------------------------------------------------------------------
df_plot <- best_ardl$model
colnames(df_plot) <- make.names(colnames(df_plot)) # Safe column names
cusum_form <- as.formula(paste(colnames(df_plot)[1], "~ ."))

# Plot 3: CUSUM Test (Parameter Stability)
cusum_test <- strucchange::efp(cusum_form, data = df_plot, type = "Rec-CUSUM")
plot(cusum_test, main = "Graph C: CUSUM Plot (Parameter Stability)", 
     col = "blue", lwd = 2, ylab = "Empirical Fluctuation")

# Plot 4: CUSUM of Squares Test (Variance Stability)
cusum_sq_test <- strucchange::efp(cusum_form, data = df_plot, type = "Rec-MOSUM")
plot(cusum_sq_test, main = "Graph D: CUSUM of Squares (Variance)", 
     col = "darkgreen", lwd = 2)

mtext("ARDL Model 2 Diagnostics: The Industrial Core", outer = TRUE, cex = 1.5, font = 2)
par(mfrow = c(1, 1))

cat("\n===== FULL ARDL PIPELINE COMPLETE =====\n")

