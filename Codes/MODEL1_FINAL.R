# ==============================================================================
# Model 1: Aggregate Biophysical Growth Nexus (Stern-Based Specification)
# Logic: GDP ~ f(Capital, Energy Intensity)
# Final Diagnostic Version
# ==============================================================================

# 0. CLEAN ENVIRONMENT & LIBRARIES
# ------------------------------------------------------------------------------
rm(list = ls())
library(dplyr)
library(tidyquant) 
library(WDI)        
library(zoo)
library(urca)
library(vars)
library(car)
library(lmtest)

# 1. DATA ACQUISITION (Live Pull: 1990 - 2024)
# ------------------------------------------------------------------------------
cat("\n===== 1. PULLING LIVE MACRO DATA (1990-2024) =====\n")

gfcf_live <- WDI(country = "IN", indicator = "NE.GDI.FTOT.KD", start = 1990, end = 2024) %>%
  dplyr::rename(gfcf = NE.GDI.FTOT.KD) %>%
  dplyr::select(year, gfcf)

# 2. DATA PREPARATION & ENERGY INTENSITY
# ------------------------------------------------------------------------------
cat("\n===== 2. CONSTRUCTING BIOPHYSICAL VARIABLES =====\n")
raw_data <- read.csv("/Users/anshuman/Downloads/qog_std_ts_jan26 (1).csv")

df_model1 <- raw_data %>%
  dplyr::filter(ccodealp == "IND") %>%
  dplyr::arrange(year) %>%
  dplyr::select(year, 
                energy_imports = wdi_eneimp, 
                gdp_per_capita = wdi_gdpcapcon2015) %>%
  left_join(gfcf_live, by = "year") %>%
  filter(year >= 1990 & year <= 2024) %>%
  stats::na.omit() %>%
  mutate(
    log_gdp       = log(gdp_per_capita),
    log_gfcf      = log(gfcf),
    log_intensity = log(energy_imports / gdp_per_capita)
  )

# 3. PRE-ESTIMATION DIAGNOSTICS: MULTICOLLINEARITY (VIF)
# ------------------------------------------------------------------------------
cat("\n--- 3. PRE-ESTIMATION: VIF DIAGNOSTICS ---\n")
vif_check <- lm(log_gdp ~ log_gfcf + log_intensity, data = df_model1)
print(car::vif(vif_check))

# 4. JOHANSEN COINTEGRATION TEST
# ------------------------------------------------------------------------------
cat("\n--- 4. JOHANSEN TEST: TRACE STATISTIC ---\n")
endo_vars <- df_model1 %>% dplyr::select(log_gdp, log_gfcf, log_intensity)
endo_ts <- ts(endo_vars, start = 1990, frequency = 1)

johansen_test <- ca.jo(endo_ts, type = "trace", ecdet = "const", K = 2)
print(summary(johansen_test))

# 5. VECM ESTIMATION
# ------------------------------------------------------------------------------
cat("\n--- 5. VECM ESTIMATION (r = 1) ---\n")
vecm_model <- cajorls(johansen_test, r = 1)

# Extracting results for interpretation
cat("\n--- Long-Run Cointegrating Vector (Beta) ---\n")
print(vecm_model$beta)

cat("\n--- Short-Run Error Correction Results ---\n")
print(summary(vecm_model$rlm))

# 6. POST-ESTIMATION RESIDUAL DIAGNOSTICS (THE TRINITY)
# ------------------------------------------------------------------------------
# Convert VECM to VAR for residual testing
var_rep <- vec2var(johansen_test, r = 1)

cat("\n\n##################################################")
cat("\n FINAL RESIDUAL DIAGNOSTICS")
cat("\n##################################################\n")

# A. Serial Correlation (Portmanteau Test)
# Null: No serial correlation. We want p > 0.05
cat("\n--- A. Serial Correlation (PT Asymptotic) ---\n")
serial_test <- serial.test(var_rep, lags.pt = 10, type = "PT.asymptotic")
print(serial_test)

# B. Normality (Jarque-Bera)
# Null: Residuals are normal. We want p > 0.05
cat("\n--- B. Normality Test (Jarque-Bera) ---\n")
norm_test <- normality.test(var_rep, multivariate.only = TRUE)
print(norm_test)

# C. Heteroskedasticity (ARCH Test)
# Null: No ARCH effects. We want p > 0.05
cat("\n--- C. Heteroskedasticity (ARCH) ---\n")
arch_test <- arch.test(var_rep, lags.multi = 5, multivariate.only = TRUE)
print(arch_test)

## ==============================================================================
# MODEL 1: THE AGGREGATE VISUALIZATION SUITE
# ==============================================================================

# 1. FIGURE 1: THE COINTEGRATION "HEARTBEAT" (ECT PLOT)
# ------------------------------------------------------------------------------
# Proof that the variables are actually anchored in the long run.
beta_vector <- vecm_model$beta
data_matrix <- as.matrix(cbind(endo_ts, 1)) 
ect_series  <- ts(data_matrix %*% beta_vector, start = 1990, frequency = 1)

plot(ect_series, type = "l", col = "darkblue", lwd = 2,
     main = "Figure 1: Cointegration Relation (ECT) Over Time",
     ylab = "Equilibrium Error", xlab = "Year")
abline(h = 0, col = "red", lty = 2)
grid()

# 2. FIGURE 2: THE SHOCK-RESPONSE (IRF)
# ------------------------------------------------------------------------------
# Visualizes how a shock to Energy Intensity (Inefficiency) hurts GDP.
irf_int_gdp <- irf(var_rep, impulse = "log_intensity", response = "log_gdp", 
                   n.ahead = 8, boot = TRUE)

plot(irf_int_gdp, main = "Figure 2: Impact of Intensity Shock on GDP",
     ylab = "Log GDP Change", xlab = "Years Ahead", col = "red")

# 3. FIGURE 3: THE CONTRIBUTION STORY (FEVD)
# ------------------------------------------------------------------------------
# Shows what % of GDP volatility is driven by Energy Efficiency vs Capital.
fevd_res <- fevd(var_rep, n.ahead = 10)
plot(fevd_res, addnames = TRUE, main = "Figure 3: Variance Decomposition (FEVD)")

# 4. FIGURE 4: THE STABILITY CHECK (CUSUM)
# ------------------------------------------------------------------------------
# Proves that your coefficients didn't "break" during the 1990-2024 period.
library(strucchange)
gdp_resid <- residuals(vecm_model$rlm)[, 1]
plot(efp(gdp_resid ~ 1, type = "OLS-CUSUM"), 
     main = "Figure 4: OLS-CUSUM Stability (GDP Equation)")

