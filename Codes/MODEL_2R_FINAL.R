# ==============================================================================
# MODEL 2: INDUSTRIAL DECOUPLING ANALYSIS
# ARDL WITH ENERGY INTENSITY + COAL
# Theory: Pesaran, Shin & Smith (2001)
# ==============================================================================


# ==============================================================================
# 0. CLEAN ENVIRONMENT
# ==============================================================================

rm(list = ls())

library(dplyr)
library(WDI)
library(urca)
library(ARDL)
library(lmtest)
library(tseries)
library(strucchange)
library(ggplot2)
library(scales)
library(car)



# ==============================================================================
# PART 1: DATA ACQUISITION
# ==============================================================================

cat("\n================================================")
cat("\n PART 1: DATA ACQUISITION")
cat("\n================================================\n")


live_macro <- WDI(
  country = "IN",
  indicator = c(
    "NE.GDI.FTOT.KD",
    "SL.TLF.CACT.ZS"
  ),
  start = 1990,
  end = 2024
) %>%
  rename(
    gfcf_val = NE.GDI.FTOT.KD,
    lab_val  = SL.TLF.CACT.ZS
  ) %>%
  select(year, gfcf_val, lab_val)



raw_data <- read.csv("/Users/anshuman/Downloads/qog_std_ts_jan26 (1).csv")



# ==============================================================================
# PART 2: DATA PREPARATION
# ==============================================================================

cat("\n================================================")
cat("\n PART 2: DATA PREPARATION")
cat("\n================================================\n")


df_model <- raw_data %>%
  
  filter(ccodealp == "IND") %>%
  
  arrange(year) %>%
  
  select(
    year,
    gdp_total_pc = wdi_gdpcapcon2015,
    ind_share = wdi_gdpind,
    energy_imports = wdi_eneimp,
    coal_share = wdi_elprodcoal
  ) %>%
  
  left_join(live_macro, by = "year") %>%
  
  filter(year >= 1990) %>%
  
  na.omit() %>%
  
  mutate(
    
    ind_gdp_abs = gdp_total_pc * (ind_share / 100),
    
    energy_intensity = energy_imports / gdp_total_pc,
    
    log_ind_gdp = log(ind_gdp_abs),
    log_gfcf = log(gfcf_val),
    log_intensity = log(energy_intensity),
    log_coal = log(coal_share),
    log_lab = log(lab_val)
    
  )


cat("Sample:", min(df_model$year), "-", max(df_model$year),
    "| N =", nrow(df_model), "\n")



# ==============================================================================
# PART 3: STRUCTURAL BREAK DETECTION
# ==============================================================================

cat("\n================================================")
cat("\n PART 3: STRUCTURAL BREAK TEST")
cat("\n================================================\n")


za_test <- urca::ur.za(df_model$log_ind_gdp, model = "both")

break_index <- za_test@bpoint

break_year  <- df_model$year[break_index]

cat("Zivot-Andrews Break Year:", break_year, "\n")


df_model <- df_model %>%
  mutate(d_break = ifelse(year >= break_year, 1, 0))



# ==============================================================================
# PART 4: ARDL MODEL ESTIMATION
# ==============================================================================

cat("\n================================================")
cat("\n PART 4: ARDL MODEL")
cat("\n================================================\n")


ardl_models <- auto_ardl(
  
  log_ind_gdp ~ log_gfcf + log_intensity + log_coal + log_lab + d_break,
  
  data = df_model,
  
  max_order = 1,
  
  selection = "BIC"
  
)

best_ardl <- ardl_models$best_model

cat("Optimal Lag Structure:", ardl_models$best_order, "\n")



# ==============================================================================
# PART 5: COINTEGRATION TEST
# ==============================================================================

cat("\n================================================")
cat("\n PART 5: BOUNDS TEST")
cat("\n================================================\n")


bounds_res <- bounds_f_test(best_ardl, case = 3)

print(bounds_res)



# ==============================================================================
# PART 6: SHORT RUN DYNAMICS
# ==============================================================================

cat("\n================================================")
cat("\n PART 6: SHORT RUN DYNAMICS")
cat("\n================================================\n")


uecm_model <- uecm(best_ardl)

summary_uecm <- summary(uecm_model)

print(summary_uecm$coefficients)



# ==============================================================================
# PART 7: LONG RUN ELASTICITIES
# ==============================================================================

cat("\n================================================")
cat("\n PART 7: LONG RUN ELASTICITIES")
cat("\n================================================\n")


lr_multipliers <- multipliers(best_ardl)

print(lr_multipliers)



# ==============================================================================
# PART 8: DIAGNOSTIC TESTS
# ==============================================================================

cat("\n================================================")
cat("\n PART 8: DIAGNOSTIC TESTS")
cat("\n================================================\n")


bgtest(best_ardl, order = 1)

bptest(best_ardl)

jarque.bera.test(residuals(best_ardl))

resettest(best_ardl, power = 2, type = "fitted")



# ==============================================================================
# PART 9: MULTICOLLINEARITY TEST
# ==============================================================================

cat("\n================================================")
cat("\n PART 9: VIF TEST")
cat("\n================================================\n")


vif_model <- lm(
  log_ind_gdp ~ log_gfcf + log_intensity + log_coal + log_lab + d_break,
  data = df_model
)

print(car::vif(vif_model))



# ==============================================================================
# PART 10: MODEL STABILITY
# ==============================================================================

cat("\n================================================")
cat("\n PART 10: MODEL STABILITY")
cat("\n================================================\n")


par(mfrow=c(2,2))


plot(df_model$year,
     df_model$log_ind_gdp,
     type="l",
     col="black",
     lwd=2,
     main="Actual vs Fitted")

lines(df_model$year,
      fitted(best_ardl),
      col="blue",
      lty=2)


plot(df_model$year,
     residuals(best_ardl),
     type="h",
     col="darkred")

abline(h=0)


cusum_test <- efp(log_ind_gdp ~ log_gfcf + log_intensity + log_coal + log_lab + d_break,
                  data = df_model,
                  type = "Rec-CUSUM")

plot(cusum_test)

par(mfrow=c(1,1))



# ==============================================================================
# PART 11: BAI-PERRON BREAK TEST
# ==============================================================================

cat("\n================================================")
cat("\n PART 11: BAI-PERRON BREAK TEST")
cat("\n================================================\n")


bp_model <- breakpoints(log_ind_gdp ~ 1, data=df_model)

summary(bp_model)

bp_indices <- bp_model$breakpoints

bp_years <- df_model$year[bp_indices]

print(bp_years)


plot(df_model$year,
     df_model$log_ind_gdp,
     type="l",
     lwd=3,
     col="darkblue")

for(b in bp_years){
  abline(v=b, col="red", lty=2)
}



# ==============================================================================
# PART 12: ECONOMIC VISUALIZATIONS
# ==============================================================================

cat("\n================================================")
cat("\n PART 12: ECONOMIC VISUALIZATIONS")
cat("\n================================================\n")


# --------------------------------------------------
# Industrial GDP vs Energy Intensity
# --------------------------------------------------

ggplot(df_model,
       aes(x=energy_intensity,
           y=ind_gdp_abs)) +
  
  geom_point(size=3, alpha=0.7, color="#1f78b4") +
  
  geom_smooth(method="lm", se=TRUE, color="red") +
  
  labs(
    title="Industrial Output and Energy Intensity in India",
    subtitle="Higher energy intensity is associated with lower industrial output",
    x="Energy Intensity (Energy Imports per Unit of GDP)",
    y="Industrial GDP (Absolute Value)"
  ) +
  
  theme_minimal(base_size=14)



# --------------------------------------------------
# Industrial GDP vs Coal Share
# --------------------------------------------------

ggplot(df_model,
       aes(x=coal_share,
           y=ind_gdp_abs)) +
  
  geom_point(size=3, alpha=0.7, color="#2ca02c") +
  
  geom_smooth(method="lm", se=TRUE, color="black") +
  
  labs(
    title="Coal Electricity Share and Industrial Output",
    subtitle="Coal dependence shows only a weak relationship with industrial GDP",
    x="Coal Share in Electricity Generation (%)",
    y="Industrial GDP (Absolute Value)"
  ) +
  
  theme_minimal(base_size=14)



# --------------------------------------------------
# Production Function (Capital vs Output)
# --------------------------------------------------

ggplot(df_model,
       aes(x=log_gfcf,
           y=log_ind_gdp)) +
  
  geom_point(aes(size=log_lab), alpha=0.7, color="#08519c") +
  
  geom_smooth(method="lm", color="red") +
  
  labs(
    title="Capital Investment and Industrial Output",
    subtitle="Capital formation is the strongest driver of industrial GDP",
    x="Log Gross Fixed Capital Formation",
    y="Log Industrial GDP",
    size="Labour Participation"
  ) +
  
  theme_minimal(base_size=14)



# --------------------------------------------------
# Long Run Elasticities
# --------------------------------------------------

lr_df <- data.frame(
  variable = rownames(lr_multipliers),
  elasticity = lr_multipliers$Estimate
)

ggplot(lr_df,
       aes(x=variable,
           y=elasticity,
           fill=variable)) +
  
  geom_bar(stat="identity") +
  
  geom_hline(yintercept=0, linetype="dashed") +
  
  labs(
    title="Long-Run Elasticities from the ARDL Model",
    subtitle="Capital investment shows the strongest positive relationship with industrial output",
    x="Model Variable",
    y="Estimated Long-Run Elasticity"
  ) +
  
  theme_minimal(base_size=14) +
  
  theme(legend.position="none")

# ==============================================================================
# PART 13: BLOG VISUALS
# ==============================================================================
grid.arrange(
  g1, g2,
  g3, g4,
  ncol = 2,
  top = textGrob(
    "India's Industrial Growth and Energy Decoupling (1990–2023)",
    gp = gpar(fontsize = 20, fontface = "bold")
  )
)

# --------------------------------------------------
# Decoupling Chart
# --------------------------------------------------

df_plot <- df_model %>%
  mutate(
    output_index = log_ind_gdp / log_ind_gdp[1] * 100,
    intensity_index = log_intensity / log_intensity[1] * 100
  )


ggplot(df_plot, aes(x=year)) +
  
  geom_line(aes(y=output_index,
                color="Industrial Output"),
            linewidth=1.6) +
  
  geom_line(aes(y=intensity_index,
                color="Energy Intensity"),
            linewidth=1.6,
            linetype="dashed") +
  
  labs(
    title="India Produces More Industrial Output With Less Energy",
    subtitle="Industrial output rises while energy intensity steadily declines",
    x="Year",
    y="Index (1990 = 100)",
    color="Indicator"
  ) +
  
  theme_minimal(base_size=14)



# --------------------------------------------------
# Investment vs Industrial Output
# --------------------------------------------------

ggplot(df_model, aes(x=year)) +
  
  geom_line(aes(y=scale(log_ind_gdp),
                color="Industrial Output"),
            linewidth=1.6) +
  
  geom_line(aes(y=scale(log_gfcf),
                color="Capital Investment"),
            linewidth=1.6,
            linetype="dashed") +
  
  labs(
    title="Capital Investment Has Become the Engine of Industrial Growth",
    subtitle="Industrial output closely tracks capital formation",
    x="Year",
    y="Standardized Index",
    color="Variable"
  ) +
  
  theme_minimal(base_size=14)



# --------------------------------------------------
# Coal vs Industrial Growth
# --------------------------------------------------

ggplot(df_model, aes(x=year)) +
  
  geom_line(aes(y=scale(log_ind_gdp),
                color="Industrial Output"),
            linewidth=1.6) +
  
  geom_line(aes(y=scale(log_coal),
                color="Coal Electricity Share"),
            linewidth=1.6,
            linetype="dashed") +
  
  labs(
    title="Industrial Growth Is No Longer Tightly Linked to Coal Use",
    subtitle="Coal electricity share fluctuates while industrial output continues rising",
    x="Year",
    y="Standardized Index",
    color="Variable"
  ) +
  
  theme_minimal(base_size=14)



# --------------------------------------------------
# Structural Break Visualization
# --------------------------------------------------

ggplot(df_model,
       aes(x=year,
           y=log_ind_gdp)) +
  
  geom_line(color="#1f78b4",
            linewidth=1.6) +
  
  geom_vline(xintercept=2003,
             color="red",
             linetype="dashed") +
  
  labs(
    title="Structural Break in India's Industrial Growth",
    subtitle="Electricity Act reforms in 2003 altered the industrial growth trajectory",
    x="Year",
    y="Log Industrial GDP"
  ) +
  
  theme_minimal(base_size=14)

install.packages("gridExtra")
# ============================================================
# LIBRARIES
# ============================================================
# ============================================================
# LIBRARIES
# ============================================================

library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(stringr)

theme_set(theme_minimal(base_size = 13))


# ============================================================
# COLORS
# ============================================================

blue <- "#1f3a5f"
grey <- "#7a7a7a"


# ============================================================
# DATA PREP
# ============================================================

df_plot <- df_model %>%
  mutate(
    output_index = log_ind_gdp / log_ind_gdp[1] * 100,
    intensity_index = log_intensity / log_intensity[1] * 100
  )


# ============================================================
# GRAPH 1
# ============================================================

g1 <- ggplot(df_plot, aes(year)) +
  
  geom_line(aes(y = output_index, color="Industrial Output"), linewidth=1.8) +
  geom_line(aes(y = intensity_index, color="Energy Intensity"), linewidth=1.8) +
  
  scale_color_manual(values=c(blue, grey)) +
  
  labs(
    title="Industrial Output vs Energy Intensity",
    subtitle="India produces more output while using less energy per unit",
    x="Year",
    y="Index (1990 = 100)",
    color=""
  ) +
  
  theme(
    legend.position="bottom",
    plot.title=element_text(face="bold"),
    plot.margin=margin(20,20,20,20)
  )


# ============================================================
# GRAPH 2
# ============================================================

g2 <- ggplot(df_model, aes(year)) +
  
  geom_line(aes(y=scale(log_ind_gdp), color="Industrial Output"), linewidth=1.8) +
  geom_line(aes(y=scale(log_gfcf), color="Capital Investment"), linewidth=1.8) +
  
  scale_color_manual(values=c(blue, grey)) +
  
  labs(
    title="Investment and Industrial Growth",
    subtitle="Capital formation tracks industrial output",
    x="Year",
    y="Standardized Index",
    color=""
  ) +
  
  theme(
    legend.position="bottom",
    plot.title=element_text(face="bold"),
    plot.margin=margin(20,20,20,20)
  )


# ============================================================
# GRAPH 3
# ============================================================

g3 <- ggplot(df_model, aes(year)) +
  
  geom_line(aes(y=scale(log_ind_gdp), color="Industrial Output"), linewidth=1.8) +
  geom_line(aes(y=scale(log_coal), color="Coal Electricity Share"), linewidth=1.8) +
  
  scale_color_manual(values=c(blue, grey)) +
  
  labs(
    title="Coal Use vs Industrial Output",
    subtitle="Industrial growth increasingly diverges from coal dependence",
    x="Year",
    y="Standardized Index",
    color=""
  ) +
  
  theme(
    legend.position="bottom",
    plot.title=element_text(face="bold"),
    plot.margin=margin(20,20,20,20)
  )


# ============================================================
# GRAPH 4
# ============================================================

g4 <- ggplot(df_model, aes(year, log_ind_gdp)) +
  
  geom_line(color=blue, linewidth=1.8) +
  
  geom_vline(
    xintercept=break_year,
    linetype="dashed",
    color=grey,
    linewidth=1.5
  ) +
  
  labs(
    title="Structural Break in Industrial Growth",
    subtitle=paste("Electricity reforms around", break_year),
    x="Year",
    y="Log Industrial GDP"
  ) +
  
  theme(
    plot.title=element_text(face="bold"),
    plot.margin=margin(20,20,20,20)
  )


# ============================================================
# EXPORT DASHBOARD (VERY LARGE CANVAS)
# ============================================================

file_path <- "/Users/anshuman/Desktop/industrial_dashboard.png"

png(file_path,
    width = 4200,
    height = 3200,
    res = 350)

grid.arrange(
  
  arrangeGrob(g1, g2, g3, g4, ncol=2),
  
  top = textGrob(
    "India's Industrial Growth and Energy Decoupling (1990–2023)",
    gp = gpar(fontsize = 20, fontface="bold")
  ),
  
  padding = unit(2, "line")
)

dev.off()

cat("Saved to:", file_path)
