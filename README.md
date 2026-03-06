
# The Dual Metabolism of India  
## Macro-Lock-in vs Industrial Decoupling in the Energy–Growth Nexus (1990–2024)

This repository presents a comprehensive **econometric investigation of India's energy–growth nexus** using a **dual-track modeling framework**. The project evaluates the structural relationship between **economic growth and energy dependence** at two distinct levels of the economy:

1. **Aggregate Macroeconomic System**
2. **Formal Industrial Sector**

The study introduces the concept of a **"Dual Metabolism"** in the Indian economy:

- A **biophysically constrained macroeconomic system** exhibiting energy lock-in.
- A **potentially adaptive industrial core** capable of partial energy decoupling.

Two econometric frameworks are used:

| Model | Econometric Framework | Target System |
|------|----------------------|--------------|
| Model 1 | Vector Error Correction Model (VECM) | Aggregate Macroeconomy |
| Model 2 | Autoregressive Distributed Lag (ARDL) | Formal Industrial Sector |

---

# Repository Structure
```text
Dual-Metabolism-India/
│
├── Codes/
│   ├── MODEL1_FINAL.R
│   └── MODEL2_FINAL.R
│
├── Dataset/
│   └── QoG_TimeSeries.csv
│
├── Graphs/
│   ├── Model1/
│   │   ├── IRF/
│   │   ├── FEVD/
│   │   └── Stability/
│   │
│   └── Model2/
│       ├── LongRun/
│       ├── Diagnostics/
│       └── IRF/
│
├── Results/
│   ├── Regression_Tables/
│   ├── Diagnostics/
│   └── Summary_Statistics/
│
└── README.md
```

---

# System Requirements

To reproduce the analysis, ensure your system meets the following requirements.

### Software
- **R ≥ 4.0**
- **RStudio (recommended)**

### Operating Systems
- Windows 10+
- macOS
- Linux

### Hardware

Minimum:

- 8GB RAM
- 5GB available storage

Recommended:

- 16GB RAM for faster IRF simulations

---

# Required R Libraries

Install the following packages before running the scripts.

```r
install.packages(c(
"tidyverse",
"urca",
"tseries",
"vars",
"ARDL",
"lmtest",
"car"
))
```
## Library Functions

| Package | Function |
|--------|---------|
| tidyverse | Data cleaning and manipulation |
| urca | Unit root and Johansen cointegration testing |
| tseries | Statistical time-series tests |
| vars | VAR and VECM estimation |
| ARDL | Bounds testing and ARDL estimation |
| lmtest | Diagnostic testing |
| car | Regression diagnostics |

---

## Dataset

The analysis uses the **Quality of Government (QoG) Standard Time-Series Dataset**.

### Variables Extracted

- Real GDP per capita  
- Gross Fixed Capital Formation  
- Energy intensity  
- Industrial output  
- Labour force participation  
- Coal electricity share  

---

## ⚠️ Disclaimer

This repository **does not claim ownership of the QoG dataset**.

Users should refer to the official repository for updates.

**QoG Institute:**  
https://www.qog.pol.gu.se

---

# Methodological Framework

The research separates the economy into **two structural subsystems**.

| System | Modeling Strategy |
|------|------------------|
| Aggregate economy | System dynamics (VECM) |
| Industrial sector | Production function (ARDL) |

This separation is motivated by the observation that **macro-level thermodynamic constraints differ from sectoral production dynamics**.

---

# Model 1 — Aggregate Biophysical Framework (VECM)

## Theoretical Foundation

Model 1 follows **biophysical macroeconomics**, which argues that economic output is fundamentally constrained by **energy availability and capital formation**.

Macroeconomic aggregates are often **non-stationary but cointegrated**, requiring the use of a **Vector Error Correction Model (VECM)**.

---

# Model 2 — Industrial Decoupling Framework (ARDL)

## Theoretical Foundation

Industrial output is modeled using a **Cobb–Douglas production function augmented with an energy input**.

The model is estimated using the **Pesaran–Shin–Smith (2001) ARDL bounds testing approach**.

### Why ARDL?

ARDL is chosen because it:

- Handles mixed integration orders **I(0)** and **I(1)**
- Performs well in **small samples**
- Allows inclusion of **structural breaks**
