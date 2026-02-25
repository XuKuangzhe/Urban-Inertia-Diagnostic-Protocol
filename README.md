# Urban Structural Inertia Diagnostic Protocol (USIDP)

![Status](https://img.shields.io/badge/Status-Active-success)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/Made%20with-R-blue.svg)](https://www.r-project.org/)
[![Stan](https://img.shields.io/badge/Powered%20by-Stan-red.svg)](https://mc-stan.org/)

> **A Bayesian framework for diagnosing structural inertia and intervention efficacy in data-scarce urban energy systems.**

This repository hosts the source code and diagnostic protocols for the research paper: **"A Bayesian Framework to Diagnose Structural Inertia in Data-Scarce Megacities: Proof-of-Concept in Beijing"**.

---

## 📖 Overview

Rapid urbanization in the Global South often "locks in" high-carbon infrastructure. Diagnosing this **Structural Inertia** is challenging due to severe data scarcity ($N < 20$) and noise from exogenous shocks (e.g., COVID-19).

This protocol adapts **anomaly detection algorithms** from signal processing (State-Space Models) to the domain of urban energy policy. By utilizing a **Bayesian Structural Time Series (BSTS)** framework, we isolate latent structural trends from transient interventions, providing a robust diagnostic tool for cities with limited historical data.

### Key Capabilities
* **Small-Sample Robustness:** Validated parameter recovery even with $N \approx 10-15$ data points via Hierarchical Bayesian LASSO.
* **Shock Isolation:** Explicitly models intervention effects (e.g., pandemic lockdowns) to prevent estimation bias.
* **Regime Classification:** Distinguishes between "Structural Saturation," "Stochastic Drift," and "Elastic Growth."
* **Policy Stress Testing:** Simulates "Regime Switch" scenarios (e.g., -40% structural shock) to identify effective intervention thresholds.

---

## 📂 Repository Structure

The codebase is organized to separate statistical models from execution scripts:

```bash
.
├── R_code/
│   ├── Main_Analysis.R        # Primary pipeline: Data loading, Model fitting (Stan), and LOO-CV
│   ├── SI_Analysis.R          # Supplementary analysis (Sensitivity tests, alternative priors)
│   └── HDI.R                  # Helper function for Highest Density Intervals
├── Stan Model/
│   ├── SSM_Full_Dynamic.stan  # Full Model: STIRPAT (Pop, GDP, Tech) + Dynamic Trend
│   └── SSM_Null_Dynamic.stan  # Null Model: Second-order Random Walk (Pure Inertia)
├── data/
|   ├── tableA4.csv
│   └── tableA5.csv 
└── README.md

```

## 🚀 Getting Started
Prerequisites
To run this protocol, you will need R (version >= 4.0) and the following packages:

```R
install.packages(c("rstan", "loo", "tidyverse", "bayestestR", "ggplot2", "patchwork"))
```

Note: This project relies on rstan for MCMC sampling. Ensure your C++ toolchain (Rtools on Windows, Xcode on Mac) is correctly configured.

Usage Guide (Generalizability)
This protocol is city-agnostic. Researchers can apply it to other data-scarce cities (e.g., Lagos, Delhi, São Paulo) by following these steps:

1. Clone the repository

```bash
git clone [https://github.com/YourUsername/Urban-Inertia-Diagnostic-Protocol.git](https://github.com/YourUsername/Urban-Inertia-Diagnostic-Protocol.git)
```

2. Prepare your data
Format your city's time-series data as a .csv file with the following columns:

Year: Integer (e.g., 2010-2023)

Consumption: Log-transformed energy consumption

Covariates: Log-transformed drivers (Population, GDP, Secondary Industry Share, etc.) - Required for Full Model

Shock: Binary (0/1) or decay vector for intervention periods (e.g., COVID-19)

3. Run the Analysis
Open R_code/Main_Analysis.R and point the data loader to your local file:

```R
# In Main_Analysis.R
# Replace with your city's dataset path
city_data <- read.csv("path/to/your_city_data.csv")

# Run the Bayesian sampler
fit <- stan(file = "Stan Model/SSM_Null_Dynamic.stan", data = city_data, ...)
```

The script will output posterior distributions for structural inertia ($\sigma_{\mu}$) and the impact magnitude of exogenous shocks ($\beta_{shock}$).

## 📊 Data Availability
While this repository hosts the diagnostic algorithms, the full replication datasets (including detailed energy balance sheets for Beijing 2010-2023) and raw outputs are archived for long-term preservation at the Open Science Framework (OSF):

👉 Access Full Dataset on https://osf.io/v7j3g/overview

## 📝 License
This project is licensed under the MIT License - see the LICENSE file for details. This ensures the protocol remains open and adaptable for the global research community.

Contact: Kuangzhe Xu Cyberspace Security University of China 


