# Aspects of Robust Regression

This repository contains code and sample data for simulations and analyses related to the **Aspects of Robust Regression** project. The goal of this project is to explore the robustness of degree of freedom ($\nu$) and regression coefficient ($\beta$) estimators under different error distributions and contamination settings.

---

## 📂 Repository Structure

### **R Codes/**
This folder contains all R scripts used to generate, run, and analyze the simulations.

| File | Description |
|------|--------------|
| `df_sim.R` | Conducts degrees-of-freedom (df) estimation simulations.|
| `df_sim_results.R` | Analyzes and summarizes df simulation outcomes. |
| `beta_tnormal_error_sim.R` | Conducts regression coefficient ($\beta$) estimation simulations on data with t-distributed and normal errors. |
| `beta_stackloss_sim.R` | Conducts regression coefficient ($\beta$) estimation simulations on real Stackloss dataset. |
| `beta_contaminated_error_sim.R` | Conducts regression coefficient ($\beta$) estimation simulations on data with contaminated non-t errors. |
| `beta_sim_results.R` | Summarizes and visualizes simulation results for regression coefficient estimates. |

---

### **Sample Data/**
This folder includes data files in `.Rda` produced from simulation runs. Each subfolder corresponds to a specific simulation setting.

| Folder | Description |
|---------|--------------|
| `df data/` | Contains simulation results about df estimation. |
| `beta terror data/` | Contains simulation results about $\beta$ estimation in data with t error. |
| `beta normal error data/` | Contains simulation results about $\beta$ estimation in data with normal error. |
| `beta stackloss data/` | Contains simulation results about $\beta$ estimation in Stackloss dataset. |
| `beta contaminated error data/` | Contains beta estimation simulation results in data with contaminated non-t error structures. |

Each subfolder (e.g., `contaminated_2pterror_1`, `contaminated_2pterror_2`, etc.) stores a small portion of sample `.Rda` files (`allsim_1.Rda`, `allsim_2.Rda`, …, `allsim_10.Rda`) corresponding to repeated runs for robustness evaluation. 

---

## 🧠 Project Overview

This project considers a **comparison of Bayesian and frequentist methods** for robust regression. It involves a mix of theoretical analysis and simulations, with an emphasis on **high-dimensional regression with regularization**.  Comparing Bayesian and frequentist approaches is of interest both for the **foundations of inference** and for the **practical assessment of the reliability of Bayesian approaches**; the latter is closely related to the asymptotic theory of likelihood-based inference. Robust regression methods are an important technique for ensuring that statistical conclusions remain valid even when the model used for inference differs from the model generating the data.  New methods and theory for high-dimensional regression are an area of active development, and **robustness in this setting** is of particular interest.

Robust regression methods aim to reduce the influence of outliers and deviations from model assumptions (e.g., normality of errors).  
In this project, we:

- Compare classical and robust regression estimators under various error distributions.  
- Evaluate sensitivity to contamination and heavy tails using simulation studies.  
- Summarize estimator performance through bias, variance, and mean squared error (MSE) metrics.  

---

## 🧩 Methods & Models

The project explores:
- **Ordinary Least Squares (OLS)** under ideal conditions.  
- **Robust alternatives** (e.g., M-estimators) under contamination or t-distributed noise.  
- The impact of sample contamination on estimator efficiency.  

Simulation setups include:
- Contaminated normal errors  
- Heavy-tailed t-distributed errors  
- Real-data-based Stackloss experiments  

---

## ⚙️ How to Run

1. Open the `R Codes/` directory.
2. Run the desired simulation file (e.g., `beta_tnormal_error_sim.R`) in R.
3. Results are automatically saved in the corresponding subfolder under `Sample Data/`.
4. Use `beta_sim_results.R` or `df_sim_results.R` to summarize or visualize simulation outcomes.

Example:
```r
# Run simulation with t-distributed errors
source("R Codes/beta_tnormal_error_sim.R")

# Summarize results
source("R Codes/beta_sim_results.R")

---

## ⚙️ Dependencies & Requirements

The code is written in **R**.  
To run the scripts successfully, make sure you have the following environment:

- **R version** ≥ 4.2.0  
- **Required R packages**:
  ```r
  install.packages(c(
    "MASS",
    "ggplot2",
    "dplyr",
    "tidyr",
    "gridExtra",
    "reshape2"
  ))

---

