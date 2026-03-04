# Aspects of Robust Regression

This repository contains R code and sample data for simulations and analyses related to the **Aspects of Robust Regression** project. 

---

## Project Overview

Linear regression estimators are known to be sensitive to outliers, and one alternative to obtain a robust and efficient estimator of the regression parameter is to model the error with Student's $t$ distribution.  In this article, we compare estimators of the degrees of freedom parameter in the $t$ distribution using frequentist and Bayesian methods, and then study properties of the corresponding estimated regression coefficient. We also include the comparison with some recommended approaches in the literature, including fixing the degrees of freedom and robust regression using the Huber loss. 
Our extensive simulations on both synthetic and real data demonstrate that estimating the degrees of freedom via the adjusted profile log-likelihood approach yields regression coefficient estimators with high accuracy, performing comparably to the maximum likelihood estimator where the degrees of freedom are fixed at their true values. These findings provide a detailed synthesis of $t$-based robust regression and underscore a key insight: the proper calibration of the degree of freedom is as crucial as the choice of the robust distribution itself for achieving optimal performance.

---

## Methods & Models

The analyses in this project include:

- **Estimation of degrees of freedom ($\nu$)** using five approaches with $\omega = \frac{1}{\nu}$ reparamatrization:  
  - Profile likelihood  
  - Adjusted profile likelihood  
  - Independence Jeffreys  
  - Marginal Jeffreys  
  - Marginal Fisher  

  These estimated $\nu$ values were then used in *t*-regression models to estimate $\beta$ coefficients.

- **Comparison of regression methods:**  
  - Traditional methods: OLS, *t*-regression (fixed $\nu$), and Huber regression  
  - Our proposed Two-stage methods using estimated $\nu$

- **Evaluation metrics:**  
  - Bias  
  - Variance  
  - Root mean squared error (RMSE)

Simulation settings included:
- Heavy-tailed *t*(2) and normal errors  
- Real-data Stackloss experiments  
- High-dimensional settings (varying $\frac{p}{n}$ ratios)  
- Errors from $N(0,9)$, *t*(2), $\chi^2(4)-4$, and two-point (-5 or 5) with contamination rates of 10%, 20%, and 30%  

--- 

## Dependencies & Requirements 

The code is written in **R**. To run the scripts successfully, make sure you have the following R packages: 

- **Required R packages** to run locally:
```r
  install.packages(c(
    "MASS",
    "ggplot2"
  ))
```

- **Additional R packages** to run on server:
```r
  install.packages(c(
    "doParallel",
    "foreach",
    "rlmDataDriven" ))
```
---

## Repository Structure

### **R Codes/**
This folder contains all R scripts used to generate, run, and analyze the simulations.

| File | Description |
|------|--------------|
| `df_sim.R` | Runs degrees-of-freedom ($\nu$) estimation simulations. |
| `df_sim_results.R` | Analyzes and summarizes df simulation results. |
| `beta_tnormal_error_sim.R` | Simulates regression coefficient ($\beta$) estimation with *t*-distributed and normal errors. |
| `beta_stackloss_sim.R` | Simulates $\beta$ estimation using the real Stackloss dataset. |
| `beta_contaminated_error_sim.R` | Simulates $\beta$ estimation under contaminated error structures. |
| `beta_sim_results.R` | Summarizes and visualizes $\beta$ estimation results. |

---

### **Sample Data/**
This folder includes `.Rda` files generated from simulation runs. Each subfolder corresponds to a specific simulation setting.

| Folder | Description |
|---------|--------------|
| `df data/` | Results from df estimation simulations. |
| `beta terror data/` | Results from $\beta$ estimation under *t*-distributed errors. |
| `beta normal error data/` | Results from $\beta$ estimation under normal errors. |
| `beta stackloss data/` | Results from $\beta$ estimation using the Stackloss dataset. |
| `beta contaminated error data/` | Results from $\beta$ estimation under contaminated errors. |

Each subfolder (e.g., `contaminated_2pterror_1`, `contaminated_2pterror_2`, etc.) contains example `.Rda` files (`allsim_1.Rda`, `allsim_2.Rda`, …, `allsim_10.Rda`) representing repeated simulation runs for robustness evaluation.

---
