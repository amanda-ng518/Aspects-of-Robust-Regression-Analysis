# Aspects of Robust Regression

This repository contains code and sample data for simulations and analyses related to the **Aspects of Robust Regression** project. The goal of this project is to evaluate the robustness of **degrees-of-freedom ($\nu$)** and **regression coefficient ($\beta$)** estimators under different error distributions and contamination settings.

---

## Project Overview

This project compares **Bayesian** and **frequentist** approaches to robust regression through a combination of theoretical analysis and simulation studies, with an emphasis on **high-dimensional regression with regularization**.

The comparison is motivated by both:
- foundational differences in inference between Bayesian and frequentist frameworks, and  
- the practical need to assess the reliability and robustness of Bayesian methods.  

Robust regression methods help ensure that statistical conclusions remain valid even when the model used for inference deviates from the data-generating process. High-dimensional regression poses additional challenges, making **robustness** a particularly important focus in this work.

In regression analysis, the primary objective is typically the estimation of the regression coefficients $\beta$, which characterize the relationship between the predictors and the response variable. In this project, we propose a two-stage procedure under the student's-*t* regression framework as follows: we first estimate the degrees of freedom parameter $\nu$, and then use this estimate to infer the regression coefficients $\beta$.

We explored five Frequentist and Bayesian approaches for estimating $\nu$ under the Student's-*t* regression model and conduct simulations to evaluate their performance. From the simulations, we demonstrate that our proposed procedure yields superior estimation of $\beta$ compared to existing methods, using both simulated and real datasets. Notably, in high-dimensional settings, the adjusted profile likelihood approach performs comparably to the widely-used robust regression method, Huber regression. We also assessed the robustness of our proposed methods in the presence of outliers.

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
| `beta_contaminated_error_sim.R` | Simulates $\beta$ estimation under contaminated non-*t* error structures. |
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
| `beta contaminated error data/` | Results from $\beta$ estimation under contaminated non-*t* errors. |

Each subfolder (e.g., `contaminated_2pterror_1`, `contaminated_2pterror_2`, etc.) contains example `.Rda` files (`allsim_1.Rda`, `allsim_2.Rda`, …, `allsim_10.Rda`) representing repeated simulation runs for robustness evaluation.

---
