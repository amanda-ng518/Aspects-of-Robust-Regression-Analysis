# Aspects of Robust Regression

This repository contains R code and sample data for simulations and analyses related to the **Robust Regression with Student's $T$: The Role of Degrees of Freedom** paper. 

Linear regression estimators are known to be sensitive to outliers, and one alternative to obtain a robust and efficient estimator of the regression parameter is to model the error with Student's $t$ distribution.  In this article, we compare estimators of the degrees of freedom parameter in the $t$ distribution using frequentist and Bayesian methods, and then study properties of the corresponding estimated regression coefficient. We also include the comparison with some recommended approaches in the literature, including fixing the degrees of freedom and robust regression using the Huber loss.  Our extensive simulations on both synthetic and real data demonstrate that estimating the degrees of freedom via the adjusted profile log-likelihood approach yields regression coefficient estimators with high accuracy, performing comparably to the maximum likelihood estimator where the degrees of freedom are fixed at their true values. These findings provide a detailed synthesis of $t$-based robust regression and underscore a key insight: the proper calibration of the degree of freedom is as crucial as the choice of the robust distribution itself for achieving optimal performance.

---

## Dependencies & Requirements 

The code is written in **R**. To run the scripts successfully, make sure you have the following R packages: 

- **Required R packages** to run locally:
```r
  install.packages(c("MASS", "ggplot2"))
```
---

## Repository Structure

### **R Codes**
This folder contains all R scripts used to generate, run, and analyze the simulations.

| File | Description |
|------|--------------|
| `data_simulation_functions.R` | Functions to simulate data with t errors, normal errors, contaminated errors. |
| `helper_functions.R` | Helper functions that should be preloaded when running $\nu$ or $\beta$ estimations. |
| `nu_estimation_functions.R` | Functions to run $\nu$ estimations. |
| `nu_sim_results.R` | Analyzes and summarizes $\nu$ simulation results about $\nu$. |
| `beta_estimation_function.R` | Functions to run $\beta$ estimations. |
| `beta_sim_results.R` | Summarizes and visualizes $\beta$ estimation results in terms of RMSE, Bias and SE. |

---

### **Sample Data**
This folder includes `.Rda` files generated from simulation runs. Each subfolder corresponds to a specific simulation setting.

| Folder | Description |
|---------|--------------|
| `df data` | Results from $\nu$ estimation simulations. |
| `beta terror data` | Results from $\beta$ estimation under *t*-distributed errors. |
| `beta normal error data` | Results from $\beta$ estimation under normal errors. |
| `beta stackloss data` | Results from $\beta$ estimation using the Stackloss dataset. |
| `beta contaminated error data` | Results from $\beta$ estimation under contaminated errors. |

Each subfolder (e.g., `contaminated_2pterror_1`, `contaminated_2pterror_2`, etc.) contains example `.Rda` files (`allsim_1.Rda`, `allsim_2.Rda`, …, `allsim_10.Rda`) representing repeated simulation runs for robustness evaluation.

---
## Functions

### $\nu$-estimation
The following are the main R functions for estimating degrees of freedom ($\nu$) using four approaches with $\omega = \frac{1}{\nu}$ reparamatrization

- `estimate_nu_profile()`: Profile likelihood  
- `estimate_nu_adj_profile()`: Adjusted profile likelihood
- `estimate_nu_IJ()`: Full Bayes
- `estimate_nu_nu_block()`: Pseudo Bayes

These $\nu$-estimation functions require inputs: `x`, `y` and starting point `omega_init`. They output a list of: 

- `omega`: optimal $\omega$
- `nu`: optimal $\nu$ (essentially 1/ optimal $\omega$)
- `convergence`: code to check BFGS convergence (0 if success)

### $\beta$-estimation
For regression coefficient ($\beta$) estimation, use `estimate_beta()`. Aside from `x`, `y` and starting point `omega_init`, you will need to specify
the method to conduct $\nu$ estimation. The available methods include:

- "OLS": $\nu$ will not be estimated
- "Huber": $\nu$ will not be estimated
- "Profile"
- "Adj profile"
- "Full Bayes"
- "Pseudo Bayes"
- any positive integer (interpreted as a fixed nu): $\hat \nu$ will be set as this fixed integer

It outputs a list of : 

- `method`: method used to estimate `nu`, same as your input `method`
- `nu_hat`: optimal $\nu$ (essentially 1/ optimal $\omega$)
- `beta_hat`: optimal $\beta$
- `sigma_hat`: optimal $\sigma$
- `success_nu`: Convergence code for $\nu$-estimation (0 if success)
- `success_beta`: Convergence code for $\beta$-estimation (0 if success)
    
---

## Example

The code below simulates data with t error. Covariates x are generated from standard normal distribution.
```r
# Simulate a set of data with 20 rows, 5 columns with true beta all set as 0 and t(2) error
sim_data <- simulate_t_error_data(n = 20, p = 5, beta = rep(0,5), sigma = 1, nu = 2, seed = 123)
x <- sim_data$x
y <- sim_data$y
```

The code below simulates data with standard normal error (controlled by `mean` and `sigma`) and 20% contaminated two-point errors.
```r
# Simulate a set of data with 20 rows, 5 columns with true beta all set as 0 and t(2) error
contam_sim_data <- simulate_contaminated_data(n = 20, p = 5, beta = rep(0, 5), mean = 0,
                                       sigma = 1, contam_type = "twopt", NULL, contam_prob = 0.2,
                                       seed = NULL)
x <- contam_sim_data$x
y <- contam_sim_data$y
```

Note: Please load all functions in `helper_functions.R` before proceeding to the estimation steps.

The code below estimates $\nu$ using the profile likelihood approach.

```r
# Estimate nu with profile likelihood approach
nu_estimation <- estimate_nu_profile(y, x, omega_init = 1/2)
est_nu <- nu_estimation$nu # Optimal nu
nu_estimation$convergence # Check convergence = 0 (if success)
```

The code below estimates $\beta$ using the adjusted profile likelihood approach. 

```r
# Estimate beta with adjusted profile likelihood approach
beta_estimation <- estimate_beta(y, x, method = "Adj profile", omega_init = 0.5) # By default we use OLS estimates as initial guess for beta and sigma
est_beta <- beta_estimation$beta_hat # Optimal beta
nu_hat <- beta_estimation$nu_hat # Optimal nu
beta_estimation$success_beta # Check convergence = 0 (if success)
```












----
