# ============================================================
# Beta estimation function
# ============================================================
# Pre-requisites: Functions in helper_functions.R and nu_estimation_functions.R must be pre-loaded.
# Pre-requisites: MASS, rlm library

# FUNCTION: Estimate Regression Coefficients beta Using Various Methods
# INPUT: y Numeric vector. The response variable.
# INPUT: x Numeric matrix. The design matrix (must include an intercept column).
# INPUT: method Character or numeric. Available methods:
#        "OLS", "Huber", "Profile", "Adj profile", "Full Bayes", "Pseudo Bayes",
#        or a positive integer (interpreted as a fixed nu).
# INPUT: omega_init Numeric. Initial value for the omega = 1/nu parameter (default is 0.5).
# INPUT: beta_init Optional numeric vector. Initial values for beta. If not specified, ordinary least squares (OLS) estimates are used as default.
# INPUT: sigma_init Optional numeric value. Initial value for sigma. If not specified, the residual standard deviation from OLS is used as the default.
# INPUT: control Optional list of control parameters passed to optim.
#
# OUTPUT: method Character. The estimation method used.
# OUTPUT: nu_hat Estimated or fixed degrees of freedom nu. NA for OLS and Huber.
# OUTPUT: omega_hat Estimated or fixed value of omega = 1/nu. NA for OLS and Huber.
# OUTPUT: beta_hat Estimated regression coefficients. 
# OUTPUT: sigma_hat Estimated scale parameter (residual standard deviation).
# OUTPUT: success_nu Convergence code for nu estimation. 0 indicates successful convergence; NA for OLS, Huber, and fixed nu.
# OUTPUT: success_beta Convergence code for beta-sigma optimization. 0 indicates success; NA for OLS and Huber.

estimate_beta <- function(y, x,
                          method = c("OLS", "Huber", "Profile", "Adj profile", "Full Bayes",
                                     "Pseudo Bayes"),
                          omega_init = 0.5,
                          beta_init = NULL,
                          sigma_init = NULL,
                          control = list(maxit = 10000)) {
  method <- as.character(method)
  n <- length(y)
  p <- ncol(x)
  
  # ---- CASE 1: OLS ----
  if (method == "OLS") {
    lmfit <- lm(y ~ x - 1)
    beta_hat <- coef(lmfit)
    sigma_hat <- sd(residuals(lmfit))
    return(list(
      method = method,
      nu_hat = NA,
      omega_hat = NA,
      beta_hat = beta_hat,
      sigma_hat = sigma_hat,
      success_nu = NA,
      success_beta = NA
    ))
  }
  
  # ---- CASE 2: Huber regression ----
  if (method == "Huber") {
    fit_lm <- lm(y ~ x - 1)
    RB <- rlm(y ~ x[, -1], psi = psi.huber, k = 1.345, maxit = 100)
    beta_huber <- rlmDD(y, x[, -1], fit_lm$coef[-2], RB$coef,
                        method = "Huber", plot = "N")$esti[["coefficients"]]
    sigma_hat <- sd(y - as.vector(x %*% beta_huber))
    return(list(
      method = method,
      nu_hat = NA,
      omega_hat = NA,
      beta_hat = beta_huber,
      sigma_hat = sigma_hat,
      success_nu = NA,
      success_beta = NA
    ))
  }
  
  # ---- CASE 3: Fixed nu or nu-estimation method ----
  
  if (is.numeric(method)) method <- as.character(method)
  valid_methods <- c("Profile", "Adj profile", "Full Bayes", "Pseudo Bayes")
  
  # 3a. Estimate nu if method is one of the 5 estimation methods
  if (method %in% valid_methods) {
    est <- switch(method,
                  "Profile" = estimate_nu_profile(y, x, omega_init),
                  "Adj profile" = estimate_nu_adj_profile(y, x, omega_init),
                  "Full Bayes" = estimate_nu_IJ(y, x, omega_init),
                  "Pseudo Bayes" = estimate_nu_nu_block(y, x, omega_init)
    )
    omega_hat <- est$omega
    nu_hat <- ifelse(omega_hat > 0, 1 / omega_hat, Inf)
    success_nu <- est$convergence
  } else {
    # 3b. Fixed nu
    if (suppressWarnings(!is.na(as.numeric(method))) &&
        as.numeric(method) > 0) {
      nu_hat <- as.numeric(method)
      omega_hat <- ifelse(nu_hat == 0, Inf, 1 / nu_hat)
    } else {
      stop("Invalid method: must be 'OLS', 'Huber', one of the 4 nu estimation methods, or a non-negative integer.")
    }
    success_nu <- NA
  }
  
  # 3c. Optimize beta and sigma given nu
  fit_res <- fit_profile_mle_fixed_nu(y, x, nu_hat, par_init = NULL, control = control)
  
  list(
    method = method,
    nu_hat = nu_hat,
    omega_hat = omega_hat,
    beta_hat = fit_res$beta,
    sigma_hat = fit_res$sigma,
    success_nu = success_nu,
    success_beta = fit_res$convergence
  )
}
