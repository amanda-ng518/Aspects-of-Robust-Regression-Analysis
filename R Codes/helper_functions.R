# ============================================================
# Helper functions 
# ============================================================

# FUNCTION: Profile MLE of sigma and beta for Fixed Degrees of Freedom
# This helper function computes the profile maximum likelihood estimates (MLEs)
# of the regression coefficients (beta) and the scale parameter (sigma)
# for a fixed degrees of freedom nu under a Student-t error model.
# It returns the estimated parameters along with the corresponding negative log-likelihood
# value evaluated at those estimates.
#
# INPUT: nu Numeric. The fixed degrees of freedom for the t-distribution.
# INPUT: y Numeric vector. The response variable.
# INPUT: x Numeric matrix. The design matrix (must include an intercept column).
# INPUT: par_init Numeric vector. Initial value for beta and sigma optimization. Default is OLS estimates.
# OUTPUT: beta_hat Profile MLE of beta
# OUTPUT: sigma_hat Profile MLE of sigma
# OUTPUT: loglik Negative log-likelihood value evaluated at the fitted parameters.
# OUTPUT: convergence Convergence code returned by the optimizer (0 indicates successful convergence).

fit_profile_mle_fixed_nu <- function(y, x, nu = Inf, par_init = NULL, control = list()) {
  n <- length(y)
  p <- ncol(x)
  if (is.null(par_init)) {
    lmfit <- lm(y ~ x - 1) # x includes intercept
    sigma_init <- stats::sd(residuals(lmfit))
    par_init <- c(coef(lmfit), sigma_init)
  }
  treg_neg_loglik_no_nu <- function(par) {
    beta <- par[1:p]
    sigma <- abs(par[p+1])
    z <- (y - x %*% beta)/sigma
    if (is.infinite(nu)) {
      ll <- sum(dnorm(z, log = TRUE)) - n * log(sigma)
    } else {
      ll <- sum(dt(x = z, df = nu, log = TRUE)) - n * log(sigma)
    }
    -ll
  }
  opt <- optim(par = par_init, fn = treg_neg_loglik_no_nu, method = "BFGS", control = c(list(maxit = 10000), control))
  par_hat <- opt$par
  list(beta = par_hat[1:p], sigma = abs(par_hat[p+1]), loglik = -opt$value, convergence = opt$convergence)
}

# FUNCTION: Profile MAP of sigma and beta for Fixed Degrees of Freedom
# This helper function computes the maximum a posteriori (MAP) estimates of the regression
# coefficients (beta) and the scale parameter (sigma) for a fixed
# value of the degrees of freedom parameter nu in a Student-t regression model.
# This function is used as the first-stage optimization step in
# estimating the full joint posterior over all model parameters, including nu,
# where the full objective includes an additional prior on nu.
#
# INPUT: nu Numeric. The fixed degrees of freedom for the t-distribution.
# INPUT: y Numeric vector. The response variable.
# INPUT: x Numeric matrix. The design matrix (must include an intercept column).
# INPUT: par_init Numeric vector. Initial value for beta and sigma optimization. Default is OLS estimates.
# OUTPUT: beta_hat Profile MAP of beta
# OUTPUT: sigma_hat Profile MAP of sigma
# OUTPUT: loglik Value of the log-posterior objective (excluding the nu prior) evaluated at the estimates.
# OUTPUT: convergence Convergence code returned by the optimizer (0 indicates successful convergence).

fit_profile_map_fixed_nu <- function(y, x, nu, par_init = NULL, control = list(maxit = 10000)) {
  n <- length(y)
  p <- ncol(x)
  
  # Default initialization via OLS
  if (is.null(par_init)) {
    lmfit <- lm(y ~ x - 1)  # x includes intercept
    beta_init <- coef(lmfit)
    sigma_init <- sd(residuals(lmfit))
    par_init <- c(beta_init, sigma_init)
  }
  
  # Negative log-likelihood + log sigma prior for fixed nu
  treg_neg_loglik_sigmaprior_no_nu <- function(par) {
    beta <- par[1:p]
    sigma <- abs(par[p + 1])
    
    z <- as.vector((y - x %*% beta) / sigma)
    if (is.infinite(nu)) {
      ll <- sum(log(dnorm(z))) - (n + 1) * log(sigma)
    } else {
      ll <- sum(log(dt(z, df = nu))) - (n + 1) * log(sigma)
    }
    -ll
  }
  opt <- optim(par = par_init, fn = treg_neg_loglik_sigmaprior_no_nu, method = "BFGS", control = c(list(maxit = 10000), control))
  par_hat <- opt$par
  list(beta = par_hat[1:p], sigma = abs(par_hat[p+1]), loglik = -opt$value, convergence = opt$convergence)
}

# FUNCTION: Optimize omega parameters for a given target function
# This helper function performs numerical optimization over the
# omega parameter vector using the BFGS algorithm. It is primarily designed
# for minimizing a omega-based objective function within the function's
# estimation routines.

# INPUT: omega_init Numeric vector of initial values for the omega.
# INPUT: target_fn A function to be minimized, taking omega as input and returning a scalar objective value.
# OUTPUT: omega The optimized omega parameter vector.
# OUTPUT: convergence Convergence code returned by optim. A value of 0 indicates successful convergence.}
# OUTPUT: value The objective function value evaluated at the optimum.

optimize_omega <- function(omega_init, target_fn) {
  out <- optim(par = omega_init, fn = target_fn, method = "BFGS", control = list(maxit = 10000))
  list(omega = out$par, convergence = out$convergence, value = out$value)
}
