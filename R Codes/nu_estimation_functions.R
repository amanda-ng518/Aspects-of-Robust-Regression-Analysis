# ============================================================
# Nu estimation functions
# ============================================================
# Pre-requisites: Functions in helper_functions.R must be pre-loaded.


# FUNCTION: Profile Log-Likelihood Estimator 
# This function computes the profile log-likelihood estimator for the degrees of freedom
# parameter nu via optimization over omega = 1/nu.
#
# INPUT: y Numeric vector. The response variable.
# INPUT:  x Numeric matrix. The design matrix (should include an intercept column).
# INPUT:  omega_init Numeric. Initial value for omega. If simulated t-error data is used, setting this to the true value of 1/nu can improve convergence speed.
# OUTPUT: omega The optimized value of omega.
# OUTPUT: nu The corresponding value of nu = 1/omega.
# OUTPUT: convergence An integer convergence code from optim. A value of 0 indicates successful convergence.}

estimate_nu_profile <- function(y, x, omega_init = 1/2) {
  n <- length(y)
  p <- ncol(x)
  # objective returns negative profile loglik given omega
  profile_neg <- function(omega) {
    if (omega < 0) return(Inf)
    if (omega == 0) {
      fit <- fit_profile_mle_fixed_nu(y = y, x = x, nu = Inf)
    } else {
      nu <- 1 / omega
      fit <- fit_profile_mle_fixed_nu(y = y, x = x, nu = nu)
    }
    -fit$loglik
  }
  out <- optimize_omega(omega_init, profile_neg)
  omega <- out$omega
  list(omega = omega, nu = ifelse(omega > 0, 1/omega, Inf),
       convergence = out$convergence, value = out$value)
}

# FUNCTION: Adjusted-profile log likelihood estimator
# This function computes the adjusted profile log-likelihood estimator for the degrees of freedom
# parameter nu via optimization over omega = 1/nu.
#
# INPUT: y Numeric vector. The response variable.
# INPUT:  x Numeric matrix. The design matrix (should include an intercept column).
# INPUT:  omega_init Numeric. Initial value for omega. If simulated t-error data is used, setting this to the true value of 1/nu can improve convergence speed.
# OUTPUT: omega The optimized value of omega.
# OUTPUT: nu The corresponding value of nu = 1/omega.
# OUTPUT: convergence An integer convergence code from optim. A value of 0 indicates successful convergence.}

estimate_nu_adj_profile <- function(y, x, omega_init = 1/2) {
  n <- length(y); p <- ncol(x)
  adj_neg <- function(omega) {
    if (omega < 0) return(Inf)
    if (omega == 0) {
      fit <- fit_profile_mle_fixed_nu(y, x, nu = Inf)
      # j_mat for normal case
      j11 <- matrix(0, nrow = p, ncol = p)
      for (ii in 1:n) j11 <- j11 + x[ii,] %*% t(x[ii,])
      j11 <- j11 / (fit$sigma^2)
      j12 <- matrix(0, ncol = 1, nrow = p)
      j22 <- 2 * n / (fit$sigma^2)
      j_mat <- rbind(cbind(j11, t(j12)), cbind(j12, j22))
      adj <- -0.5 * log(abs(det(j_mat)))
      - (fit$loglik + adj)
    } else {
      nu <- 1/omega
      fit <- fit_profile_mle_fixed_nu(y, x, nu = nu)
      # compute j_mat
      j11 <- matrix(0, nrow = p, ncol = p)
      j12 <- rep(0, p)
      j22 <- 0
      resid <- as.numeric(y - x %*% fit$beta)
      for (ii in 1:n) {
        tmp1 <- resid[ii]
        tmp2 <- nu * fit$sigma^2 + tmp1^2
        j11 <- j11 + (x[ii,] %*% t(x[ii,])) * (tmp2 - 2*tmp1^2) / (tmp2^2)
        j12 <- j12 + x[ii,] * 2 * nu * fit$sigma * tmp1 / (tmp2^2)
        j22 <- j22 + 2 * nu * tmp1^2 / ( (nu * fit$sigma^2 + tmp1^2)^2 )
      }
      j11 <- (nu + 1) * j11
      j12 <- (nu + 1) * j12
      j22 <- (nu + 1) * j22
      j_mat <- rbind(cbind(j11, matrix(j12, ncol = 1)), cbind(matrix(j12, nrow = 1), j22))
      adj <- -0.5 * log(abs(det(j_mat)))
      -(fit$loglik + adj)
    }
  }
  out <- optimize_omega(omega_init, adj_neg)
  omega <- out$omega
  list(omega = omega, nu = ifelse(omega > 0, 1/omega, Inf),
       convergence = out$convergence, value = out$value)
}

# FUNCTION: Full Bayes estimator
# This function computes the Full Bayes approach estimator for the degrees of freedom
# parameter nu via optimization over omega = 1/nu.
#
# INPUT: y Numeric vector. The response variable.
# INPUT:  x Numeric matrix. The design matrix (should include an intercept column).
# INPUT:  omega_init Numeric. Initial value for omega. If simulated t-error data is used, setting this to the true value of 1/nu can improve convergence speed.
# OUTPUT: omega The optimized value of omega.
# OUTPUT: nu The corresponding value of nu = 1/omega.
# OUTPUT: convergence An integer convergence code from optim. A value of 0 indicates successful convergence.}

estimate_nu_IJ <- function(y, x, omega_init = 1/2) {
  n <- length(y)
  p <- ncol(x)
  obj_neg <- function(omega) {
    if (omega < 0) return(Inf)
    if (omega == 0) {
      fit <- fit_profile_map_fixed_nu(y, x, nu = Inf)
      joint_logprior <- 0
      -(fit$loglik + joint_logprior)
    } else {
      nu <- 1/omega
      fit <- fit_profile_map_fixed_nu(y, x, nu = nu)
      # Nu Jeffreys prior approximation
      jj <- (nu/(nu+3))*(trigamma(nu/2)-trigamma(nu/2+1/2)-2*(nu+3)/(nu*(nu+1)^2))
      if (is.na(jj) || jj <= 0) prior <- 1 else prior <- sqrt(jj)*nu^2
      -(fit$loglik + log(prior))
    }
  }
  out <- optimize_omega(omega_init, obj_neg)
  list(omega = out$omega, nu = ifelse(out$omega > 0, 1/out$omega, Inf),
       convergence = out$convergence, value = out$value)
}

# FUNCTION: Pseudo Bayes estimator
# This function computes the Pseudo Bayes estimator for the degrees of freedom
# parameter nu via optimization over omega = 1/nu.
#
# INPUT: y Numeric vector. The response variable.
# INPUT:  x Numeric matrix. The design matrix (should include an intercept column).
# INPUT:  omega_init Numeric. Initial value for omega. If simulated t-error data is used, setting this to the true value of 1/nu can improve convergence speed.
# OUTPUT: omega The optimized value of omega.
# OUTPUT: nu The corresponding value of nu = 1/omega.
# OUTPUT: convergence An integer convergence code from optim. A value of 0 indicates successful convergence.}

estimate_nu_nu_block <- function(y, x, omega_init = 1/2) {
  n <- length(y)
  p <- ncol(x)
  obj_neg <- function(omega) {
    if (omega < 0) return(Inf)
    if (omega == 0) {
      fit <- fit_profile_mle_fixed_nu(y, x, nu = Inf)
      -(fit$loglik)
    } else {
      nu <- 1/omega
      fit <- fit_profile_mle_fixed_nu(y, x, nu = nu)
      # Nu fisher block prior
      jj <- (trigamma(nu/2) - trigamma(nu/2 + 1/2) - 2*(nu + 5)/(nu*(nu+1)*(nu+3)))
      if (is.na(jj) || jj <= 0) prior <- 1 else prior <- sqrt(jj) * nu^2
      -(fit$loglik + log(prior))
    }
  }
  out <- optimize_omega(omega_init, obj_neg)
  list(omega = out$omega, nu = ifelse(out$omega > 0, 1/out$omega, Inf), convergence = out$convergence, value = out$value)
}
