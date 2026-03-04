# ============================================================
# Data Simulation Functions
# ============================================================

# FUNCTION: Simulate Linear Regression Data with t-Distributed Errors
# INPUT: n Integer. Number of observations to simulate. Default is 300.
# INPUT: p Integer. Number of predictors (including the intercept). Default is 1.
# INPUT: beta Numeric vector of regression coefficients of length \code{p}. The first element corresponds to the intercept.
# INPUT: sigma Numeric. Scale parameter for the t-distributed errors. Default is 1.
# INPUT: nu Numeric. Degrees of freedom for the t-distribution. Default is 2.
# INPUT: seed Optional integer. Random seed for reproducibility.
# OUTPUT:y A numeric vector of simulated response values.
# OUTPUT: x A numeric matrix of predictor values (including intercept), drawn from a standard normal distribution.

simulate_t_error_data <- function(n = 300, p = 1, beta = rep(0,1), sigma = 1, nu = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(beta) != p) stop("length(beta) must equal p")
  ncov <- p - 1
  x <- matrix(0, nrow = n, ncol = ncov)
  if (ncov > 0) for (i in 1:ncov) x[,i] <- rnorm(n)
  x <- cbind(1, x)
  mu <- x %*% beta
  errors <- sigma * rt(n = n, df = nu)
  y <- as.numeric(mu + errors)
  list(y = y, x = x)
}


# FUNCTION: Simulate Linear Regression Data with Normal Errors
# INPUT: n Integer. Number of observations to simulate. Default is 300.
# INPUT: p Integer. Number of predictors (including the intercept). Default is 1.
# INPUT: beta Numeric vector of regression coefficients of length \code{p}. The first element corresponds to the intercept.
# INPUT: mean Numeric. Mean parameter for the normal errors. Default is 0.
# INPUT: sigma Numeric. Scale parameter for the t-distributed errors. Default is 1.
# INPUT: seed Optional integer. Random seed for reproducibility.
# OUTPUT:y A numeric vector of simulated response values.
# OUTPUT: x A numeric matrix of predictor values (including intercept), drawn from a standard normal distribution.

simulate_n_error_data <- function(n = 300, p = 1, beta = rep(0,1), mean = 0, sigma = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(beta) != p) stop("length(beta) must equal p")
  ncov <- p - 1
  x <- matrix(0, nrow = n, ncol = ncov)
  if (ncov > 0) for (i in 1:ncov) x[,i] <- rnorm(n)
  x <- cbind(1, x)
  mu <- x %*% beta
  errors <- rnorm(n = n, mean = mean, sd = sigma)
  y <- as.numeric(mu + errors)
  list(y = y, x = x)
}

# Prerequisite: load stackloss data from MASS library
library(MASS)
stackloss_data = stackloss
# FUNCTION: Simulate Linear Regression Data with Stackloss data
# Data is generated such that the first 21 rows and first 3 columns are extracted from the stackloss data.
# Extra observations and columns are simulated using MVN with mean and covariance structure replicated from
# stackloss data.
#
# INPUT: n Integer. Number of observations to simulate. Default is 300.
# INPUT: p Integer. Number of predictors (including the intercept). Default is 1.
# INPUT: nu Numeric. Degrees of freedom for the t-distribution. Default is 2.
# INPUT: seed Optional integer. Random seed for reproducibility.
# OUTPUT:y A numeric vector of simulated response values.
# OUTPUT: x A numeric matrix of predictor values (including intercept), drawn from a standard normal distribution.

simulate_stackloss_data <- function(n = 21, p = 4, nu = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  data_fit_lm <- lm(as.matrix(stackloss_data[4]) ~ as.matrix(stackloss_data[1:3])) 
  true_beta <- data_fit_lm$coefficients
  sigma_true <- sigma(data_fit_lm)
  
  ncov <- p - 1
  if(p == 4){
    x <- as.matrix(stackloss_data[1:3])
  }else if (p > 4 | n > 21){
    X_orig <- as.matrix(stackloss_data[1:3])
    mu_X <- colMeans(X_orig)
    cov_X <- cov(X_orig)
    new_X <- mvrnorm(n = n - nrow(X_orig), mu = mu_X, Sigma = cov_X)
    X_simulated <- rbind(X_orig, new_X)
    x <- matrix(nrow = n, ncol = ncov)
    x[, 1:3] <- X_simulated
    for (i in 4:ncov) {
      x[, i] <- rnorm(n)  # Fill remaining columns with standard normal 
    }
    true_beta <- append(true_beta, rep(0, p-4))
  }
  
  x <- cbind(1, x)
  mu <- x %*% true_beta
  errors <- sigma_true*rt(n = n, df = nu)
  y <- mu + errors
  list(y = y, x = x)
}


# FUNCTION: Simulate Linear Regression Data with Contaminated Errors
# INPUT: n Integer. Number of observations to simulate. Default is 300.
# INPUT: p Integer. Number of predictors (including the intercept). Default is 1.
# INPUT: beta Numeric vector of regression coefficients of length \code{p}. The first element corresponds to the intercept.
# INPUT: mean Numeric. Mean parameter for the normal errors. Default is 0.
# INPUT: sigma Numeric. Scale parameter for the t-distributed errors. Default is 1.
# INPUT: contam_type Character. Type of contaminated error. Available options include: "N_0_9", "t_2", "chisq", "twopt". 
#       Set NULL (and contam_prob = 0) to ignore contamination specification.  
# INPUT: contam_prob Numeric. Proportion of data with contaminated error. Default is 0.
# INPUT: seed Optional integer. Random seed for reproducibility.
# OUTPUT:y A numeric vector of simulated response values.
# OUTPUT: x A numeric matrix of predictor values (including intercept), drawn from a standard normal distribution.

simulate_contaminated_data <- function(n = 300, p = 1, beta = rep(0, 1), mean = 0,
                                       sigma = 1, contam_type = c("N_0_9", "t_2", "chisq", "twopt", NULL), contam_prob = 0,
                                       seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  if (length(beta) != p) stop("length(beta) must equal p")
  if (is.null(contam_type) & contam_prob != 0) stop("Please specify contam_type")
  if (contam_prob < 0 | contam_prob > 1) stop("contam_prob must be between 0 and 1")
  
  # Design matrix
  ncov <- p - 1
  x <- matrix(0, nrow = n, ncol = ncov)
  if (ncov > 0) for (i in 1:ncov) x[, i] <- rnorm(n)
  x <- cbind(1, x)
  
  mu <- as.numeric(x %*% beta)
  
  # Initialize errors
  errors <- numeric(n)
  
  # Number of contaminated points
  m <- floor(n * contam_prob)
  
  # Contamination types
  # No contamination
  if (is.null(contam_type)){
    errors <- rnorm(n, mean = mean, sd = sigma)
  }
  # N(0,9)
  else if (contam_type == "N_0_9") {
    errors[1:m] <- rnorm(m, mean = 0, sd = 3)
    errors[(m+1):n] <- rnorm(n - m, mean = mean, sd = sigma)
    
    # t(2)
  } else if (contam_type == "t_2") {
    errors[1:m] <- sigma * rt(m, df = 2)
    errors[(m+1):n] <- rnorm(n - m, mean = mean, sd = sigma)
    
    # chisq(4) - 4
  } else if (contam_type == "chisq") {
    errors[1:m] <- rchisq(m, df = 4) - 4
    errors[(m+1):n] <- rnorm(n - m, mean = mean, sd = sigma)
    
    # two-point contamination
  } else if (contam_type == "twopt") {
    half <- floor(m / 2)
    y <- mu
    y[1:half] <- 5
    y[(half+1):m] <- -5
    y[(m+1):n] <- mu[(m+1):n] + rnorm(n - m, mean = mean, sd = sigma)
    return(list(y = y, x = x))
    
  } else {
    stop("Please specify correct contiminated error, or leave it as NULL")
  }
  
  y <- mu + errors
  list(y = y, x = x)
}
