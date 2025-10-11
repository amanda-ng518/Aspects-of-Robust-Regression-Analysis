options(digits = 10)
RNGkind(sample.kind = "Rounding") #to keep the set.seed() consistent

#------------- Write in parallel computing style--------------#
# Pass the seed argument from the batch script 
args <- commandArgs(TRUE) 
args <- as.array(args) 
seed <- as.numeric(args[1]) 

#------------- Parameters-----------------#
set.seed(2025)
nu_true <- 2 # c(2, 5, 10)
size <- 300 # c(300, 2500, 4500)
dim <- 1 # c(1, 2, 5, 10, 20, 40, 60, 80)
ncov <- dim-1
beta_true <- rep(0, ncov+1)
sigma_true <- 1
x <- matrix(nrow = size, ncol = ncov)
if (ncov != 0) for (i in 1:ncov) x[,i] <- rnorm(n = size, mean = 0, sd = 1)
x <- cbind(1, x)
mu_true <- x %*% beta_true

#-------------Simulation function-----------------#
all_methods_sim_omega <- function(thechosenseed){ 
  
  set.seed(thechosenseed) 
  
  #------------- Simulate data-----------------#
  error_true <- sigma_true*rt(n = size, df = nu_true)
  y <- mu_true + error_true
  
  #------------Initial guess from lm-----------#
  fit_lm <- lm(y ~ x) 
  beta_initial <- fit_lm$coefficients[-2]
  sigma_initial <- sigma(fit_lm)
  omega_initial <- 1/nu_true
  # omega_initial <- 0
  
  #----------Objective functions using omega parametrization----------#
  
  # Profile log likelihood
  profile_loglik <- function(omega){
    if (omega < 0){
      return(Inf)
    } else {
      if (omega == 0){
        treg_neg_loglik_no_nu <- function(para){
          beta <- para[1:(ncov+1)]
          sigma <- para[ncov+2]
          sigma <- abs(sigma)
          z = (y - x %*% beta)/sigma
          loglik <- sum(log(dnorm(x = z))) - size*log(sigma)
          return(-loglik)
        }
      } else {
        nu = 1/omega
        treg_neg_loglik_no_nu <- function(para){
          beta <- para[1:(ncov+1)]
          sigma <- para[ncov+2]
          sigma <- abs(sigma)
          z = (y - x %*% beta)/sigma
          loglik <- sum(log(dt(x = z, df = nu))) - size*log(sigma)
          return(-loglik)
        }
      }
      optim_loglik <- optim(par = c(beta_initial, sigma_initial), fn = treg_neg_loglik_no_nu,
                            method = "BFGS", control = list(maxit = 10000))
      beta_mle <- optim_loglik$par[1:(ncov+1)]
      sigma_mle <- abs(optim_loglik$par[ncov+2])
      loglik_nu <- -optim_loglik$value
      profile <- loglik_nu
      return(-profile)
    } 
  }
  
  # Adjusted profile log likelihood
  adj_profile_loglik <- function(omega){
    if (omega < 0){
      return(Inf)
    } else {
      if (omega == 0){
        treg_neg_loglik_no_nu <- function(para){
          beta <- para[1:(ncov+1)]
          sigma <- para[ncov+2]
          sigma <- abs(sigma)
          z = (y - x %*% beta)/sigma
          loglik <- sum(log(dnorm(x = z))) - size*log(sigma)
          return(-loglik)
        }
        optim_loglik <- optim(par = c(beta_initial, sigma_initial), fn = treg_neg_loglik_no_nu,
                              method = "BFGS", control = list(maxit = 10000))
        beta_mle <- optim_loglik$par[1:(ncov+1)]
        sigma_mle <- abs(optim_loglik$par[ncov+2])
        loglik_nu <- -optim_loglik$value
        # adjustment term
        j_11 <- 0
        for (ii in 1:size) j_11 <- j_11 + x[ii,] %*% t(x[ii,])
        j_11 <- j_11/sigma_mle^2
        j_12 <- matrix(0, ncol=1, nrow=nrow(j_11))
        j_22 <- 2*size/sigma_mle^2
        j_mat <- cbind(rbind(j_11, t(j_12)), rbind(matrix(j_12,ncol=1), j_22))
      } else {
        nu = 1/omega
        treg_neg_loglik_no_nu <- function(para){
          beta <- para[1:(ncov+1)]
          sigma <- para[ncov+2]
          sigma <- abs(sigma)
          z = (y - x %*% beta)/sigma
          loglik <- sum(log(dt(x = z, df = nu))) - size*log(sigma)
          return(-loglik)
        }
        optim_loglik <- optim(par = c(beta_initial, sigma_initial), fn = treg_neg_loglik_no_nu,
                              method = "BFGS", control = list(maxit = 10000))
        beta_mle <- optim_loglik$par[1:(ncov+1)]
        sigma_mle <- abs(optim_loglik$par[ncov+2])
        loglik_nu <- -optim_loglik$value
        # adjustment term
        j_11 <- 0
        j_12 <- 0
        j_22 <- 0
        for (ii in 1:size){
          tmp1 <- (y - x %*% beta_mle)[ii]
          tmp2 <- nu*sigma_mle^2+tmp1^2
          
          j_11 <- j_11 + (x[ii,] %*% t(x[ii,]))*(tmp2-2*tmp1^2)/tmp2^2
          j_12 <- j_12 + x[ii,]*2*nu*sigma_mle*tmp1/tmp2^2
          #j_22 <- j_22 + (1/sigma_mle^2)*tmp1^2/tmp2 + 2*nu*tmp1^2/tmp2^2
          j_22 <- j_22 + 2*nu*tmp1^2/(sigma_mle^2 *nu + tmp1^2)^2 # simplified
        }
        j_11 <- (nu+1)*j_11
        j_12 <- (nu+1)*j_12
        #j_22 <- -size/sigma_mle^2 + (nu+1)*j_22
        j_22 <- (nu+1)*j_22 # simplified
        j_mat <- cbind(rbind(j_11, t(j_12)), rbind(matrix(j_12,ncol=1), j_22))
      }
      adjusted_profile <- loglik_nu - log(abs(det(j_mat)))/2
      return(-adjusted_profile)
    } 
  }
  
  # Independent Jeffrey's 
  profile_loglik_IJ <- function(omega){
    if (omega < 0){
      return(Inf)
    } else {
      if (omega == 0){
        treg_neg_loglik_no_nu <- function(para){
          beta <- para[1:(ncov+1)]
          sigma <- para[ncov+2]
          sigma <- abs(sigma)
          z = (y - x %*% beta)/sigma
          loglik <- sum(log(dnorm(x = z))) - size*log(sigma)
          loglik <- loglik - log(sigma)
          return(-loglik) 
        }
        optim_loglik <- optim(par = c(beta_initial, sigma_initial), fn = treg_neg_loglik_no_nu,
                              method = "BFGS", control = list(maxit = 10000))
        beta_mle <- optim_loglik$par[1:(ncov+1)]
        sigma_mle <- abs(optim_loglik$par[ncov+2])
        loglik_nu <- -optim_loglik$value
        Jeffreys_prior <- 1
      } else {
        nu = 1/omega
        treg_neg_loglik_no_nu <- function(para){
          beta <- para[1:(ncov+1)]
          sigma <- para[ncov+2]
          sigma <- abs(sigma)
          z = (y - x %*% beta)/sigma
          loglik <- sum(log(dt(x = z, df = nu))) - size*log(sigma)
          loglik <- loglik - log(sigma)
          return(-loglik)
        }
        optim_loglik <- optim(par = c(beta_initial, sigma_initial), fn = treg_neg_loglik_no_nu,
                              method = "BFGS", control = list(maxit = 10000))
        beta_mle <- optim_loglik$par[1:(ncov+1)]
        sigma_mle <- abs(optim_loglik$par[ncov+2])
        loglik_nu <- -optim_loglik$value
        Jeffreys_prior <- sqrt((nu/(nu+3))*(trigamma(nu/2)-trigamma(nu/2+1/2)-2*(nu+3)/(nu*(nu+1)^2)))
        Jeffreys_prior <- Jeffreys_prior*nu^2
      }
      joint_loglik <- loglik_nu + 1*log(Jeffreys_prior)
      return(-joint_loglik)
    }
  }
  
  # Marginal independent Jeffrey's
  profile_loglik_mar_IJ <- function(omega){
    if (omega < 0){
      return(Inf)
    } else {
      if (omega == 0){
        treg_neg_loglik_no_nu <- function(para){
          beta <- para[1:(ncov+1)]
          sigma <- para[ncov+2]
          sigma <- abs(sigma)
          z = (y - x %*% beta)/sigma
          loglik <- sum(log(dnorm(x = z))) - size*log(sigma)
          return(-loglik) 
        }
        optim_loglik <- optim(par = c(beta_initial, sigma_initial), fn = treg_neg_loglik_no_nu,
                              method = "BFGS", control = list(maxit = 10000))
        beta_mle <- optim_loglik$par[1:(ncov+1)]
        sigma_mle <- abs(optim_loglik$par[ncov+2])
        loglik_nu <- -optim_loglik$value
        Jeffreys_prior <- 1
      } else {
        nu = 1/omega
        treg_neg_loglik_no_nu <- function(para){
          beta <- para[1:(ncov+1)]
          sigma <- para[ncov+2]
          sigma <- abs(sigma)
          z = (y - x %*% beta)/sigma
          loglik <- sum(log(dt(x = z, df = nu))) - size*log(sigma)
          return(-loglik)
        }
        optim_loglik <- optim(par = c(beta_initial, sigma_initial), fn = treg_neg_loglik_no_nu,
                              method = "BFGS", control = list(maxit = 10000))
        beta_mle <- optim_loglik$par[1:(ncov+1)]
        sigma_mle <- abs(optim_loglik$par[ncov+2])
        loglik_nu <- -optim_loglik$value
        Jeffreys_prior <- sqrt((nu/(nu+3))*(trigamma(nu/2)-trigamma(nu/2+1/2)-2*(nu+3)/(nu*(nu+1)^2)))
        Jeffreys_prior <- Jeffreys_prior*nu^2
      }
      joint_loglik <- loglik_nu + 1*log(Jeffreys_prior)
      return(-joint_loglik)
    }
  }
  
  # Marginal Fisher
  profile_loglik_nu_block <- function(omega){
    if (omega < 0){
      return(Inf)
    } else {
      if (omega == 0){
        treg_neg_loglik_no_nu <- function(para){
          beta <- para[1:(ncov+1)]
          sigma <- para[ncov+2]
          sigma <- abs(sigma)
          z = (y - x %*% beta)/sigma
          loglik <- sum(log(dnorm(x = z))) - size*log(sigma)
          return(-loglik) 
        }
        optim_loglik <- optim(par = c(beta_initial, sigma_initial), fn = treg_neg_loglik_no_nu,
                              method = "BFGS", control = list(maxit = 10000))
        beta_mle <- optim_loglik$par[1:(ncov+1)]
        sigma_mle <- abs(optim_loglik$par[ncov+2])
        loglik_nu <- -optim_loglik$value
        Jeffreys_prior <- 1
      } else {
        nu = 1/omega
        treg_neg_loglik_no_nu <- function(para){
          beta <- para[1:(ncov+1)]
          sigma <- para[ncov+2]
          sigma <- abs(sigma)
          z = (y - x %*% beta)/sigma
          loglik <- sum(log(dt(x = z, df = nu))) - size*log(sigma)
          return(-loglik)
        }
        optim_loglik <- optim(par = c(beta_initial, sigma_initial), fn = treg_neg_loglik_no_nu,
                              method = "BFGS", control = list(maxit = 10000))
        beta_mle <- optim_loglik$par[1:(ncov+1)]
        sigma_mle <- abs(optim_loglik$par[ncov+2])
        loglik_nu <- -optim_loglik$value
        Jeffreys_prior <- sqrt(trigamma(nu/2)-trigamma(nu/2+1/2)-2*(nu+5)/(nu*(nu+1)*(nu+3)))
        Jeffreys_prior <- Jeffreys_prior*nu^2
      }
      joint_loglik <- loglik_nu + 1*log(Jeffreys_prior)
      return(-joint_loglik)
    }
  }
  
  
  #------------------Optimization using BFGS-----------------#
  # optim: run optimization with BFGS for at most 10000 iterations
  # omega: return optimized omega
  # success: indicator of successful convergence
  # sanity: check if optimized omega is at the boundary (0 or Inf)
  
  begin_time <- Sys.time()
  optim_profile <- optim(par = omega_initial, fn = profile_loglik,
                         method = "BFGS", control = list(maxit = 10000))
  omega_mle <- optim_profile$par
  success_mle <- optim_profile$convergence
  sanity_mle <- as.numeric(omega_mle == 0 || omega_mle == Inf)
  end_time <- Sys.time()
  print(end_time-begin_time)
  
  begin_time <- Sys.time()
  optim_adj_profile <- optim(par = omega_initial, fn = adj_profile_loglik,
                             method = "BFGS", control = list(maxit = 10000))
  omega_adj <- optim_adj_profile$par
  success_adj <- optim_adj_profile$convergence
  sanity_adj <- as.numeric(omega_adj == 0 || omega_adj == Inf)
  end_time <- Sys.time()
  print(end_time-begin_time)
  
  begin_time <- Sys.time()
  optim_profile_IJ <- optim(par = omega_initial, fn = profile_loglik_IJ,
                            method = "BFGS", control = list(maxit = 10000))
  omega_IJ <- optim_profile_IJ$par
  success_IJ <- optim_profile_IJ$convergence
  sanity_IJ <- as.numeric(omega_IJ == 0 || omega_IJ == Inf)
  end_time <- Sys.time()
  print(end_time-begin_time)
  
  begin_time <- Sys.time()
  optim_profile_mar_IJ <- optim(par = omega_initial, fn = profile_loglik_mar_IJ,
                                method = "BFGS", control = list(maxit = 10000))
  omega_mar_IJ <- optim_profile_mar_IJ$par
  success_mar_IJ <- optim_profile_mar_IJ$convergence
  sanity_mar_IJ <- as.numeric(omega_mar_IJ == 0 || omega_mar_IJ == Inf)
  end_time <- Sys.time()
  print(end_time-begin_time)
  
  begin_time <- Sys.time()
  optim_profile_nu_block <- optim(par = omega_initial, fn = profile_loglik_nu_block,
                                  method = "BFGS", control = list(maxit = 10000))
  omega_nu_block <- optim_profile_nu_block$par
  success_nu_block <- optim_profile_nu_block$convergence
  sanity_nu_block <- as.numeric(omega_nu_block == 0 || omega_nu_block == Inf)
  end_time <- Sys.time()
  print(end_time-begin_time)
  
  return(list(omega_mle = omega_mle, success_mle = success_mle, sanity_mle = sanity_mle,
              omega_adj = omega_adj, success_adj = success_adj, sanity_adj = sanity_adj,
              omega_IJ = omega_IJ, success_IJ = success_IJ, sanity_IJ = sanity_IJ,
              omega_mar_IJ = omega_mar_IJ, success_mar_IJ = success_mar_IJ, sanity_mar_IJ = sanity_mar_IJ,
              omega_nu_block = omega_nu_block, success_nu_block = success_nu_block, sanity_nu_block = sanity_nu_block))
}

#---------------------Run simulations------------------#
para_result <- all_methods_sim_omega(seed) 

#---------------------Save outputs------------------#
# Save outputs to a specified directory
setwd(" ")
# setwd(paste("./n", size, "_p", ncov+1, "_nu", nu_true, sep=""))
save(para_result, file=paste("allsim_", seed, ".Rda", sep = ""))  

