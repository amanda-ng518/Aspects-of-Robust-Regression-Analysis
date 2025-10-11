options(digits = 10)
RNGkind(sample.kind = "Rounding") #to keep the set.seed() consistent

#------------- Write in parallel computing style--------------#
# Pass the seed argument from the batch script 
args <- commandArgs(TRUE) 
args <- as.array(args) 
seed <- as.numeric(args[1]) 


#-----------------------Load library--------------#
library(doParallel)
library(foreach)
library(rlmDataDriven)
library(MASS)

#------------- Parameters-----------------#
set.seed(2025)
size <- 20 # c(20, 500, 100)
ncov <- 3 # c(3, 3, 49)
p <- ncov + 1 
nu_true <- 2
true_beta <- rep(0, p)
sigma_true <- 1

#-------------Generate design matrix X-----------------#
if(p == 1){
  x <- matrix(data = 1, nrow = size, ncol = 1)
}else{
  x <- matrix(nrow = size, ncol = ncov)
  for (i in 1:ncov) x[,i] <- rnorm(n = size, mean = 0, sd = 1)
  x <- cbind(1, x)
}
mu_true <- x %*% true_beta

#-------------Simulation function-----------------#
sim_omegaparam_beta <- function(thechosenseed){ 
  
  set.seed(thechosenseed) 
  
  #------------- Simulate data-----------------#
  # t(2) error
  error_true <- sigma_true*rt(n = size, df = nu_true)
  
  # N(0,1) error
  #error_true <- rnorm(n = size, mean = 0, sd = sigma_true)
  
  y <- mu_true + error_true
  
  #------------Initial guess from lm-----------#
  fit_lm <- lm(y ~ x) 
  beta_initial <- fit_lm$coefficients[-2]
  sigma_initial <- sigma(fit_lm)
  omega_initial <- 1/nu_true
  
  #------------Huber regression-----------#
  RB <- rlm(y ~ x[,-1], psi = psi.huber, k = 1.345, maxit = 100)
  beta_huber <- rlmDD(y, x[,-1], fit_lm$coef[-2], RB$coef, method = "Huber", plot = "N")$esti[["coefficients"]]
  
  #----------Objective functions using omega parametrization----------#
  
  # Profile log likelihood
  profile_loglik <- function(omega){
    if (omega < 0){
      return(1e10)
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
      return(1e10)
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
          j_22 <- j_22 + (1/sigma_mle^2)*tmp1^2/tmp2 + 2*nu*tmp1^2/tmp2^2
          #j_22 <- j_22 + 2*nu*tmp1^2/(sigma_mle^2 *nu + tmp1^2)^2 # simplified
        }
        j_11 <- (nu+1)*j_11
        j_12 <- (nu+1)*j_12
        j_22 <- -size/sigma_mle^2 + (nu+1)*j_22
        #j_22 <- (nu+1)*j_22 # simplified
        j_mat <- cbind(rbind(j_11, t(j_12)), rbind(matrix(j_12,ncol=1), j_22))
      }
      adjusted_profile <- loglik_nu - log(abs(det(j_mat)))/2
      return(-adjusted_profile)
    } 
  }
  
  # Independent Jeffrey's 
  profile_loglik_IJ <- function(omega){
    if (omega < 0){
      return(1e10)
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
  
  # Log likelihood
  treg_neg_loglik<- function(para, omega){
    
    beta <- para[1:(ncov+1)]
    sigma <- abs(para[ncov+2])
    z <- (y - x %*% beta)/sigma
    if (omega == 0){
      loglik <- sum(log(dnorm(x = z, mean = 0, sd = 1))) - size*log(sigma)                 
    }else if (omega < 0){
      return (Inf)
    }else{
      loglik <- sum(log(dt(x = z, df = 1/omega))) - size*log(sigma)
    }
    return(-loglik)
  }
  
  #------------------Optimization using BFGS-----------------#
  # optim_nu: run optimization with BFGS on nu 
  # omega: return optimized omega
  # success: indicator of successful convergence
  # sanity: check if optimized omega is at the boundary (0 or Inf)
  # optim_beta: run optimization with BFGS on beta with nu fixed at optim_nu
  # beta: return optimized beta
  
  # Profile
  begin_time <- Sys.time()
    optim_mle_profile <- optim(par = omega_initial, fn = profile_loglik,
                               method = "BFGS", control = list(maxit = 10000))
    omega_mle <- optim_mle_profile$par
    success_mle <- optim_mle_profile$convergence
    sanity_mle <- as.numeric(omega_mle == 0 || omega_mle == 1e10)
    mle_optim_loglik_beta <- optim(par = c(beta_initial, sigma_initial), fn = treg_neg_loglik,
                                   omega = omega_mle, 
                                   method = "BFGS", control = list(maxit = 10000))
    beta_mle<- mle_optim_loglik_beta$par[2:(ncov+1)]
    beta_success_mle <- mle_optim_loglik_beta$convergence
  end_time <- Sys.time()
  print(end_time-begin_time)
  
  # Adj
  begin_time <- Sys.time()
    optim_adj_profile <- optim(par = omega_initial, fn = adj_profile_loglik,
                               method = "BFGS", control = list(maxit = 10000))
    omega_adj <- optim_adj_profile$par
    success_adj <- optim_adj_profile$convergence
    sanity_adj <- as.numeric(omega_adj == 0 || omega_adj == 1e10)
    adj_optim_loglik_beta <- optim(par = c(beta_initial, sigma_initial), fn = treg_neg_loglik,
                                   omega = omega_adj, 
                                   method = "BFGS", control = list(maxit = 10000))
    beta_adj <- adj_optim_loglik_beta$par[2:(ncov+1)]
    beta_success_adj <- adj_optim_loglik_beta$convergence
  end_time <- Sys.time()
  print(end_time-begin_time)
  
  # IJ
  begin_time <- Sys.time()
    optim_profile_IJ <- optim(par = omega_initial, fn = profile_loglik_IJ,
                              method = "BFGS", control = list(maxit = 10000))
    omega_IJ <- optim_profile_IJ$par
    success_IJ <- optim_profile_IJ$convergence
    sanity_IJ <- as.numeric(omega_IJ == 0 || omega_IJ == 1e10)
    IJ_optim_loglik_beta <- optim(par = c(beta_initial, sigma_initial), fn = treg_neg_loglik,
                                  omega = omega_IJ, 
                                  method = "BFGS", control = list(maxit = 10000))
    beta_IJ <- IJ_optim_loglik_beta$par[2:(ncov+1)]
    beta_success_IJ <- IJ_optim_loglik_beta$convergence
  end_time <- Sys.time()
  print(end_time-begin_time)
  
  # True nu
  begin_time <- Sys.time()
  true_nu = nu_true
  true_optim_loglik_beta <- optim(par = c(beta_initial, sigma_initial), fn = treg_neg_loglik,
                                  omega = 1/true_nu, 
                                  method = "BFGS", control = list(maxit = 10000))
  beta_true <- true_optim_loglik_beta$par[2:(ncov+1)]
  beta_success_true <- true_optim_loglik_beta$convergence
  end_time <- Sys.time()
  print(end_time-begin_time)
  
  # Fixed nu = 1:10, 20, 30
  begin_time <- Sys.time()
  for (fixed_nu in c(1:10, 20, 30)){
    fixed_optim_loglik_beta <- optim(par = c(beta_initial, sigma_initial), fn = treg_neg_loglik,
                                     omega = 1/fixed_nu, 
                                     method = "BFGS", control = list(maxit = 10000))
    assign(paste0("beta_fixed_", fixed_nu, sep = ""), fixed_optim_loglik_beta$par[2:(ncov+1)])
    assign(paste0("beta_success_fixed_", fixed_nu, sep = ""), fixed_optim_loglik_beta$convergence)
  }
  end_time <- Sys.time()
  print(end_time-begin_time)
  
  return(list(nu_adj = 1/omega_adj, success_adj = success_adj, sanity_adj = sanity_adj, 
              beta_adj = beta_adj, beta_success_adj = beta_success_adj,
              nu_mle = 1/omega_mle, success_mle = success_mle, sanity_mle = sanity_mle, 
              beta_mle = beta_mle, beta_success_mle = beta_success_mle,
              nu_IJ = 1/omega_IJ, success_IJ = success_IJ, sanity_IJ = sanity_IJ,
              beta_IJ = beta_IJ, beta_success_IJ = beta_success_IJ,
              beta_true = beta_true, beta_success_true= beta_success_true,
              beta_fixed_1 = beta_fixed_1, beta_success_fixed_1 = beta_success_fixed_1,
              beta_fixed_2 = beta_fixed_2, beta_success_fixed_2 = beta_success_fixed_2,
              beta_fixed_3 = beta_fixed_3, beta_success_fixed_3 = beta_success_fixed_3,
              beta_fixed_4 = beta_fixed_4, beta_success_fixed_4 = beta_success_fixed_4,
              beta_fixed_5 = beta_fixed_5, beta_success_fixed_5 = beta_success_fixed_5,
              beta_fixed_6 = beta_fixed_6, beta_success_fixed_6 = beta_success_fixed_6,
              beta_fixed_7 = beta_fixed_7, beta_success_fixed_7 = beta_success_fixed_7,
              beta_fixed_8 = beta_fixed_8, beta_success_fixed_8 = beta_success_fixed_8,
              beta_fixed_9 = beta_fixed_9, beta_success_fixed_9 = beta_success_fixed_9,
              beta_fixed_10 = beta_fixed_10, beta_success_fixed_10 = beta_success_fixed_10,
              beta_fixed_20 = beta_fixed_20, beta_success_fixed_20 = beta_success_fixed_20,
              beta_fixed_30 = beta_fixed_30, beta_success_fixed_30 = beta_success_fixed_30,
              beta_ols = beta_initial[2:(ncov+1)],
              beta_huber = beta_huber[2:(ncov+1)]))
}

#---------------------Run simulations------------------#
result <- sim_omegaparam_beta(seed) 

#---------------------Save outputs------------------#
# Save outputs to a specified directory
setwd(" ")
# setwd(paste("./n", size, "_p", ncov+1, "_nu", nu_true, sep=""))
save(result, file=paste("allsim_", seed, ".Rda", sep = ""))   
