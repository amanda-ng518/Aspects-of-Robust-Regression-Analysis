# Read the output downloaded from Compute Canada 
setwd(" ") # Set to data dir

#-------------Analyze the simulation output under a single dimension-------------# 

# Parameters
nu_true <- 2 # c(2, 5, 10)
omega_true <- 1/nu_true
size <- 300 # c(300, 2500, 4500)
dim <- 1 # c(1, 2, 5, 10, 20, 40, 60, 80)
ncov <- dim-1

# Set to data dir for this setting
setwd(paste("./n", size, "_p", ncov+1, "_nu", nu_true, sep="")) 

Nsim <- 500
success <- numeric(Nsim)
omega_mle <- c()
success_mle <- c()
omega_adj <- c()
success_adj <- c()
omega_IJ <- c()
success_IJ <- c()
omega_mar_IJ <- c() 
success_mar_IJ <- c()
omega_nu_block <- c()
success_nu_block <- c()

sanity_mle <- c()
sanity_adj <- c()
sanity_IJ <- c()
sanity_mar_IJ <- c()
sanity_nu_block <- c()

nu_mle <- c()
nu_adj <- c()
nu_IJ <- c()
nu_mar_IJ <- c()
nu_nu_block <-c()

for (k in 1:Nsim){ ## k is the index for replications 
  
  # Read the output with Nsim number of replications 
  possibleError <- tryCatch(
    # Read the output with Nsim number of replications 
    load(paste("allsim_", k, ".Rda", sep = "")),
    error=function(e) e)
  
  if(inherits(possibleError, "error")){
    success[k] <- 0 # failure: fill in with 0 
  } else{
    success[k] <- 1 # success: fill in with 1
    
    # Omega scale
    omega_mle <- c(omega_mle, para_result$omega_mle)
    success_mle <- c(success_mle, para_result$success_mle)
    sanity_mle <-c(sanity_mle, para_result$sanity_mle)
    omega_adj <- c(omega_adj, para_result$omega_adj)
    success_adj <- c(success_adj, para_result$success_adj)
    sanity_adj <-c(sanity_adj, para_result$sanity_adj)
    omega_IJ <- c(omega_IJ, para_result$omega_IJ)
    success_IJ <- c(success_IJ, para_result$success_IJ)
    sanity_IJ <-c(sanity_IJ, para_result$sanity_IJ)
    omega_mar_IJ <- c(omega_mar_IJ, para_result$omega_mar_IJ)
    success_mar_IJ <- c(success_mar_IJ, para_result$success_mar_IJ)
    sanity_mar_IJ <-c(sanity_mar_IJ, para_result$sanity_mar_IJ)
    omega_nu_block <- c(omega_nu_block, para_result$omega_nu_block)
    success_nu_block <- c(success_nu_block, para_result$success_nu_block)
    sanity_nu_block <-c(sanity_nu_block, para_result$sanity_nu_block)
    
    # Nu scale
    nu_mle = 1/omega_mle
    nu_adj = 1/omega_adj
    nu_IJ = 1/omega_IJ
    nu_mar_IJ = 1/omega_mar_IJ
    nu_nu_block = 1/omega_nu_block
    
  }
} 

# Convergence check
sum(success)
sum(success_mle)
sum(success_adj)
sum(success_IJ)
sum(success_mar_IJ)
sum(success_nu_block)

# Sanity check
sum(sanity_mle)
sum(sanity_adj)
sum(sanity_IJ)
sum(sanity_mar_IJ)
sum(sanity_nu_block)

# Range of returned values (nu scale)
range(nu_mle)
range(nu_adj)
range(nu_IJ)
range(nu_mar_IJ)
range(nu_nu_block)

# Nu bias
bias_nu_mle = mean(nu_mle) - nu_true
bias_nu_adj = mean(nu_adj) - nu_true
bias_nu_IJ = mean(nu_IJ) - nu_true
bias_nu_mar_IJ = mean(nu_mar_IJ) - nu_true
bias_nu_nu_block = mean(nu_nu_block) - nu_true

# Nu se
se_nu_mle <- sd(nu_mle)
se_nu_adj <- sd(nu_adj)
se_nu_IJ <- sd(nu_IJ)
se_nu_mar_IJ <- sd(nu_mar_IJ)
se_nu_nu_block <- sd(nu_nu_block)

# Rmse (nu scale)
sqrt(mean((nu_mle-nu_true)^2))
sqrt(mean((nu_adj-nu_true)^2))
sqrt(mean((nu_IJ-nu_true)^2))
sqrt(mean((nu_mar_IJ-nu_true)^2))
sqrt(mean((nu_nu_block-nu_true)^2))

#--------------------RMSE Boxplot for a single dimension--------------------#
plot_dir <- " "  # Set to your plot dir
file_name <- paste("n", size, "_p", ncov+1, "_nu", nu_true, "boxplot.png", sep="")
full_path <- file.path(plot_dir, file_name)
png(full_path, width = 800, height = 600)

boxplot(nu_mle, nu_adj, nu_IJ, nu_mar_IJ, nu_nu_block,
        names=c("MLE", "Adj", "IJ", "Mar IJ", "Mar Fisher"),
        ylab="Nu", main=paste("n=", size, ", p=", ncov+1, ", nu=", nu_true, sep=""))
abline(h=nu_true, col="red", lty=2)

dev.off()


#-------------Analyze the simulation output across dimensions-------------# 
# Set your data dir
setwd(" ")

# Parameters
nu_true <- 2 # c(2, 5, 10)
omega_true <- 1/nu_true
size <-300 # c(300, 2500, 4500)
ncovs <- c(1,2,5,10,20,40,60,80)-1
rmse_nu_mle <- c()
rmse_nu_adj <- c()
rmse_nu_IJ <- c()
rmse_nu_mar_IJ <- c()
rmse_nu_nu_block <- c()

for (ncov in ncovs){
  
  i = which(ncovs == ncov)
  
  # Set to data dir for this setting
  setwd(paste("./n", size, "_p", ncov+1, "_nu", nu_true, sep="")) 
  
  Nsim <- 500
  
  omega_mle <- c()
  omega_adj <- c()
  omega_IJ <- c()
  omega_mar_IJ <- c()
  omega_nu_block <- c()
  
  for (k in 1:Nsim){ ## k is the index for replications
    load(paste("allsim_", k, ".Rda", sep = ""))
    
    omega_mle <- c(omega_mle, para_result$omega_mle)
    omega_adj <- c(omega_adj, para_result$omega_adj)
    omega_IJ <- c(omega_IJ, para_result$omega_IJ)
    omega_mar_IJ <- c(omega_mar_IJ, para_result$omega_mar_IJ)
    omega_nu_block <- c(omega_nu_block, para_result$omega_nu_block)
  }
  rmse_nu_mle <- c(rmse_nu_mle, sqrt(mean((1/omega_mle-nu_true)^2)))
  rmse_nu_adj <- c(rmse_nu_adj, sqrt(mean((1/omega_adj-nu_true)^2)))
  rmse_nu_IJ <- c(rmse_nu_IJ, sqrt(mean((1/omega_IJ-nu_true)^2)))
  rmse_nu_mar_IJ <- c(rmse_nu_mar_IJ, sqrt(mean((1/omega_mar_IJ-nu_true)^2)))
  rmse_nu_nu_block <- c(rmse_nu_nu_block, sqrt(mean((1/omega_nu_block-nu_true)^2)))
}

# RMSE
rmse_nu_mle
rmse_nu_adj
rmse_nu_IJ
rmse_nu_mar_IJ
rmse_nu_nu_block

#--------------------RMSE line plot across dimensions----------#
plot_dir <- " "  # Set your plot dir
file_name <- paste("nu", nu_true,"rmseplot.png", sep="")
full_path <- file.path(plot_dir, file_name)
png(full_path, width = 800, height = 600)

plot(x=ncovs+1, y=rmse_nu_mle,
     ylim=range(rmse_nu_mle,rmse_nu_adj,rmse_nu_IJ,rmse_nu_mar_IJ,rmse_nu_nu_block),
     type="l", col="black",
     main=paste("n=", size, ", nu=", nu_true, sep=""), xlab="p", ylab="Nu RMSE")
lines(x=ncovs+1, y=rmse_nu_adj, col="blue")
lines(x=ncovs+1, y=rmse_nu_IJ, col="red")
lines(x=ncovs+1, y=rmse_nu_mar_IJ, col="green")
lines(x=ncovs+1, y=rmse_nu_nu_block, col="purple")
legend("topleft", title="method",
       legend=c("MLE", "Adj", "IJ", "Mar IJ", "Mar Fisher"),
       col=c("black","blue","red","green","purple"), lty=1, bty="n", xpd=T)

dev.off()

