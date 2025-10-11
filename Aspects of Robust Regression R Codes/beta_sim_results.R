#-----------------------Load library--------------#
library(ggplot2)

#---------------------Simulation setting--------------#
size = 100 
p <- 3
ncov <- p - 1
Nsim <- 500
error <- "no contamination" # no contamination, N(0, 9), chi^2(4) − 4, t(2), two-point contamination

#---------------------True beta-----------------#
# Completely simulated data
true_beta_exclude_intercept = rep(0, p - 1)

# Stackloss
# stackloss_data <- stackloss
# x <- as.matrix(stackloss_data[1:3])
# x <- cbind(1, x)
# data_fit_lm <- lm(as.matrix(stackloss_data[4]) ~ x)
# true_beta <- data_fit_lm$coefficients[-2]
# true_beta <- append(true_beta, rep(0, p - 4))
# true_beta_exclude_intercept <- true_beta[2:p]

# Contaminated error
# true_beta_exclude_intercept = c(5, 10)

#----------------------Initialize containers for each method------------------#
methods <- c(paste0("fixed_", c(1:10, 20, 30)), "OLS","Huber","mle","adj", "IJ")
bias_list <- list()
se_list <- list()
rmse_list<- list()

for (m in methods) {
  bias_list[[m]] <- matrix(NA, nrow = 1, ncol = ncov)
  se_list[[m]] <- matrix(NA, nrow = 1, ncol = ncov)
  rmse_list[[m]]<- matrix(NA, nrow = 1, ncol = ncov)
}

beta_store <- list(adj = NULL, IJ = NULL, mle = NULL, true = NULL, 
                   fixed_1 = NULL,fixed_2 = NULL,fixed_3 = NULL,fixed_4 = NULL,
                   fixed_5 = NULL,fixed_6 = NULL,fixed_7 = NULL,fixed_8 = NULL,
                   fixed_9 = NULL,fixed_10 = NULL,fixed_20 = NULL,fixed_30 = NULL,
                   OLS = NULL, Huber = NULL)

x_labels <- c(as.character(c(1:10, 20, 30)), "OLS", "Huber","Profile","ADJ", "Full IJ")

#-----------------------Analyze the simulation output------------------#
# Enter appropriate data directory
setwd(" ")

# Extract beta estimates
for (k in 1:Nsim) {
  fname <- paste0("allsim_", k, ".Rda")
  
  
  possibleError <- tryCatch(load(fname), error = function(e) e)
  
  if (!inherits(possibleError, "error")) {
    beta_store$adj   <- rbind(beta_store$adj,   result$beta_adj)
    beta_store$mle   <- rbind(beta_store$mle,   result$beta_mle)
    beta_store$IJ    <- rbind(beta_store$IJ,    result$beta_IJ)
    beta_store$true  <- rbind(beta_store$true,  result$beta_true)
    for (i in c(1:10, 20, 30)) {
      beta_store[[paste0("fixed_", i)]] <- rbind(beta_store[[paste0("fixed_", i)]], result[[paste0("beta_fixed_", i)]])
    }
    beta_store$OLS <- rbind(beta_store$OLS,  result$beta_ols)
    beta_store$Huber <- rbind(beta_store$Huber,  result$beta_huber)
  }
}

# Compute Bias, SE and RMSE for each method
for (m in methods) {
  beta_mat <- beta_store[[m]]
  if (!is.null(beta_mat)) {
    bias_list[[m]] <- colMeans(abs(beta_mat - true_beta_exclude_intercept))
    se_list[[m]]   <- apply(beta_mat, 2, sd, na.rm = TRUE)
    rmse_list[[m]] <- sqrt(colMeans((beta_mat - 
                                       matrix(true_beta_exclude_intercept, 
                                              nrow = nrow(beta_mat), 
                                              ncol = length(true_beta_exclude_intercept), 
                                              byrow = TRUE))^2))
  }
}

#---------------------Overall RMSE Plot-------------------------#
# Enter appropriate plot dir
plot_dir <- " "

# Calculate overall RMSE for each method
overall_rmse_vals <- sapply(methods, function(m) {
  sqrt(mean(rmse_list[[m]]^2))  # square all RMSEs per beta, average, then sqrt
})

# Create a data frame for plotting
df_overall_rmse <- data.frame(
  Nu = factor(x_labels, levels = x_labels),
  RMSE = overall_rmse_vals
)

# Create a flag for the est nu approaches
df_overall_rmse$color_flag <- "darkblue"
df_overall_rmse$color_flag[(nrow(df_overall_rmse)-2):nrow(df_overall_rmse)] <- "red"

# Split into two groups for line and no line
df_overall_rmse <- df_overall_rmse %>%
  mutate(group = ifelse(Nu %in% c("Huber", "Profile", "ADJ", "Full IJ"), "method", "nu_ols"))

# Create a descriptive flag instead of raw colors
df_overall_rmse$method_type <- ifelse(df_overall_rmse$color_flag == "darkblue", 
                                      "Non v-estimation approaches", "v-estimation approaches")

# Plot
p_overall_rmse <- ggplot(df_overall_rmse, aes(x = Nu, y = RMSE)) +
  # Line only for fixed nu approaches + OLS
  geom_line(data = subset(df_overall_rmse, group == "nu_ols"), 
            aes(group = 1, color = method_type), size = 1) +
  # Points with legend
  geom_point(aes(color = method_type), size = 3) +
  scale_color_manual(values = c("Non v-estimation approaches" = "darkblue", 
                                "v-estimation approaches" = "red")) +
  labs(
    x = "Fixed nu",
    y = "Overall RMSE",
    color = "Approaches"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = c(0.01, 0.99),   # top left
    legend.justification = c("left", "top"),
    legend.text  = element_text(size = 10),
    legend.title  = element_text(size = 10),
    legend.background = element_rect(fill = alpha('white', 0.6))
  )

ggsave(filename = file.path(plot_dir, paste0("overall_rmse_beta_", error,".png")),
       plot = p_overall_rmse,
       width = 8,
       height = 4,
       bg = "white")

#---------------------Overall Bias Plot-------------------------#
# Calculate overall Bias for each method
overall_bias_vals <- sapply(methods, function(m) {
  mean(bias_list[[m]])  
})

# Create a data frame for plotting
df_overall_bias <- data.frame(
  Nu = factor(x_labels, levels = x_labels),
  Bias = overall_bias_vals
)

# Create a flag for the est nu approaches
df_overall_bias$color_flag <- "darkblue"
df_overall_bias$color_flag[(nrow(df_overall_bias)-2):nrow(df_overall_bias)] <- "red"

# Split into two groups for line and no line
df_overall_bias <- df_overall_bias %>%
  mutate(group = ifelse(Nu %in% c("Huber", "Profile", "ADJ", "Full IJ"), "method", "nu_ols"))

# Create a descriptive flag instead of raw colors
df_overall_bias$method_type <- ifelse(df_overall_bias$color_flag == "darkblue", 
                                      "Non v-estimation approaches", "v-estimation approaches")

# Create the plot
p_overall_bias <- ggplot(df_overall_bias, aes(x = Nu, y = Bias, group = 1)) +
  # Line only for fixed nu approaches + OLS
  geom_line(data = subset(df_overall_bias, group == "nu_ols"), 
            aes(group = 1, color = method_type), size = 1) +
  # Points with legend
  geom_point(aes(color = method_type), size = 3) +
  scale_color_manual(values = c("Non v-estimation approaches" = "darkblue", 
                                "v-estimation approaches" = "red")) +
  labs(
    x = "Fixed nu",
    y = "Overall Absolute Bias",
    color = "Approaches"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = c(0.01, 0.99),   # top left
    legend.justification = c("left", "top"),
    legend.text  = element_text(size = 10),
    legend.title  = element_text(size = 10),
    legend.background = element_rect(fill = alpha('white', 0.6))
  )


ggsave(filename = file.path(plot_dir, paste0("overall_bias_beta_",error,".png")),
       plot = p_overall_bias,
       width = 8,
       height = 4,
       bg = "white")

#---------------------Overall SE Plot-------------------------#
# Calculate overall SE for each method
overall_se_vals <- sapply(methods, function(m) {
  mean(se_list[[m]])  
})

# Create a data frame for plotting
df_se_bias <- data.frame(
  Nu = factor(x_labels, levels = x_labels),
  SE = overall_se_vals
)

# Create a flag for the est nu approaches
df_se_bias$color_flag <- "darkblue"
df_se_bias$color_flag[(nrow(df_se_bias)-2):nrow(df_se_bias)] <- "red"

# Split into two groups for line and no line
df_se_bias <- df_se_bias %>%
  mutate(group = ifelse(Nu %in% c("Huber", "Profile", "ADJ", "Full IJ"), "method", "nu_ols"))

# Create a descriptive flag instead of raw colors
df_se_bias$method_type <- ifelse(df_se_bias$color_flag == "darkblue", 
                                      "Non v-estimation approaches", "v-estimation approaches")

# Create the plot
p_overall_se <- ggplot(df_se_bias, aes(x = Nu, y = SE, group = 1)) +
  # Line only for fixed nu approaches + OLS
  geom_line(data = subset(df_se_bias, group == "nu_ols"), 
            aes(group = 1, color = method_type), size = 1) +
  # Points with legend
  geom_point(aes(color = method_type), size = 3) +
  scale_color_manual(values = c("Non v-estimation approaches" = "darkblue", 
                                "v-estimation approaches" = "red")) +
  labs(
    x = "Fixed nu",
    y = "Overall Standard Error",
    color = "Approaches"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = c(0.01, 0.99),   # top left
    legend.justification = c("left", "top"),
    legend.text  = element_text(size = 10),
    legend.title  = element_text(size = 10),
    legend.background = element_rect(fill = alpha('white', 0.6))
  )

ggsave(filename = file.path(plot_dir, paste0("overall_se_beta_",error,".png")),
       plot = p_overall_se,
       width = 8,
       height = 4,
       bg = "white")

