
library(dplyr)
library(randomForest)
library(ggplot2)
library(egg)
library(latex2exp)


# Source external functions if any
source("basic functions for TLDP.R")

source("dynamic pricing functions for_auto Scenario_1.1.R")



# Define parameters

n_P_fix <- 20000
n_Q_fix <- 10000
gamma_fix <- 1  
kappa_fix <- 1 
C_I_fix <- 1
C_r_fix <- 1/4

theta_0 <- 0.2
theta <- c(0.3, 0.2)  
tilde_theta <- -0.1
nu <- 0.1

r_star_fix <- 1/4  
tilde_p_fix <- 1/2

# Number of simulation runs
num_runs <- 100

# Case 5: Vary C_r
C_r_values <- c(0.15, 0.20, 0.25, 0.3, 0.35)


# Initialize an empty data frame to store all results in the specified format
results_summary <- data.frame(C_r = numeric(), Method = character(), value = numeric(), stringsAsFactors = FALSE)

# Loop over each parameter configuration and run multiple simulations
for (i in 1:length(C_r_values)) {
  C_r <- C_r_values[i]
  
  for (run in 1:num_runs) {
    # Run the simulation and get empirical regrets
    regrets <- run_simulation(n_P_fix, n_Q_fix, gamma_fix, kappa_fix, C_r, r_star_fix, tilde_p_fix, theta_0, theta, tilde_theta, nu, revenue_function, C_I_fix, seed = 125 + 10000 * (i + 1) + 1000 * run)
    
    # Append each method's result to the results_summary data frame
    results_summary <- rbind(results_summary,
                             data.frame(C_r = C_r, Method = 'TLDP', value = regrets[1]),
                             data.frame(C_r = C_r, Method = 'ABE', value = regrets[2]),
                             data.frame(C_r = C_r, Method = 'ExUCB', value = regrets[3]))
    
    # Print progress
    print(paste("Iteration i:", i, "Run j:", run))
  }
}

results_summary$Method <- factor(results_summary$Method, levels = c("TLDP", "ABE", "ExUCB"))





