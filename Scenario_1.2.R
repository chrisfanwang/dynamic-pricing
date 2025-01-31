
library(dplyr)
library(randomForest)
library(ggplot2)
library(egg)
library(latex2exp)


# Source external functions if any
source("basic functions for TLDP.R")

source("dynamic pricing functions for_auto Scenario_1.2.R")




# Define parameters


n_P_fix <- 20000
n_Q_fix <- 10000
gamma_fix <- 1  
kappa_fix <- 1 
C_I_fix <- 1
C_r_fix <- 1/4

theta_0 <- 0.3
theta <- c(0.1, 0.4, -0.2)  
tilde_theta <- - 0.1
nu <- 0.1

r_star_fix <- 1/4  
tilde_p_fix <- 1/2
# Number of simulation runs
num_runs <- 100

# Case 1: Vary kappa
kappa_values <- c(0.2, 0.4, 0.6, 0.8, 1)

# Initialize an empty data frame to store all results in the specified format
results_summary <- data.frame(kappa = numeric(), Method = character(), value = numeric(), stringsAsFactors = FALSE)

# Loop over each parameter configuration and run multiple simulations
for (i in 1:length(kappa_values)) {
  kappa <- kappa_values[i]
  
  for (run in 1:num_runs) {
    # Run the simulation and get empirical regrets
    regrets <- run_simulation(n_P_fix, n_Q_fix, gamma_fix, kappa, C_r_fix, r_star_fix, tilde_p_fix, theta_0, theta, tilde_theta, nu, revenue_function, C_I_fix, seed = 124 + 10000 * (i + 1) + 1000 * run)

    # Append each method's result to the results_summary data frame
    results_summary <- rbind(results_summary,
                             data.frame(kappa = kappa, Method = 'TLDP', value = regrets))
    
    # Print progress
    print(paste("Iteration i:", i, "Run j:", run))
  }
}

results_summary$Method <- factor(results_summary$Method, levels = c("TLDP"))





