

base_centers <- matrix(c(
  1/4, 1/4, 1/4,
  3/4, 1/4, 1/4,
  1/4, 3/4, 1/4,
  3/4, 3/4, 1/4,
  1/4, 1/4, 3/4,
  3/4, 1/4, 3/4,
  1/4, 3/4, 3/4,
  3/4, 3/4, 3/4
), ncol = 3, byrow = TRUE)



centers_df <- as.data.frame(base_centers)
p_list  <- rep(c(1/4, 3/4), each = 4)






# Define parameters

sigma <- 0.05

n_P_fix <- 20000
n_Q_fix <- 10000
gamma_fix <- 1  
kappa_fix <- 1 
C_I_fix <- 1
C_r_fix <- 1/4

r_star_fix <- 1/4  
tilde_p_fix <- 1/2

# Number of simulation runs
num_runs <- 50

# Case 1: Vary n_P
n_P_values <- c(0, 10000, 20000, 30000, 40000)

# Initialize an empty data frame to store all results in the specified format
results_summary <- data.frame(n_P = numeric(), Method = character(), value = numeric(), stringsAsFactors = FALSE)

# Loop over each parameter configuration and run multiple simulations
for (i in 1:length(n_P_values)) {
  n_P <- n_P_values[i]
  
  for (run in 1:num_runs) {
    # Run the simulation and get empirical regrets
    regrets <- run_simulation(n_P, n_Q_fix, gamma_fix, kappa_fix, C_r_fix, r_star_fix, tilde_p_fix, revenue_function, C_I_fix, sigma, seed = 121 + 10000 * (i) + 1000 * run)
    
    # Append each method's result to the results_summary data frame
    results_summary <- rbind(results_summary,
                             data.frame(n_P = n_P, Method = 'TLDP', value = regrets[1]),
                             data.frame(n_P = n_P, Method = 'ABE', value = regrets[2]),
                             data.frame(n_P = n_P, Method = 'ExUCB', value = regrets[3]))
    
    # Print progress
    print(paste("Iteration i:", i, "Run j:", run))
  }
}


results_summary$Method <- factor(results_summary$Method, levels = c("TLDP", "ABE", "ExUCB"))






