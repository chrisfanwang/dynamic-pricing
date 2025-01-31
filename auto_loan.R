library(dplyr)
library(randomForest)
library(ggplot2)
library(egg)
library(latex2exp)


# Source external functions if any
source("basic functions for TLDP.R")

source("dynamic pricing functions for auto loan data.R")


# Step 1: Read and Preprocess the Data
data0 <- read.csv("CPRM_AutoLOan_OnlineAutoLoanData.csv", stringsAsFactors = TRUE)

# Convert date columns to Date type
data0$Approve_Date <- as.Date(data0$Approve_Date, format = '%m/%d/%Y')
data0$Apply_Date <- as.Date(data0$Apply_Date, format = '%m/%d/%Y')
data0$Fund_Date <- as.Date(data0$Fund_Date, format = '%m/%d/%Y')

# Order data by Apply_Date
data0 <- data0[order(data0$Apply_Date), ]

# Display the count of samples per state
print("Sample Counts by State:")
print(table(data0$State))

# Step 2: Define State to Region and Division Mapping
state_mapping <- data.frame(
  State = c("CT", "ME", "MA", "NH", "RI", "VT",
            "NJ", "NY", "PA",
            "IL", "IN", "MI", "OH", "WI",
            "IA", "KS", "MN", "MO", "NE", "ND", "SD",
            "DE", "FL", "GA", "MD", "NC", "SC", "VA", "DC", "WV",
            "AL", "KY", "MS", "TN",
            "AR", "LA", "OK", "TX",
            "AZ", "CO", "ID", "MT", "NV", "NM", "UT", "WY",
            "AK", "CA", "HI", "OR", "WA"),
  Region = c(rep("Northeast", 6),      # New England (6 states)
             rep("Northeast", 3),      # Middle Atlantic (3 states)
             rep("Midwest", 5),         # East North Central (5 states)
             rep("Midwest", 7),         # West North Central (7 states)
             rep("South", 9),           # South Atlantic (9 states)
             rep("South", 4),           # East South Central (4 states)
             rep("South", 4),           # West South Central (4 states)
             rep("West", 8),            # Mountain (8 states)
             rep("West", 5)),           # Pacific (5 states)
  Division = c(
    # Northeast
    rep("New England", 6),
    rep("Middle Atlantic", 3),
    # Midwest
    rep("East North Central", 5),
    rep("West North Central", 7),
    # South
    rep("South Atlantic", 9),
    rep("East South Central", 4),
    rep("West South Central", 4),
    # West
    rep("Mountain", 8),
    rep("Pacific", 5)
  ),
  stringsAsFactors = FALSE
)

# Step 3: Merge Mapping with data0 and Remove Unmatched States
data_merged <- data0 %>%
  left_join(state_mapping, by = "State") %>%
  filter(!is.na(Region))  # Remove records with unmatched states

# Notify about removed records
unmatched_count <- nrow(data0) - nrow(data_merged)
if (unmatched_count > 0) {
  cat("Removed", unmatched_count, "records with unmatched states.\n")
}

# Step 4: Split Data into Divisions and Assign to Separate data.frames
division_names <- unique(data_merged$Division)
for (div in division_names) {
  # Create a valid R variable name by replacing spaces with underscores
  div_var_name <- paste0("Division_", gsub(" ", "_", div))
  
  # Assign the filtered data to the new data.frame
  assign(div_var_name, filter(data_merged, Division == div))
}

# Verify that data.frames have been created
division_df_names <- paste0("Division_", gsub(" ", "_", division_names))
print("Created Division Data Frames:")
print(division_df_names)

# Step 5: Define Transformation Function
compute_price_and_normalize <- function(df, libor_rate, vars_to_normalize) {
  # Compute Price of the Loan (Premium)
  df <- df %>%
    mutate(
      Price = mp * (1 - (1 + libor_rate)^(-Term)) / libor_rate - Amount_Approved
    )
  
  # Normalize specified variables
  df <- df %>%
    mutate(across(all_of(vars_to_normalize), ~ (.-min(., na.rm = TRUE)) / (max(.) - min(., na.rm = TRUE))))
  
  # Compute y
  df <- df %>%
    mutate(
      y = apply * Price
    )
  
  # Select relevant columns and remove missing values
  df_clean <- df %>%
    select(y, all_of(vars_to_normalize)) %>%
    na.omit()
  
  return(df_clean)
}

# Define LIBOR rate and variables to normalize
libor <- 0.0012
variables_used <- c('Price', 'Primary_FICO', 'Competition_rate', 'Amount_Approved', 'onemonth', 'Term')

# Step 6: Apply Transformation to All Divisions
for (div in division_df_names) {
  # Access the division data.frame
  df <- get(div)
  
  # Apply the transformation
  transformed_df <- compute_price_and_normalize(df, libor_rate = libor, vars_to_normalize = variables_used)
  
  # Assign the transformed data.frame back with the same name
  assign(div, transformed_df)
}

# Step 7: Count Samples by Region and Division
# Count samples by Region
region_counts <- data_merged %>%
  group_by(Region) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))

cat("Sample Counts by Region:\n")
print(region_counts)

# Count samples by Division
division_counts <- data_merged %>%
  group_by(Division) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))

cat("\nSample Counts by Division:\n")
print(division_counts)

# Example: View the first few rows of Division_Pacific after transformation
if (exists("Division_Pacific")) {
  cat("\nFirst few rows of Division_Pacific after transformation:\n")
  print(head(Division_Pacific))
} else {
  cat("\nDivision_Pacific does not exist.\n")
}



target_data_set <- Division_East_South_Central

source_data_ava_set <- Division_Pacific

formula <- y ~ Price + Primary_FICO + Competition_rate + Amount_Approved + onemonth + Term

find_optimal_price_optimize <- function(observation, model) {
  # Define the function to be minimized (negative predicted y)
  objective_function <- function(price) {
    # Ensure 'price' is within bounds
    price <- min(max(price, 0), 1)
    # Create a copy of the observation
    obs <- observation
    # Set 'Price' to the current value
    obs$Price <- price
    # Predict 'y' using the model
    pred_y <- predict(model, newdata = obs)
    # Return negative of predicted y because optimize minimizes by default
    return(-as.numeric(pred_y))
  }
  
  # Use optimize function to find the price that minimizes negative predicted y
  opt_result <- optimize(f = objective_function, interval = c(0, 1), maximum = FALSE)
  # Return the optimal price
  optimal_price <- opt_result$minimum
  return(optimal_price)
}

compute_n_re_P <- function(B, D_P) {
  # B is a list with center and radius
  # D_P is the source dataset with columns X (covariate), p (price), Y (revenue)
  # Combine X and p into a matrix
  if (is.null(D_P)){
    return(list(n_P = 0, re_P = 0))
  }else{
    D_P_Xp <- as.matrix(cbind(D_P$X.1,D_P$X.2, D_P$X.3, D_P$X.4, D_P$X.5, D_P$p))
    # Compute distances using ℓ_∞-norm
    distances <- apply(D_P_Xp, 1, function(row) {
      max(abs(row - B$center))
    })
    in_ball <- distances <= B$radius
    n_P <- sum(in_ball)
    re_P <- sum(D_P$Y[in_ball])
    return(list(n_P = n_P, re_P = re_P))
  }
}

variable_covariate <- c('Primary_FICO','Competition_rate','Amount_Approved','onemonth','Term')

d <- length(variable_covariate)

# Define the simulation function
run_simulation <- function(D_P_list, target_data_set, rf_model, f_optimal, C_r, kappa, d, gamma, seed) {
  
  set.seed(seed)
  
  # Get the number of data frames in D_P_list
  D_P_size <- length(D_P_list)
  
  # Sample target data
  target_indices <- sample(1:nrow(target_data_set), size = floor(0.9 * nrow(target_data_set)))
  target_data <- target_data_set[target_indices, ]
  
  # Train random forest model
 
  
  # Prepare target covariate matrix
  X_target <- as.matrix(target_data[, variable_covariate])
  
  # Compute optimal prices using the trained model
  optimal_prices <- sapply(1:nrow(target_data), function(i) {
    obs <- target_data[i, ]
    find_optimal_price_optimize(obs, rf_model)
  })
  
  # Update the target data with the optimal prices and predictions
  target_data$Price <- optimal_prices
  p_optimal <- optimal_prices
  f_optimal <- as.vector(predict(rf_model, newdata = target_data))
  
  # Initialize constants
  gamma <- 1
  kappa <- 1
  C_I <- 1
  C_r <- 1 / 4
  
  # Calculate empirical regret for each element in D_P_list
  emp_regret <- numeric(D_P_size)
  
  for (i in 1:D_P_size) {
    results <- transfer_learning_pricing_real(target_data, X_target, D_P_list[[i]], C_r, kappa, d, gamma, C_I, rf_model)
    emp_regret[i] <- sum(f_optimal - results$f_xt_pt)
    print(paste("Run i:", i))
  }
  
  # Calculate empirical regret for other algorithms
  results_ABE <- ABE_algorithm_real(target_data, X_target, d, rf_model)
  emp_regret_ABE <- sum(f_optimal - results_ABE$f_xt_pt)
  
  results_ExUCB <- ExUCB_algorithm_real(target_data, X_target, d, rf_model)
  emp_regret_ExUCB <- sum(f_optimal - results_ExUCB$f_xt_pt)
  
  # Return empirical regrets for all algorithms
  return(c(emp_regret, emp_regret_ABE, emp_regret_ExUCB))
}

# Simulation parameters
num_run <- 100
n_P_all <- nrow(source_data_ava_set)
n_P_list <- c(0, floor(0.25 * n_P_all), floor(0.5 * n_P_all), floor(0.75 * n_P_all), n_P_all)

# Initialize matrices for storing regrets
regret <- matrix(NA, nrow = length(n_P_list), ncol = num_run)
regret_ABE <- rep(NA, num_run)
regret_ExUCB <- rep(NA, num_run)

# Prepare D_P_list with subsets of source data
D_P_list <- lapply(n_P_list, function(n) {
  source_data <- tail(source_data_ava_set, n)
  X <- as.matrix(source_data[, variable_covariate])
  p <- as.vector(source_data$Price)
  Y <- as.vector(source_data$y)
  data.frame(X = X, p = p, Y = Y)
})

set.seed(123)
rf_model <- randomForest(formula, data = target_data_set, ntree = 500, importance = TRUE)

# Run simulations
for (j in 1:num_run){
result <- run_simulation(D_P_list, target_data_set, rf_model, f_optimal, C_r, kappa, d, gamma, 123 + 100 * j)
print(paste("Run j:", j))
regret[, j] <- result[1:length(n_P_list)]
regret_ABE[j] <- result[length(n_P_list) + 1]
regret_ExUCB[j] <- result[length(n_P_list) + 2]
}



