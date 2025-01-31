
# Function to compute omega(B)
compute_omega <- function(B, log_term) {
  omega <- ceiling(log_term / (B$radius^2))
  return(omega)
}

# Function to compute T_B^Q
compute_n_BQ <- function(B) {
  if (B$omega < B$n_P) {
    return(0)
  } else {
    return(B$omega)
  }
}

# Function to check if a ball B contains X_t
contains_X_t <- function(B, X_t, d) {
  return(all(abs(X_t - B$center[1:d]) <= B$radius))
}

# Function to compute feasible p intervals for a ball B, excluding overlaps from smaller balls
compute_feasible_p_intervals <- function(B, A_t, X_t, d) {
  # Define initial feasible p range for B
  p_min <- max(0, B$center[d + 1] - B$radius)
  p_max <- min(1, B$center[d + 1] + B$radius)
  feasible_intervals <- list(c(p_min, p_max))
  
  # Iterate over all smaller balls to exclude their p ranges
  for (B_prime in A_t) {
    # Skip if B_prime is not smaller or does not contain X_t
    if (B_prime$radius >= B$radius) next
    if (!contains_X_t(B_prime, X_t, d)) next
    
    # Define p range for B_prime
    p_prime_min <- max(0, B_prime$center[d + 1] - B_prime$radius)
    p_prime_max <- min(1, B_prime$center[d + 1] + B_prime$radius)
    
    # Exclude [p_prime_min, p_prime_max] from feasible_intervals
    new_feasible_intervals <- list()
    for (interval in feasible_intervals) {
      a <- interval[1]
      b <- interval[2]
      
      # No overlap
      if (b <= p_prime_min || a >= p_prime_max) {
        new_feasible_intervals <- c(new_feasible_intervals, list(c(a, b)))
      } else {
        # Overlap exists, split interval if necessary
        if (a < p_prime_min) {
          new_feasible_intervals <- c(new_feasible_intervals, list(c(a, p_prime_min)))
        }
        if (b > p_prime_max) {
          new_feasible_intervals <- c(new_feasible_intervals, list(c(p_prime_max, b)))
        }
      }
    }
    feasible_intervals <- new_feasible_intervals
  }
  
  return(feasible_intervals)
}

# Function to identify relevant balls at time t
find_relevant_balls <- function(X_t, A_t, d) {
  relevant_balls <- list()
  
  # Filter balls that contain X_t
  containing_balls <- Filter(function(B) contains_X_t(B, X_t, d), A_t)
  
  # Sort containing_balls by decreasing radius to handle larger balls first
  containing_balls <- containing_balls[order(sapply(containing_balls, function(B) -B$radius))]
  
  for (B in containing_balls) {
    # Compute feasible p intervals for B
    feasible_p <- compute_feasible_p_intervals(B, A_t, X_t, d)
    
    # Check if there's at least one feasible p interval
    if (length(feasible_p) > 0 && any(sapply(feasible_p, function(interval) interval[2] > interval[1]))) {
      relevant_balls <- c(relevant_balls, list(B))
    }
  }
  
  return(relevant_balls)
}


# Function to compute feasible p_t using compute_feasible_p_intervals
compute_p_t <- function(X_t, B_sel, A_t, d) {
  # Compute feasible p intervals for B_sel
  feasible_intervals <- compute_feasible_p_intervals(B_sel, A_t, X_t, d)
  
  # If no feasible intervals remain, select p_t as the center price of B_sel
  if (length(feasible_intervals) == 0) {
    p_t <- B_sel$center[d + 1]
    warning("No feasible p_t intervals after exclusion. Using center price of B_sel.")
    return(p_t)
  }
  
  # Filter intervals with positive length
  feasible_intervals <- feasible_intervals[sapply(feasible_intervals, function(interval) interval[2] > interval[1])]
  
  # If still no feasible intervals, select p_t as the center price of B_sel
  if (length(feasible_intervals) == 0) {
    p_t <- B_sel$center[d + 1]
    warning("No feasible p_t intervals with positive length. Using center price of B_sel.")
    return(p_t)
  }
  
  # Calculate the total length of feasible intervals
  lengths <- sapply(feasible_intervals, function(x) x[2] - x[1])
  total_length <- sum(lengths)
  
  # Generate a random point in [0, total_length]
  u <- runif(1, 0, total_length)
  
  # Determine which interval u falls into
  cum_lengths <- cumsum(lengths)
  interval_index <- which(u <= cum_lengths)[1]
  
  # Position within the selected interval
  if (interval_index == 1) {
    pos_in_interval <- u
  } else {
    pos_in_interval <- u - cum_lengths[interval_index - 1]
  }
  
  # Sample p_t accordingly
  p_t <- feasible_intervals[[interval_index]][1] + pos_in_interval
  
  return(p_t)
}





