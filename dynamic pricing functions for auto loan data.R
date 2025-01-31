

transfer_learning_pricing_real <- function(target_data, X_target, D_P, C_r, kappa_tilde_r , d, gamma, C_I, rf_model){
  # Initialize
  t <- 1
  n_P_total <- ifelse(is.null(D_P), 0, nrow(D_P))
  
  n_Q <- nrow(target_data)
  
  without_log_term <- max(n_Q, (kappa_tilde_r * n_P_total)^((d + 3) / (d + 3 + gamma)))
  
  log_term <- log(without_log_term)
  
  tilde_r <- C_r *  (log_term/without_log_term)^(1/(d+3))
  
  # Initialize ball ID counter
  ball_id_counter <- 1
  
  # Initialize the initial ball B0
  B0 <- list()
  B0$id <- ball_id_counter
  ball_id_counter <- ball_id_counter + 1
  B0$center <- rep(0.5, d + 1)  # Assuming X and p are in [0,1]
  B0$radius <- 0.5  # Maximum possible radius in [0,1]^{d+1} for ℓ_∞-ball
  n_re_P <- compute_n_re_P(B0, D_P)
  B0$n_t <- n_re_P$n_P
  B0$re_t <- n_re_P$re_P
  B0$n_P <- n_re_P$n_P
  B0$re_P <- n_re_P$re_P
  B0$omega <- compute_omega(B0, log_term)
  B0$n_BQ <- compute_n_BQ(B0)
  
  # Initialize Balls_list and A_t
  Balls_list <- list(B0)
  A_t <- Balls_list
  
  results <- data.frame(
    t = integer(n_Q),
    p_t = numeric(n_Q),
    f_xt_pt = numeric(n_Q),
    Y_t = integer(n_Q)
  )
  
  # Main loop
  while (t <= n_Q) {
    # Observe X_t
    X_t <- X_target[t, ]  # Assuming X_t is drawn uniformly from [0,1]^d
    # Step 1: Find relevant_t
    relevant_balls <- find_relevant_balls(X_t, A_t, d)
    # If no relevant balls found, select the smallest ball containing X_t
    if (length(relevant_balls) == 0) {
      warning(paste("No relevant balls found at time", t, ". Skipping step."))
      t <- t + 1
      next
    }
    
    # Step 2: Compute I_t^{pre}(B) for all B in A_t
    for (i in seq_along(A_t))  {
      B <- A_t[[i]]
      B$v_t <- ifelse(B$n_t > 0, B$re_t / B$n_t, 0)
      B$conf_t <- ifelse(B$n_t > 0, 2 * sqrt(log_term / B$n_t), Inf)
      B$I_t_pre <- B$v_t + C_I * B$radius + B$conf_t
      # Update B in A_t using its id
      A_t[[i]] <- B
    }
    
    # Step 3: Compute I_t(B) for B in relevant_balls
    I_t_values <- numeric(length(relevant_balls))
    for (i in seq_along(relevant_balls)) {
      B <- relevant_balls[[i]]
      idx <- which(sapply(A_t, function(x) x$id == B$id))
      B <- A_t[[idx]]
      min_value <- Inf
      for (B_prime in A_t) {
        dist_centers <- max(abs(B$center - B_prime$center))
        value <- B_prime$I_t_pre + C_I * dist_centers
        if (value < min_value) {
          min_value <- value
        }
      }
      B$I_t <- C_I * B$radius + min_value
      # Update B in relevant_balls and A_t
      relevant_balls[[i]] <- B
      idx <- which(sapply(A_t, function(x) x$id == B$id))
      A_t[[idx]] <- B
      I_t_values[i] <- B$I_t
    }
    
    # Step 4: Select B^{sel}
    index_sel <- which.max(I_t_values)
    B_sel <- relevant_balls[[index_sel]]
    
    
    # Step 5: Select p_t
    # Select any p such that (X_t, p_t) ∈ dom(B_sel, A_t)ct p_t as the center's p-coordinate
    p_t <- compute_p_t(X_t, B_sel, A_t, d)
    
    target_data$Price[t]<- p_t
    
    f_xt_pt <- predict(rf_model, newdata = target_data[t, ])
    
    prob_p <-  f_xt_pt  / p_t
    prob_p[!is.finite(prob_p)] <- 0
    prob_p <- pmax(pmin(prob_p, 1), 0)
    
    Y_t <- ifelse(runif(1) < prob_p, p_t, 0)
    
    # Record results
    results$t[t] <- t
    results$p_t[t] <- p_t
    results$f_xt_pt[t] <- f_xt_pt
    results$Y_t[t] <- Y_t
    
    
    
    # Step 6: While loop to split balls
    # Use current n_t(B_sel) and re_t(B_sel) before updating with new observation
    while (B_sel$n_t >= B_sel$n_BQ + B_sel$n_P && B_sel$radius >= 2 * tilde_r) {
      # Create new ball B'
      B_new <- list()
      B_new$id <- ball_id_counter
      ball_id_counter <- ball_id_counter + 1
      B_new$center <- c(X_t, p_t)
      B_new$radius <- B_sel$radius / 2
      n_re_P <- compute_n_re_P(B_new, D_P)
      B_new$n_t <- n_re_P$n_P  # Initialize with source data
      B_new$re_t <- n_re_P$re_P
      B_new$n_P <- n_re_P$n_P
      B_new$re_P <- n_re_P$re_P
      B_new$omega <- compute_omega(B_new, log_term)
      B_new$n_BQ <- compute_n_BQ(B_new)
      B_new$v_t <- ifelse(B_new$n_t > 0, B_new$re_t / B_new$n_t, 0)
      B_new$conf_t <- ifelse(B_new$n_t > 0, 2 * sqrt(log_term / B_new$n_t), Inf)
      B_new$I_t_pre <- B_new$v_t + C_I * B_new$radius + B_new$conf_t
      # Add B_new to Balls_list and A_t
      Balls_list <- c(Balls_list, list(B_new))
      A_t <- c(A_t, list(B_new))
      # Update B_sel to B_new
      B_sel <- B_new
    }
    
    # Step 7: For all B in A_t except B_sel, update n_{t+1}(B) and re_{t+1}(B)
    for (i in seq_along(A_t)) {
      B <- A_t[[i]]
      if (B$id != B_sel$id) {
        # Other balls remain the same
        B$n_t <- B$n_t
        B$re_t <- B$re_t
        A_t[[i]] <- B
      }
    }
    
    # Step 8: Update n_{t+1}(B^{sel}) and re_{t+1}(B^{sel}), then increment t
    idx_sel <- which(sapply(A_t, function(x) x$id == B_sel$id))
    A_t[[idx_sel]]$n_t <- A_t[[idx_sel]]$n_t + 1
    A_t[[idx_sel]]$re_t <- A_t[[idx_sel]]$re_t + Y_t
    B_sel <- A_t[[idx_sel]]
    
    
    
    # Increment time
    t <- t + 1
  }
  # Return results
  return(results)
  
}



ABE_algorithm_real<- function(target_data, X_target, d,rf_model, M = 0.1) {
  # T: Total time periods
  # d: Dimension of covariates
  # revenue_function: Function to compute revenue given X_t and p_t
  # M_1 to M_4: Constants specific to the problem context
  
  T<-nrow(X_target)
  # Parameters
  K <- floor(log(T) / ((d + 4) * log(2)))
  N_k <- ceiling(log(T))
  
  # Functions to compute Delta_k and n_k
  Delta_k <- function(k) {
    2^(-k) * log(T)
  }
  
  n_k <- function(k) {
    max(0, ceiling(((2^(4 * k + 15)) / (M^2 * log(T)^3)) * 
                     (log(T) + log(log(T)) - (d + 2) * k * log(2))))
  }
  
  # Initialize root bin
  B_root <- list(
    a = rep(0, d),   # Lower bounds of the bin in each dimension
    b = rep(1, d),   # Upper bounds of the bin in each dimension
    level = 0,
    N = 0,           # Number of customers observed in the bin
    N_total = 0,     # Total number of times the bin was selected
    p_l = 0,         # Lower bound of price interval
    p_u = 1,         # Upper bound of price interval
    delta = 1 / (N_k - 1),   # Grid size for prices
    Y_j = rep(0, N_k),       # Average revenue for each price
    N_j = rep(0, N_k),       # Number of times each price was used
    price_set = seq(0, 1, length.out = N_k),  # Set of prices to explore
    p_star = NULL,   # Empirically-optimal price
    children = list()  # List to store child bins after splitting
  )
  
  
  # Initialize partition with the root bin
  P_t <- list(B_root)
  
  # Function to find the index of the bin containing X_t
  find_bin_index <- function(P_t, X_t) {
    for (i in seq_along(P_t)) {
      B <- P_t[[i]]
      if (all(X_t >= B$a) && all(X_t <= B$b)) {
        return(i)
      }
    }
    return(NULL)
  }
  
  # Function to split a bin into child bins
  split_bin <- function(B) {
    d <- length(B$a)
    level <- B$level + 1
    k <- level
    Delta_k_val <- Delta_k(k)
    N_k_val <- N_k
    
    # Generate all combinations for splitting (2^d child bins)
    comb <- expand.grid(rep(list(0:1), d))
    children <- list()
    
    for (i in 1:nrow(comb)) {
      i_vec <- as.integer(comb[i, ])
      a_child <- B$a
      b_child <- B$b
      for (j in 1:d) {
        mid_point <- (B$a[j] + B$b[j]) / 2
        if (i_vec[j] == 0) {
          b_child[j] <- mid_point
        } else {
          a_child[j] <- mid_point
        }
      }
      # Set price interval for child bins based on the parent's empirically-optimal price
      p_star <- B$p_star
      p_l <- max(p_star - Delta_k_val / 2, 0)
      p_u <- min(p_star + Delta_k_val / 2, 1)
      delta <- (p_u - p_l) / (N_k_val - 1)
      price_set <- seq(p_l, p_u, length.out = N_k_val)
      
      child <- list(
        a = a_child,
        b = b_child,
        level = level,
        N = 0,
        N_total = 0,
        p_l = p_l,
        p_u = p_u,
        delta = delta,
        Y_j = rep(0, N_k_val),
        N_j = rep(0, N_k_val),
        price_set = price_set,
        p_star = NULL,
        children = list()
      )
      children[[length(children) + 1]] <- child
    }
    return(children)
  }
  
  results <- data.frame(
    t = integer(T),
    p_t = numeric(T),
    f_xt_pt = numeric(T),
    Y_t = integer(T)
  )
  
  # Main loop over time periods
  for (t in 1:T) {
    # Observe covariate X_t (assuming uniform distribution for demonstration)
    X_t <-X_target[t, ]
    
    # Find the bin containing X_t
    B_index <- find_bin_index(P_t, X_t)
    if (is.null(B_index)) {
      stop(paste("No bin found for X_t at t =", t))
    }
    B <- P_t[[B_index]]
    
    # Update bin statistics
    B$N <- B$N + 1
    k <- B$level
    
    if (k < K) {
      if (B$N < n_k(k)) {
        # Not enough data observed in B, continue exploring prices
        j <- (B$N - 1) %% N_k
        p_t <- B$p_l + j * B$delta
        
        # Observe revenue using the provided revenue function
        target_data$Price[t]<- p_t
        
        f_xt_pt <- predict(rf_model, newdata = target_data[t, ])
        
        prob_p <-  f_xt_pt  / p_t
        prob_p[!is.finite(prob_p)] <- 0
        prob_p <- pmax(pmin(prob_p, 1), 0)
        prob_p <- pmax(pmin(prob_p, 1), 0)
        
       Y_t <- ifelse(runif(1) < prob_p, p_t, 0)
        
        
        
        # Update counts and average revenues
        B$N_total <- B$N_total + 1
        B$Y_j[j + 1] <- ( B$N_j[j + 1]*B$Y_j[j + 1] + Y_t ) / (B$N_j[j + 1]+1)
        B$N_j[j + 1] <- B$N_j[j + 1] + 1
        
        # Update the bin in the partition
        P_t[[B_index]] <- B
      } else {
        # Sufficient data observed, find empirically-optimal price
        p_star_index <- which.max(B$Y_j)
        p_star <- B$price_set[p_star_index]
        B$p_star <- p_star
        
        # Split the bin into child bins and initialize them
        children <- split_bin(B)
        
        # Update the partition: remove B and add its children
        P_t[[B_index]] <- NULL
        P_t <- c(P_t, children)
      }
    } else {
      # Reached maximal level, use mid-point price
      p_t <- (B$p_l + B$p_u) / 2
      
      target_data$Price[t]<- p_t
      
      f_xt_pt <- predict(rf_model, newdata = target_data[t, ])
      
      prob_p <-  f_xt_pt  / p_t
      prob_p[!is.finite(prob_p)] <- 0
      prob_p <- pmax(pmin(prob_p, 1), 0)
      prob_p <- pmax(pmin(prob_p, 1), 0)
      
      Y_t <- ifelse(runif(1) < prob_p, p_t, 0)
      
      # Update bin statistics
      B$N_total <- B$N_total + 1
      
      # Update the bin in the partition
      P_t[[B_index]] <- B
    }
    
    results$t[t] <- t
    results$p_t[t] <- p_t
    results$f_xt_pt[t] <- f_xt_pt
    results$Y_t[t] <- Y_t
  }
  
  # Return the final partition and any other relevant information
  return(results)
}



# Updated ExUCB Algorithm Function with User's Step
ExUCB_algorithm_real <- function(target_data,  X_target, d, rf_model) {
  # Parameters for ExUCB
  n_Q<-nrow(X_target)
  
  beta <- 2/3   # Exploration phase exponent
  gamma <- 1/6  # UCB phase exponent
  C1 <- 1       # Exploration phase constant
  C2 <- 20      # Discretization constant
  lam <- 0.1    # Regularization parameter in UCB
  CU <- 1/40    # UCB constant
  pmax <- 1     # Maximum price
  B <- 1        # Maximum possible revenue
  
  
  
  # Fix the number of episodes to 10
  epi <- 10
  
  # Calculate the initial episode length based on n_Q and epi
  ini <- ceiling(n_Q / (2^epi - 1))  # 2^10 -1 = 1023
  T_total <- ini * (2^epi - 1)
  
  # Adjust T_total and ini to match n_Q exactly
  if (T_total > n_Q) {
    # The last episode might need to be shorter
    # We'll handle this in the loop by ensuring we don't exceed n_Q
    T <- n_Q
  } else {
    T <- T_total
  }
  
  # Initialize variables
  S <- numeric(T)      # Prices selected
  er <- numeric(T)     # Expected revenues (based on algorithm's selection)
  y <- numeric(T)      # Observed revenues
  r <- numeric(T)      # Revenue from current selection
  x1 <- X_target[1:T, ]  # Use the first T samples from X_target
  thx <- numeric(T)    # Estimated theta x
  
  
  
  
  cur_end <- 0
  for (j in 1:epi){
    subT <- ini * 2^(j - 1)
    # Adjust subT if it exceeds the remaining time
    if (cur_end + subT > T) {
      subT <- T - cur_end
      if (subT <= 0) break
    }
    
    # --------------------
    # Exploration Phase
    # --------------------
    explore_l <- ceiling(C1 * subT^(beta))
    min_explore <- d + 1
    if (explore_l < min_explore) {
      explore_l <- min_explore
    }
    
    if (explore_l >= subT) {
      explore_l <- subT
    }
    
    for (i in 1:explore_l){
      idx <- i + cur_end
      S[idx] <- runif(1, 0, B)
      
      # Compute optimal price and expected revenue
      target_data$Price[idx]<- S[idx]
      er[idx] <-predict(rf_model, newdata = target_data[idx, ])
      
      prob_p <-  er[idx]  / S[idx]
      prob_p[!is.finite(prob_p)] <- 0
      prob_p <- pmax(pmin(prob_p, 1), 0)
      prob_p <- pmax(pmin(prob_p, 1), 0)
      
      
      # Simulate observed revenue
      y[idx] <-ifelse(runif(1) < prob_p, S[idx], 0)
      
      r[idx] <- y[idx]/ S[idx]
      
      
    }
    
    # --------------------
    # Parameter Estimation
    # --------------------
    # Prepare data for regression: assuming linear relationship
    dat <- data.frame(y = r[(cur_end + 1):(cur_end + explore_l)], x1[(cur_end + 1):(cur_end + explore_l), ])
    zero_var_cols <- sapply(dat[, -1], function(x) var(x) == 0)
    
    glmfit <- glm(y ~ ., data = dat)
    thetahat2 <- coef(glmfit)
    thetahat <- thetahat2[-1]  
    thetahat[is.na(thetahat)] <- 0 
    # Store absolute estimation error for each component
    
    # --------------------
    # UCB Phase
    # --------------------
    exploit_l <- subT - explore_l
    if (exploit_l <= 0) {
      cur_end <- cur_end + subT
      next
    }
    intv <- ceiling(C2 * exploit_l^(gamma))
    price_set <- seq(0, pmax, length.out = intv)
    me0 = rep(0,intv)
    ti0 = rep(0,intv)
    u1 = pmax
    u2 = sum(abs(thetahat))
    u = u1 + 2*u2
    ku = u/intv
    
    for (i in (explore_l + 1):subT){
      idx <- i + cur_end
      
      # Compute expected theta x
      thx[idx] <- sum(thetahat * x1[idx, ])
      # Compute beta_t
      beta_t <- CU * max(1, ((lam * intv)^(1/2) / pmax + sqrt(2 * log(exploit_l) + intv * log((lam * intv + (i - 1) * pmax^(2)) / (lam * intv))))^(2))
      
      # Compute dex1 and dex2 based on thx[idx]
      dex1 <- floor((-thx[idx] + u2 + ku / 2) / ku) + 1
      dex2 <- floor((u1 - thx[idx] + u2 + ku / 2) / ku)
      num <- dex2 - dex1 + 1
      if (num <= 0){
        next  # Skip if no valid arms
      }
      rma <- (2 * dex1 - 1) * ku / 2 - u2 + thx[idx]
      
      if (i == explore_l + 1){
        bc <- sample(1:num, 1)
      } else {
        me <- me0[dex1:dex2]
        ti <- ti0[dex1:dex2]
        if (any(ti == 0)){
          zero_indices <- which(ti == 0)
          bc <- zero_indices[sample(length(zero_indices), 1)]
        } else {
          inde <- numeric(num)
          for (i1 in 1:num){
            inde[i1] <- ((i1 - 1) * ku + rma) * (me[i1] + sqrt(beta_t / (lam + ti[i1])))
          }
          bc_candidates <- which(inde == max(inde))
          bc <- bc_candidates[sample(length(bc_candidates), 1)]
        }
      }
      S[idx] <- (bc - 1) * ku + rma
      
      target_data$Price[idx]<- S[idx]
      er[idx] <-predict(rf_model, newdata = target_data[idx, ])
      
      prob_p <-  er[idx]  / S[idx]
      prob_p[!is.finite(prob_p)] <- 0
      prob_p <- pmax(pmin(prob_p, 1), 0)
      prob_p <- pmax(pmin(prob_p, 1), 0)
      
      
      # Simulate observed revenue
      y[idx] <-ifelse(runif(1) < prob_p, S[idx], 0)
      
      r[idx] <- y[idx]/ S[idx]

      
      me0[dex1-1+bc] = (me0[dex1-1+bc]*(lam+ti0[dex1-1+bc])+S[i+cur_end]*r[i+cur_end])/(lam+ti0[dex1-1+bc]+(S[i+cur_end])^(2))
      ti0[dex1-1+bc] = ti0[dex1-1+bc] + (S[i+cur_end])^(2)
    }
    cur_end <- cur_end + subT
    if (cur_end >= T) break
  }
  
  # --------------------
  # Compute Cumulative Regret
  # --------------------
  # Calculate Empirical Regret
  # emp_regret <- sum(ber[1:T] - er[1:T])
  # emp_regret_pre <- sum(Y_optimal[1:T] - y[1:T])
  
  # However, Y_optimal is computed outside this function
  # We'll compute it in the run_simulation function
  
  # Prepare results similar to other algorithms
  results <- data.frame(
    t = 1:T,
    p_t = S[1:T],
    f_xt_pt = er[1:T],
    Y_t = y[1:T]
  )
  
  return(results)
}



