fit_RecSBART_martingale <- function(X_train,
                                    X_test,
                                    recurrent_train,
                                    recurrent_test,
                                    terminal_train,
                                    terminal_test,
                                    num_burn = 2000,
                                    num_save = 500,
                                    num_thin = 1) {
  
  # --- Basic setup ---
  n <- nrow(X_train)
  p <- ncol(X_train)
  num_recurrent <- sapply(recurrent_train, length)
  N <- sum(num_recurrent)
  
  # --- Scale time to [0, 1] ---
  recurrent_train <- lapply(recurrent_train, function(i) i / max(terminal_train))
  terminal_train <- terminal_train / max(terminal_train)
  
  # --- Empirical intensity estimate ---
  safe_max <- function(x) if (length(x) == 0) 0 else max(x, na.rm = TRUE)
  lambda <- mean(num_recurrent) / mean(sapply(recurrent_train, safe_max))
  
  shape_lambda_prior <- 1
  rate_lambda_prior <- 1 / lambda
  
  # ---hypers for trees---
  hypers <- Hypers(cbind(sapply(recurrent_train, safe_max), X_train),
                   num_recurrent / terminal_train, sigma_hat = 1, num_tree = 50)
  opts <- Opts(update_sigma = TRUE, update_s = FALSE,
               update_alpha = FALSE, update_sigma_mu = FALSE)
  my_forest <- MakeForest(hypers, opts, warn = FALSE)
  
  # --- Frailty prior ---
  shape_W_prior <- 40
  W <- rgamma(n, shape = shape_W_prior, rate = shape_W_prior)
  
  loglik_shape_W_prior <- function(shape_W_prior) {
    sum(dgamma(W, shape = shape_W_prior, rate = shape_W_prior, log = TRUE))
  }
  
  num_iter <- num_burn + num_save * num_thin
  
  # --- Running mean storage (memory-efficient) ---
  martingale_mean <- vector("list", n)
  for (i in 1:n) {
    martingale_mean[[i]] <- numeric(num_recurrent[i] + 1)
  }
  
  sort_time <- sort(unique(unlist(recurrent_train)))
  partition <- c(0, sort_time)
  delta <- diff(partition)
  
  start.time <- Sys.time()
  
  for (iter in 1:num_iter) {
    z_list <- vector("list", n)
    X_list <- vector("list", n)
    time_list <- vector("list", n)
    q <- integer(n)
    m <- integer(n)
    
    # --- Data augmentation ---
    for (i in seq_len(n)) {
      q[i] <- rpois(1, 2 * lambda * W[i] * terminal_train[i])
      if (q[i] > 0) {
        c_vals <- runif(q[i], 0, 2 * lambda * W[i] * terminal_train[i])
        a <- c_vals / (2 * lambda * W[i])
        a_x_a <- cbind(a, matrix(rep(X_train[i, ], length(a)), ncol = p, byrow = TRUE))
        u <- runif(q[i])
        l <- my_forest$do_predict(a_x_a)
        g_vals <- a[which(u < 1 - pnorm(l))]
        m[i] <- length(g_vals)
      } else {
        m[i] <- 0
        g_vals <- numeric(0)
      }
      
      # Build augmented data
      if (num_recurrent[i] == 0 && m[i] == 0) {
        z_list[[i]] <- numeric(0)
        X_list[[i]] <- matrix(nrow = 0, ncol = p)
        time_list[[i]] <- numeric(0)
        next
      }
      
      # Recurrent and/or augmented draws
      z_t <- z_g <- numeric(0)
      if (num_recurrent[i] > 0) {
        z_t <- sapply(recurrent_train[[i]], function(t_ij) {
          mu <- my_forest$do_predict(cbind(t_ij, t(X_train[i, ])))
          msm::rtnorm(1, mean = mu, sd = 1, lower = 0, upper = Inf)
        })
      }
      if (m[i] > 0) {
        z_g <- sapply(g_vals, function(gj) {
          mu <- my_forest$do_predict(cbind(gj, t(X_train[i, ])))
          msm::rtnorm(1, mean = mu, sd = 1, lower = -Inf, upper = 0)
        })
      }
      z_list[[i]] <- c(z_t, z_g)
      X_list[[i]] <- matrix(rep(X_train[i, ], num_recurrent[i] + m[i]),
                            ncol = p, byrow = TRUE)
      time_list[[i]] <- c(recurrent_train[[i]], g_vals)
    }
    
    # --- Combine once per iteration ---
    z_all <- unlist(z_list, use.names = FALSE)
    X_store <- do.call(rbind, X_list)
    time_points <- unlist(time_list, use.names = FALSE)
    
    # --- Parameter updates ---
    sum_points <- N + sum(m)
    lambda <- rgamma(1, shape = shape_lambda_prior + sum_points,
                     rate = rate_lambda_prior + 2 * sum(W * terminal_train))
    
    for (i in 1:n) {
      W[i] <- rgamma(1, shape = shape_W_prior + m[i] + num_recurrent[i],
                     rate = shape_W_prior + lambda * terminal_train[i])
    }
    
    mu_hat <- my_forest$do_gibbs(cbind(time_points, X_store), z_all,
                                 cbind(time_points, X_store), 1)
    
    shape_W_prior <- (diversitree::mcmc(
      lik = loglik_shape_W_prior, nsteps = 1, w = 1, x.init = c(shape_W_prior),
      prior = function(shape_W_prior)
        dgamma(shape_W_prior, shape = 70, rate = 0.5, log = TRUE),
      lower = 20, upper = Inf))$pars
    
    # --- Save (update running mean only) ---
    if (iter > num_burn && ((iter - num_burn - 1) %% num_thin == 0)) {
      save_idx <- (iter - num_burn) / num_thin
      frac <- 1 / save_idx  # online mean update weight
      
      for (i in 1:n) {
        x_i <- X_train[i, ]
        W_i <- W[i]
        n_rec <- num_recurrent[i]
        
        time_full <- unique(c(sort_time, terminal_train[i]))
        XT_full <- matrix(rep(x_i, length(time_full)), ncol = p, byrow = TRUE)
        phi_full <- pnorm(as.numeric(my_forest$do_predict(cbind(time_full, XT_full))))
        diff_full <- diff(c(0, time_full))
        integral_full <- sum(diff_full * phi_full)
        cumulative_full <- integral_full * 2 * lambda * W_i
        
        # terminal component
        old_val <- martingale_mean[[i]][n_rec + 1]
        martingale_mean[[i]][n_rec + 1] <- old_val + frac * (cumulative_full - old_val)
        
        if (n_rec == 0) next
        
        for (j in 1:n_rec) {
          t_ij <- recurrent_train[[i]][j]
          k <- sum(sort_time <= t_ij)
          integral <- 0
          if (k > 0) {
            XT_k <- matrix(rep(x_i, k), ncol = p, byrow = TRUE)
            phi_k <- pnorm(as.numeric(my_forest$do_predict(cbind(sort_time[1:k], XT_k))))
            integral <- sum(delta[1:k] * phi_k)
          }
          cumulative <- integral * 2 * lambda * W_i
          
          old_val <- martingale_mean[[i]][j]
          martingale_mean[[i]][j] <- old_val + frac * (cumulative - old_val)
        }
      }
    }
    
    if (iter %% 50 == 0) gc()  # encourage memory release
  }
  
  stop.time <- Sys.time()
  cat("Runtime: ", round(difftime(stop.time, start.time, units = "mins"), 2), " minutes\n")
  
  return(list(martingale_mean = martingale_mean))
}
