library(Matrix)
library(RGCCA)
library(Rlab)
library(RConics)
library(irlba)

# Function to solve quadratic equations
quad <- function(a, b, c) {
  a <- as.complex(a)
  answer <- c((-b + sqrt(b^2 - 4 * a * c)) / (2 * a),
              (-b - sqrt(b^2 - 4 * a * c)) / (2 * a))
  return(answer)
}

# Function to update matrix L
update_L <- function(P_k, S_k, Y_k, alpha, mu, m1, m2) {
  m <- min(m1, m2)
  X = P_k - S_k + (1/mu) * Y_k
  
  # Perform SVD on X
  svd_x <- svd(X)
  
  tau = alpha/mu
  S_d<- matrix(rep(NA, m*m), nrow=m, ncol=m)
  d<- diag(svd_x$d)
  
  for (i in 1: m){
    for (j in 1:m){
      S_d[i,j]<- sign(d[i,j])*max(abs(d[i,j])-tau, 0)
    }
  }
  L_new<- svd_x$u %*% S_d %*% t(svd_x$v)
  return(L_new)
}

# Function to update matrix S
update_S <- function(P_k, L_k, Y_k, beta, mu, m1, m2) {
  X = P_k - L_k + (1/mu) * Y_k
  tau = beta/mu
  
  # Apply soft thresholding
  S_new <- matrix(apply(X, c(1, 2), function(x) sign(x) * max(abs(x) - tau, 0)), nrow = m1, ncol = m2)
  return(S_new)
}

# Function to update matrix P
update_P <- function(L_k, S_k, Y_k, obs_matrix, m1, m2, mu) {
  P_new = matrix(NA, nrow = m1, ncol = m2)
  for (i in 1:m1) {
    for (j in 1:m2) {
      P_seq <- seq(0.00001, 1, by = 0.05)
      arg_func <- -mean(obs_matrix[i, j, ]) * log(P_seq) - (1 - mean(obs_matrix[i, j, ])) * log(1 - P_seq) + 
        (mu/2) * (P_seq - L_k[i, j] - S_k[i, j] + (1/mu) * Y_k[i, j])^2
      P_new[i, j] <- P_seq[which.min(arg_func)]
    }
  }
  return(P_new)
}

# Function to update matrix Mu
update_Mu <- function(L_k, S_k, Y_k, obs_matrix, m1, m2, mu) {
  Mu_new = matrix(NA, nrow = m1, ncol = m2)
  for (i in 1:m1) {
    for (j in 1:m2) {
      Mu_seq <- seq(0.00001, max(obs_matrix[i, j, ]) + 1, by = 0.1)
      arg_func <- ((mean(obs_matrix[i, j, ]) - Mu_seq) / 0.35)^2 + 
        (mu/2) * (Mu_seq - L_k[i, j] - S_k[i, j] + (1/mu) * Y_k[i, j])^2
      Mu_new[i, j] <- Mu_seq[which.min(arg_func)]
    }
  }
  return(Mu_new)
}

# Monte Carlo sampling for ERPCA
mc_samp <- function(n, P_all, S_all, Y_all, L_all, obs_test, alpha, beta, mu) {
  for (i in 1:(n - 1)) {
    # Update S matrix
    S_new <- update_S(P_all[,,i], L_all[,,i], Y_all[,,i], beta, mu, m1 = nrow(obs_test[,,1]), m2 = ncol(obs_test[,,1]))
    S_all[,,i + 1] <- S_new
    
    # Update L matrix
    L_new <- update_L(P_all[,,i], S_all[,,i + 1], Y_all[,,i], alpha, mu, m1 = nrow(obs_test[,,1]), m2 = ncol(obs_test[,,1]))
    L_all[,,i + 1] <- L_new
    
    # Update P matrix
    P_new <- update_P(L_all[,,i + 1], S_all[,,i + 1], Y_all[,,i], obs_test, m1 = nrow(obs_test[,,1]), m2 = ncol(obs_test[,,1]), mu)
    P_all[,,i + 1] <- P_new
    
    # Update Y matrix
    Y_new <- Y_all[,,i] + mu * (P_all[,,i + 1] - L_all[,,i + 1] - S_all[,,i + 1])
    Y_all[,,i + 1] <- Y_new
  }
  return(list(L_all = L_all, S_all = S_all, P_all = P_all, Y_all = Y_all))
}

# Main function for ERPCA
mc_function_mle <- function(obs_test, runs, mu_val, alpha_val, beta_val) {
  m1 <- nrow(obs_test[,,1])
  m2 <- ncol(obs_test[,,1])
  obs_test_mean <- apply(obs_test, 1:2, mean, na.rm = TRUE)
  svd_obs <- svd(obs_test_mean)
  
  # Initialization
  L <- svd_obs$u %*% diag(c(svd_obs$d[1:2], rep(0, (min(m1, m2) - 2)))) %*% t(svd_obs$v)
  Y <- matrix(0, m1, m2)
  P <- matrix(1, m1, m2)
  
  # Create empty arrays to store the results
  P_all <- array(NA, dim = c(m1, m2, runs))
  S_all <- array(NA, dim = c(m1, m2, runs))
  L_all <- array(NA, dim = c(m1, m2, runs))
  Y_all <- array(NA, dim = c(m1, m2, runs))
  
  # Initialize first iteration
  L_all[,,1] <- L
  P_all[,,1] <- P
  Y_all[,,1] <- Y
  
  # Run Monte Carlo sampling
  mc_result <- mc_samp(runs, P_all, S_all, Y_all, L_all, obs_test, alpha_val, beta_val, mu_val)
  
  return(mc_result)
}

# Monte Carlo sampling for RPCA
mc_samp_ls <- function(n, Mu_all, S_all, Y_all, L_all, obs_test, alpha, beta, mu) {
  for (i in 1:(n - 1)) {
    # Update S matrix
    S_new <- update_S(Mu_all[,,i], L_all[,,i], Y_all[,,i], beta, mu, m1 = nrow(obs_test[,,1]), m2 = ncol(obs_test[,,1]))
    S_all[,,i + 1] <- S_new
    
    # Update Mu matrix
    Mu_new <- update_Mu(L_all[,,i], S_all[,,i + 1], Y_all[,,i], obs_test, m1 = nrow(obs_test[,,1]), m2 = ncol(obs_test[,,1]), mu)
    Mu_all[,,i + 1] <- Mu_new
    
    # Update Y matrix
    Y_new <- Y_all[,,i] + mu * (Mu_all[,,i + 1] - L_all[,,i] - S_all[,,i + 1])
    Y_all[,,i + 1] <- Y_new
    
    # Update L matrix
    L_new <- update_L(Mu_all[,,i + 1], S_all[,,i + 1], Y_all[,,i + 1], alpha, mu, m1 = nrow(obs_test[,,1]), m2 = ncol(obs_test[,,1]))
    L_all[,,i + 1] <- L_new
    
    # Convergence check
    if (norm(Mu_all[,,i + 1] - L_all[,,i + 1] - S_all[,,i + 1], "F") <= 1e-5 * norm(Mu_all[,,i + 1], "F")) {
      break
    }
  }
  return(list(L_all = L_all, S_all = S_all, Mu_all = Mu_all, Y_all = Y_all, number = i))
}

# Main function for RPCA
mc_function_ls <- function(obs_test, runs, mu, alpha, beta) {
  m1 <- nrow(obs_test[,,1])
  m2 <- ncol(obs_test[,,1])
  obs_test_mean <- apply(obs_test, 1:2, mean, na.rm = TRUE)
  svd_obs <- svd(obs_test_mean)
  
  # Initialization
  L <- svd_obs$u %*% diag(c(svd_obs$d[1:2], rep(0, (min(m1, m2) - 2)))) %*% t(svd_obs$v)
  Y <- matrix(0, m1, m2)
  Mu <- matrix(1, m1, m2)
  
  # Create empty arrays to store the results
  Mu_all <- array(NA, dim = c(m1, m2, runs))
  S_all <- array(NA, dim = c(m1, m2, runs))
  L_all <- array(NA, dim = c(m1, m2, runs))
  Y_all <- array(NA, dim = c(m1, m2, runs))
  
  # Initialize first iteration
  L_all[,,1] <- L
  Mu_all[,,1] <- Mu
  Y_all[,,1] <- Y
  
  # Run Monte Carlo sampling
  mc_result <- mc_samp_ls(runs, Mu_all, S_all, Y_all, L_all, obs_test, alpha, beta, mu)
  
  return(mc_result)
}
