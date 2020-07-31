# Data Generation -------------------------------------------------------------
# Description: Generate simulated data
# Date: July, 2020
# -----------------------------------------------------------------------------
# Input:
# n    - Number of nodes
# k    - Dimension of latent position
# seed - Numeric seed
# -----------------------------------------------------------------------------
# Output: update!!!!!!!
# U_true           - True U
# alpha_list_true  - True alpha
# Lambda_list_true - True Lambda
# A_list_true      - One realization of adjacency matrices
# -----------------------------------------------------------------------------

generate_data = function(n, k, seed){
  # Setup for Z_true ----------------------------------------------------------
  library(pracma)
  set.seed(1)
  Z = matrix(rnorm(n * k, 0, 1), ncol = k)
  Z = scale(Z, center = TRUE, scale = FALSE)
  scale_fac = sqrt(n / norm(Z%*%t(Z), "F"))
  Z_true = Z * scale_fac

  # Setup for alpha_true ---------------------------------------------------
  set.seed(1)
  alpha = runif(n, 1, 3)
  alpha_true = - alpha / sum(alpha)
  
  # Setup for w_true, u_true, c_true ---------------------------------------
  w_true = rep(1/sqrt(k),k)
  u_true = 0
  
  # Setup for A_true -------------------------------------------------------
  generate_A = function(Z, alpha, w, u){
    n = dim(Z)[1]
    
    # generate edges
    Theta <- alpha %*% t(rep(1,n)) + rep(1, n) %*% t(alpha) + Z %*% t(Z)
    P <- sigmoid(Theta)
    A_abs <- matrix(0, nrow = n, ncol = n)
    P_ut <- P[upper.tri(P)]
    A_ut <- rbinom(length(P_ut), 1, P_ut)
    A_abs[upper.tri(A_abs)] <- A_ut
    A_abs <- A_abs + t(A_abs)
    
    # generate signs
    Eta <- (Z %*% w + u * rep(n,1)) %*% t(Z %*% w + u * rep(n,1))
    Q <- sigmoid(Eta)
    A_sign <- matrix(0, nrow = n, ncol = n)
    Q_ut <- Q[upper.tri(Q)]
    A_sign_ut <- rbinom(length(Q_ut), 1, Q_ut)
    A_sign[upper.tri(A_sign)] <- A_sign_ut
    A_sign <- A_sign + t(A_sign)
    A_sign <- 2 * A_sign - 1
    
    A = A_abs * A_sign
    
    return(A)
  }
  
  set.seed(seed)
  A_true = generate_A(Z_true, alpha_true, w_true, u_true)
  
  return(list(Z_true = Z_true,
              alpha_true = alpha_true,
              w_true = w_true,
              u_true = u_true,
              A_true = A_true))
}
