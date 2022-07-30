# -----------------------------------------------------------------------------
# Description: Projected Gradient Descent Algorithm for Separate Estimation
# -----------------------------------------------------------------------------
# Input:
# A_true    - Observed signed adjacency matrix
# M         - Missing index matrix
# Z_0       - Initialization of Z
# alpha_0   - Initialization of alpha
# v_0       - Initialization of v
# tau_z     - Learning rate for Z
# tau_alpha - Learning rate for alpha
# tau_v     - Learning rate for v
# num_iter  - Maximum iteration number
# rel_tol   - Tolerance for algorithm convergence
# -----------------------------------------------------------------------------
# Output: 
# Z_0       - Estimate for Z from separate estimation
# alpha_0   - Estimate for alpha from separate estimation
# v_0       - Estimate for v from separate estimation
# nll       - Negative log likelihood for all iterates
# -----------------------------------------------------------------------------

Z_estimation = function(A_true, M, Z_0, alpha_0, tau_z, tau_alpha, num_iter, rel_tol){
  n = dim(Z_0)[1]
  k = dim(Z_0)[2]
  ones = rep(1, n)
  A_abs = abs(A_true)
  
  Theta = alpha_0 %*% t(ones) + ones %*% t(alpha_0) + Z_0 %*% t(Z_0)
  
  stopflag = FALSE
  iter = 0
  
  nll_collect = rep(0, num_iter) # negative log-likelihood
  tsr = (A_abs * Theta + log(1 - sigmoid(Theta))) * M
  nll_collect[1] = -sum(tsr[upper.tri(tsr)]) / n^2
  
  while(!stopflag){
    # (1) Update Z
    Z_1 = Z_0 + 2 * tau_z * ((M * (A_abs - sigmoid(Theta))) %*% Z_0)  # Gradient descent step
    # (2) Update alpha
    alpha_1 = alpha_0 + 2 * tau_alpha * ((M * (A_abs - sigmoid(Theta))) %*% ones)
    # (3) Projection
    Z_1 = scale(Z_1, center = TRUE, scale = FALSE)  # Centering Z
    
    # (4) Prepare for the next iteration    
    grad_norm_z = (norm(Z_1 - Z_0, "F"))^2
    grad_norm_alpha = (norm(alpha_1 - alpha_0, "F"))^2
    Z_0 <- Z_1 
    alpha_0 <- alpha_1
    Theta = alpha_0 %*% t(ones) + ones %*% t(alpha_0) + Z_0 %*%  t(Z_0)
    
    # (5) Compute Negative Log-likelihood
    tsr = (A_abs * Theta + log(1 - sigmoid(Theta))) * M 
    nll_1 = -sum(tsr[upper.tri(tsr)]) / n^2
    # print(nll_1)
    rel_nll = abs(nll_collect[iter + 1] - nll_1) / abs(nll_collect[iter + 1])
    
    # (6) Iteration number controller
    iter <- iter + 1
    if (max(grad_norm_z, grad_norm_alpha, rel_nll) <= rel_tol){
      stopflag = TRUE
      print(paste("Z: separate algorithm converged at iteration", iter,"!"))
    } else if (iter >= num_iter){
      stopflag = TRUE
      print("Z: separate algorithm achieved maximum iteration number!")
    } else{
      nll_collect[iter + 1] = nll_1
    }
    # print(iter)
  }
  
  return (list(Z_0 = Z_0,
               alpha_0 = alpha_0,
               nll = nll_collect))
}

eta_estimation = function(A_true, M, v_0, tau_v, num_iter, rel_tol){
  n = dim(v_0)[1]
  ones = rep(1, n)
  A_abs = abs(A_true)
  A_pos = matrix(0, nrow = n, ncol = n)
  A_pos[A_true > 0] = A_true[A_true > 0]
  
  Eta <- v_0 %*% t(v_0)
  
  stopflag = FALSE
  iter = 0
  
  nll_collect = rep(0, num_iter) # negative log-likelihood
  tsr = (A_pos * Eta + A_abs * log(1 - sigmoid(Eta))) * M
  nll_collect[1] = -sum(tsr[upper.tri(tsr)]) / n^2
  
  while(!stopflag){
    # (1) Update v
    v_1 = v_0 + 2 * tau_v * (((A_pos - A_abs * sigmoid(Eta)) * M) %*% v_0)  # Gradient descent step
    
    # (2) Prepare for the next iteration
    grad_norm = sum((v_1 - v_0)^2)
    v_0 <- v_1
    Eta <- v_0 %*% t(v_0)
    
    # (3) Compute Negative Log-likelihood
    tsr = (A_pos * Eta + A_abs * log(1 - sigmoid(Eta))) * M
    nll_1 = -sum(tsr[upper.tri(tsr)]) / n^2
    # print(nll_1)
    rel_nll = abs(nll_1 - nll_collect[iter + 1]) / abs(nll_collect[iter + 1])
    
    # (4) Iteration number controller
    iter <- iter + 1
    if (rel_nll <= rel_tol){
      stopflag = TRUE
      print(paste("v: separate algorithm converged at iteration", iter,"!"))
    } else if (iter >= num_iter){
      stopflag = TRUE
      print("v: separate algorithm achieved maximum iteration number!")
    } else{
      nll_collect[iter + 1] = nll_1
    }
    # print(iter)
  }
  
  return (list(v_0 = v_0,
               nll = nll_collect))
}

