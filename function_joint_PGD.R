# -----------------------------------------------------------------------------
# Description: Projected Gradient Descent Algorithm for Joint Estimation 
# -----------------------------------------------------------------------------
# Input:
# A_true    - Observed signed adjacency matrix
# M         - Missing index matrix
# Z_0       - Initialization of Z
# alpha_0   - Initialization of alpha
# w_0       - Initialization of w
# u_0       - Initialization of u
# tau_z     - Learning rate for Z
# tau_alpha - Learning rate for alpha
# tau_w     - Learning rate for w
# tau_u     - Learning rate for u
# lambda    - Learning rate for u
# num_iter  - Maximum iteration number
# rel_tol   - Tolerance for algorithm convergence
# -----------------------------------------------------------------------------
# Output: 
# Z_0       - Estimate for Z from joint estimation
# alpha_0   - Estimate for alpha from joint estimation
# w_0       - Estimate for w from joint estimation
# u_0       - Estimate for u from joint estimation
# nll       - Negative log likelihood for all iterates
# -----------------------------------------------------------------------------

PGD_wu = function(A_true, M, Z_0, alpha_0, w_0, u_0,  
               tau_z, tau_alpha, tau_w, tau_u, lambda, num_iter, rel_tol){
  n = dim(Z_0)[1]
  k = dim(Z_0)[2]
  ones = rep(1, n)
  A_abs = abs(A_true)
  A_pos = A_true * (A_true + 1) / 2

  Theta = alpha_0 %*% t(ones) + ones %*% t(alpha_0) + Z_0 %*%  t(Z_0)
  cls = Z_0 %*% w_0 + rep(u_0,n)
  Eta <- cls %*% t(cls)
  
  stopflag = FALSE
  iter = 0
  
  nll_collect = rep(0, num_iter) # negative log-likelihood
  tsr = (A_abs * Theta + log(1 - sigmoid(Theta)) + A_pos * Eta + A_abs * log(1 - sigmoid(Eta))) * M
  nll_collect[1] = -sum(tsr[upper.tri(tsr)]) / n^2
  
  while(!stopflag){
    # (1) Update Z
    temp = ((A_pos - A_abs * sigmoid(Eta)) * M) %*% cls
    Z_1 = Z_0 + 2 * tau_z * ((1-lambda)*((A_abs - sigmoid(Theta)) * M) %*% Z_0
                             + lambda* temp %*% t(w_0) )  # Gradient descent step
    # (2) Projection
    Z_1 = scale(Z_1, center = TRUE, scale = FALSE)  # Centering Z
    
    # (3) Update alpha
    alpha_1 = alpha_0 + 2 * tau_alpha * (1-lambda) * (((A_abs - sigmoid(Theta)) * M) %*% ones)
    
    # (4) Update w
    w_1 = w_0 + 2 * tau_w * lambda * (t(Z_0) %*% temp)
    
    # (5) Update u
    u_1 = u_0 + 2 * tau_u * lambda * (t(ones) %*% temp)
    
    # (6) Prepare for the next iteration   
    cls_1 = Z_1 %*% w_1 + rep(u_1,n)
    grad_norm_z = (norm(Z_1 - Z_0, "F"))^2
    grad_norm_alpha = (norm(alpha_1 - alpha_0, "F"))^2
    grad_norm_v = sum((cls_1 - cls)^2)
    Z_0 <- Z_1
    alpha_0 <- alpha_1
    w_0 <- w_1
    u_0 <- u_1
    Theta = alpha_0 %*% t(ones) + ones %*% t(alpha_0) + Z_0 %*%  t(Z_0)
    cls = cls_1
    Eta <- cls %*% t(cls)
    
    # (7) Compute Negative Log-likelihood
    tsr = (A_abs * Theta + log(1 - sigmoid(Theta)) + A_pos * Eta + A_abs * log(1 - sigmoid(Eta))) * M
    nll_1 = -sum(tsr[upper.tri(tsr)]) / n^2
    rel_nll = abs(nll_collect[iter + 1] - nll_1) / abs(nll_collect[iter + 1])
    
    # (8) Iteration number controller
    iter <- iter + 1
    if (max(rel_nll, grad_norm_z, grad_norm_alpha) <= rel_tol){
      stopflag = TRUE
      print(paste("joint estimation (w) algorithm converged at iteration", iter,"!"))
    } else if (iter >= num_iter){
      stopflag = TRUE
      print("joint estimation (w) algorithm achieved maximum iteration number!")
    } else{
      nll_collect[iter + 1] = nll_1
    }
    # print(iter)
  }
  
  return (list(Z_0 = Z_0,
               alpha_0 = alpha_0,
               w_0 = w_0,
               u_0 = u_0,
               nll = nll_collect))
}
