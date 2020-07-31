# Projected Gradient Descent --------------------------------------------------
# Description: PGD for the full model
# Date: February 16, 2020
# -----------------------------------------------------------------------------
# Input:
# A_list_true  - List of adjacency matrices
# U_0          - Initial U
# alpha_0_list - Initial aplha
# Lambda_list  - Initial Lambda
# eta_u        - Stepsize of U
# eta_alpha    - Stepsize of alpha
# eta_Lambda   - Stepsize of eta
# num_iter     - Number of iterations
# -----------------------------------------------------------------------------
# Output:
# U_0          - Estimated U
# alpha_list   - Estimated alpha
# Lambda_list  - Estimated Lambda
# nll          - List of negative log-likelihood per iteration
# -----------------------------------------------------------------------------

PGD = function(A_true, Z_0, alpha_0, w_0, u_0,  
               tau_z, tau_alpha, tau_w, tau_u, num_iter, rel_tol){
  n = dim(Z_0)[1]
  k = dim(Z_0)[2]
  ones = rep(1, n)
  A_abs = abs(A_true)
  A_pos = matrix(0, nrow = n, ncol = n)
  A_pos[A_true > 0] = A_true[A_true > 0]
  
  Theta = alpha_0 %*% t(ones) + ones %*% t(alpha_0) + Z_0 %*%  t(Z_0)
  cls = Z_0 %*% w_0 + rep(u_0, n)
  Eta <- cls %*% t(cls)
  
  stopflag = FALSE
  iter = 0
  
  nll_collect = rep(0, num_iter) # negative log-likelihood
  tsr = A_abs * Theta + log(1 - sigmoid(Theta)) + A_pos * Eta + A_abs * log(1 - sigmoid(Eta))
  nll_collect[1] = -sum(tsr[upper.tri(tsr)])
  
  while(!stopflag){
    # (2) Update Z
    Z_1 = Z_0 + 2 * tau_z * ((A_abs - sigmoid(Theta)) %*% Z_0
                             + (A_pos - A_abs * sigmoid(Eta)) %*% cls %*% t(w_0) )  # Gradient descent step
    # (3) Update w
    w_1 = w_0 + 2 * tau_w * t(Z_0) %*% (A_pos - A_abs * sigmoid(Eta)) %*% cls
    
    # (4) Update alpha
    alpha_1 = alpha_0 + 2 * tau_alpha * (A_abs - sigmoid(Theta)) %*% ones
    
    # (6) Update u
    u_1 = u_0 + 2 * tau_u * t(ones) %*% (A_pos - A_abs * sigmoid(Eta)) %*% cls
    
    # (7) Projection
    Z_1 = scale(Z_1, center = TRUE, scale = FALSE)  # Centering Z
    
    # (6) Prepare for the next iteration    
    Z_0 <- Z_1
    alpha_0 <- alpha_1
    w_0 <- w_1
    u_0 <- u_1
    Theta = alpha_0 %*% t(ones) + ones %*% t(alpha_0) + Z_0 %*%  t(Z_0)
    cls = Z_0 %*% w_0 + rep(u_0, n)
    Eta <- cls %*% t(cls)
    
    # (1) Compute Negative Log-likelihood
    tsr = A_abs * Theta + log(1 - sigmoid(Theta)) + A_pos * Eta + A_abs * log(1 - sigmoid(Eta))
    nll_1 = -sum(tsr[upper.tri(tsr)])
    rel_nll = abs(nll_collect[iter + 1] - nll_1) / abs(nll_collect[iter + 1])
    
    # (7) Iteration number controller
    iter <- iter + 1
    if (rel_nll <= rel_tol){
      stopflag = TRUE
      print(paste("joint estimation algorithm converged at iteration", iter,"!"))
    } else if (iter >= num_iter){
      stopflag = TRUE
      print("joint estimation algorithm achieved maximum iteration number!")
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
