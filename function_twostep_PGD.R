Z_estimation = function(A_true, Z_0, alpha_0, tau_z, tau_alpha, num_iter, rel_tol){
  n = dim(Z_0)[1]
  k = dim(Z_0)[2]
  ones = rep(1, n)
  A_abs = abs(A_true)
  
  Theta = alpha_0 %*% t(ones) + ones %*% t(alpha_0) + Z_0 %*% t(Z_0)
  
  stopflag = FALSE
  iter = 0
  
  nll_collect = rep(0, num_iter) # negative log-likelihood
  tsr = A_abs * Theta + log(1 - sigmoid(Theta))
  nll_collect[1] = -sum(tsr[upper.tri(tsr)])
  
  while(!stopflag){
    # (2) Update Z
    Z_1 = Z_0 + 2 * tau_z * (A_abs - sigmoid(Theta)) %*% Z_0  # Gradient descent step
    
    # (3) Update alpha
    alpha_1 = alpha_0 + 2 * tau_alpha * (A_abs - sigmoid(Theta)) %*% ones
    
    # (4) Projection
    Z_1 = scale(Z_1, center = TRUE, scale = FALSE)  # Centering Z
    
    # (5) Prepare for the next iteration    
    Z_0 <- Z_1 
    alpha_0 <- alpha_1
    Theta = alpha_0 %*% t(ones) + ones %*% t(alpha_0) + Z_0 %*%  t(Z_0)
    
    # (1) Compute Negative Log-likelihood
    tsr = A_abs * Theta + log(1 - sigmoid(Theta))
    nll_1 = -sum(tsr[upper.tri(tsr)])
    rel_nll = abs(nll_collect[iter + 1] - nll_1) / abs(nll_collect[iter + 1])
    
    # (6) Iteration number controller
    iter <- iter + 1
    if (rel_nll <= rel_tol){
      stopflag = TRUE
      print(paste("first step algorithm converged at iteration", iter,"!"))
    } else if (iter >= num_iter){
      stopflag = TRUE
      print("first step algorithm achieved maximum iteration number!")
    } else{
      nll_collect[iter + 1] = nll_1
    }
    # print(iter)
  }
  
  return (list(Z_0 = Z_0,
               alpha_0 = alpha_0,
               nll = nll_collect))
}

eta_estimation = function(A_true, v_0, tau_v, num_iter, rel_tol){
  n = dim(v_0)[1]
  ones = rep(1, n)
  A_abs = abs(A_true)
  A_pos = matrix(0, nrow = n, ncol = n)
  A_pos[A_true > 0] = A_true[A_true > 0]
  
  Eta <- v_0 %*% t(v_0)
  
  stopflag = FALSE
  iter = 0
  
  nll_collect = rep(0, num_iter) # negative log-likelihood
  tsr = A_pos * Eta + A_abs * log(1 - sigmoid(Eta))
  nll_collect[1] = -sum(tsr[upper.tri(tsr)])
  
  while(!stopflag){
    # (2) Update v
    v_1 = v_0 + 2 * tau_v * (A_pos - A_abs * sigmoid(Eta)) %*% v_0  # Gradient descent step
    
    # (6) Prepare for the next iteration    
    v_0 <- v_1
    Eta <- v_0 %*% t(v_0)
    
    # (1) Compute Negative Log-likelihood
    tsr = A_pos * Eta + A_abs * log(1 - sigmoid(Eta))
    nll_1 = -sum(tsr[upper.tri(tsr)])
    rel_nll = abs(nll_1 - nll_collect[iter + 1]) / abs(nll_collect[iter + 1])
    
    # (6) Iteration number controller
    iter <- iter + 1
    if (rel_nll <= rel_tol){
      stopflag = TRUE
      print(paste("second-step algorithm converged at iteration", iter,"!"))
    } else if (iter >= num_iter){
      stopflag = TRUE
      print("second-step algorithm achieved maximum iteration number!")
    } else{
      nll_collect[iter + 1] = nll_1
    }
    # print(iter)
  }
  
  return (list(v_0 = v_0,
               nll = nll_collect))
}

