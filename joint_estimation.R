setwd("~/Dropbox (Personal)/Research/SignedNetwork/code")
source("function_generate_data.R")
source("function_PGD.R")
source("function_oppo.R")
source("function_twostep_PGD.R")
source("util.R")

# setup 
n_list = c(250, 500, 1000, 2000, 4000)
k = 2
num_iter = 200
seed_list = 1:10
rel_tol = 1e-5
random_init = FALSE

joint_estimation = function(n, k, num_iter, seed, random_init){
  #################################### generate data
  df = generate_data(n, k, seed)
  ones = rep(1,n)
  
  #################################### joint estimation
  # initialization
  if (random_init == TRUE){
    set.seed(100)
    Z_0 = matrix(rnorm(n * k, 0, 1), ncol = k)
    Z_0 = scale(Z_0, center = TRUE, scale = FALSE)
    alpha_0 = runif(n, -2, -1)
    w_0 = runif(k, -1, 1)
    u_0 = runif(1, -1, 1)
  } else{
    # first step: estimate z
    set.seed(100)
    Z_0 = matrix(rnorm(n * k, 0, 1), ncol = k)
    Z_0 = scale(Z_0, center = TRUE, scale = FALSE)
    alpha_0 = runif(n, -2, -1)
    
    tau = 1.0
    tau_z = tau / norm(Z_0, "2")^2
    tau_alpha = tau / (2*n)
    
    res = Z_estimation(df$A_true, Z_0, alpha_0, tau_z, tau_alpha, num_iter, rel_tol)
    Z_hat = res$Z_0
    alpha_hat = res$alpha_0
    
    # second step: estimate eta
    set.seed(100)
    v_0 = matrix(rnorm(n * 1, 0, 1), ncol = 1)
    tau = 1.0
    tau_v = tau / norm(v_0, "2")^2
    
    res = eta_estimation(df$A_true, v_0, tau_v, num_iter, rel_tol)
    v_hat = res$v_0
    
    # third step: estimate w, u through linear regression
    reg = lm(v_hat~1+Z_hat)
    u_hat = coef(reg)[[1]]
    w_hat = as.vector(coef(reg)[2:(k+1)])
    
    # use three_step estimation as the initialization for joint estimation
    Z_0 = Z_hat
    alpha_0 = alpha_hat
    w_0 = w_hat
    u_0 = u_hat
  }
  
  # step size
  if (random_init == TRUE){
    tau = 1.0
  } else{
    tau = 0.1
  }
  tau_z = tau / (norm(Z_0, "2")^2 + sum((Z_0 %*% w_0 + u_0 * rep(1, n))^2))
  tau_alpha = tau / (2*n)
  tau_w = tau / sum((t(Z_0) %*% (Z_0 %*% w_0 + u_0 * rep(1, n)))^2)
  tau_u = tau / (n * (n-1))
  
  res = PGD(df$A_true, Z_0, alpha_0, w_0, u_0, 
            tau_z, tau_alpha, tau_w, tau_u, num_iter, rel_tol)
  Z_hat = res$Z_0
  alpha_hat = res$alpha_0
  w_hat = res$w_0
  u_hat = res$u_0
  nll = res$nll
  # plot(1:num_iter, nll)

  # estimation error
  Theta = df$alpha_true %*% t(ones) + ones %*% t(df$alpha_true) + df$Z_true %*%  t(df$Z_true)
  Theta_hat = alpha_hat %*% t(ones) + ones %*% t(alpha_hat) + Z_hat %*%  t(Z_hat)
  e_theta = norm(Theta_hat - Theta, "F") / norm(Theta, "F")
  otg = oppo(Z_hat, df$Z_true) 
  e_Z = otg$error / norm(df$Z_true, "F")
  
  cls = df$Z_true %*% df$w_true + rep(df$u_true, n)
  Eta <- cls %*% t(cls)
  cls_hat = Z_hat %*% w_hat + rep(u_hat, n)
  Eta_hat <- cls_hat %*% t(cls_hat)
  e_eta = norm(Eta_hat - Eta, "F") / norm(Eta, "F")
  temp = oppo(w_hat, t(otg$adjust) %*% df$w_true)
  e_w = temp$error / norm(as.matrix(df$w_true), "F")
  e_u = abs(u_hat - df$u_true * temp$adjust)
  return (list(e_theta = e_theta^2,
               e_eta = e_eta^2,
               e_Z = e_Z^2, 
               e_w = e_w, 
               e_u= e_u))
}

results = NULL
for (n in n_list){
  for (seed in seed_list){
    print(seed)
    res = joint_estimation(n, k, num_iter, seed, random_init)
    results = rbind(results, c(res$e_theta, res$e_eta, res$e_Z, res$e_w, res$e_u, n, seed))
  }
  if (random_init == TRUE){
    save(results, file = paste("./res/joint_estimation_random_init_", n, ".RData"))
  } else{
    save(results, file = paste("./res/joint_estimation_three_step_init_", n, ".RData"))
  }
}
if (random_init == TRUE){
  save(results, file = "./res/joint_estimation_random_init.RData")
} else{
  save(results, file = "./res/joint_estimation_three_step_init.RData")
}
