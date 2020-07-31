setwd("~/Dropbox (Personal)/Research/SignedNetwork/code")
source("function_generate_data.R")
source("function_oppo.R")
source("function_twostep_PGD.R")
source("util.R")

# setup 
n_list = c(250, 500, 1000, 2000, 4000)
k = 2
num_iter = 200
seed_list = 1:10
rel_tol = 1e-5

three_step_estimation = function(n, k, num_iter, seed){
  #################################### generate data
  df = generate_data(n, k, seed)
  ones = rep(1,n)
  
  #################################### three-step estimation
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
  nll = res$nll
  # plot(1:num_iter, nll)
  
  Theta = df$alpha_true %*% t(ones) + ones %*% t(df$alpha_true) + df$Z_true %*%  t(df$Z_true)
  Theta_hat = alpha_hat %*% t(ones) + ones %*% t(alpha_hat) + Z_hat %*%  t(Z_hat)
  e_theta = norm(Theta_hat - Theta, "F") / norm(Theta, "F")
  otg = oppo(Z_hat, df$Z_true) 
  e_Z = otg$error / norm(df$Z_true, "F")
  
  # second step: estimate eta
  set.seed(100)
  v_0 = matrix(rnorm(n * 1, 0, 1), ncol = 1)
  tau = 1.0
  tau_v = tau / norm(v_0, "2")^2
  
  res = eta_estimation(df$A_true, v_0, tau_v, num_iter, rel_tol)
  v_hat = res$v_0
  nll = res$nll
  # plot(1:num_iter, nll)
  
  # cls = df$Z_true %*% df$w_true + rep(df$u_true, n)
  # Eta <- df$c_true * cls %*% t(cls)
  # Eta_hat <- v_hat %*% t(v_hat)
  # e_eta = norm(Eta_hat - Eta, "F") / norm(Eta, "F")
  
  # third step: estimate w, u through linear regression
  reg = lm(v_hat~1+Z_hat)
  u_hat = coef(reg)[[1]]
  w_hat = as.vector(coef(reg)[2:(k+1)])
  
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
    res = three_step_estimation(n, k, num_iter, seed)
    results = rbind(results, c(res$e_theta, res$e_eta, res$e_Z, res$e_w, res$e_u, n, seed))
  }
  save(results, file = paste("./res/three_step_estimation_", n, ".RData"))
}
save(results, file = "./res/three_step_estimation.RData")
