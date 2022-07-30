setwd("~/SignedLSM")
source("function_generate_data.R")
source("function_separate_PGD.R")
source("util.R")
library(filling)

# setup 
n_list = c(500, 1000, 2000, 4000)
k_list = c(2, 4, 8)
num_iter_coef = 100
seed_list = 1:20
rel_tol = 1e-5
mode = "usvt_init" # "usvt_init" or "random_init"

results = NULL
for (n in n_list){
  ones = rep(1,n)
  M = matrix(1, n, n)
  num_iter = ceiling(num_iter_coef * log(n))
  for (k in k_list){
    for (seed in seed_list){
      file_seq = paste("_n_", n, "_k_", k, "_seed_", seed, ".RData", sep="")
      print(file_seq)
      # generate data
      data_f = paste("./data/synthetic/data", file_seq, sep="")
      if (file.exists(data_f)){
        load(data_f)
      } else{
        df = generate_data(n, k, seed, 0, 0)
        save(df, file = data_f)
      }
      
      # first step: estimate z
      if (file.exists(paste("./res/synthetic/separate_z_", mode, file_seq, sep=""))){
        load(paste("./res/synthetic/separate_z_", mode, file_seq, sep=""))
      } else{
        # initialization
        if (mode == "random_init"){
          set.seed(100)
          Z_0 = matrix(rnorm(n * k, 0, 1), ncol = k)
          Z_0 = scale(Z_0, center = TRUE, scale = FALSE)
          alpha_0 = runif(n, -2, -1)
        } else if (mode == "usvt_init"){
          P_hat = fill.USVT(abs(df$A_true), eta = 0.01)$X
          P_hat[which(P_hat>0.49)] = 0.49# project Q into [0.05, 0.49]
          P_hat[which(P_hat<0.05)] = 0.05
          P_hat = (P_hat + t(P_hat)) / 2
          Theta_hat = log(P_hat / (1-P_hat))
          alpha_mean = mean(Theta_hat) / 2
          alpha_0 = rowMeans(Theta_hat) - alpha_mean
          J = diag(n)- matrix(-1/n, n, n)
          R = J %*% (Theta_hat - alpha_0 %*% t(ones) - ones %*% t(alpha_0)) %*% J
          svd_result = svd(R)
          Z_0 = svd_result$u[,1:k] %*% diag(sqrt(svd_result$d[1:k]))
        }
        # learning rate
        tau = 1.0
        tau_z = tau / norm(Z_0, "2")^2
        tau_alpha = tau / (2*n)
        # fitting
        start_time = proc.time()
        res = Z_estimation(df$A_true, M, Z_0, alpha_0, tau_z, tau_alpha, num_iter, rel_tol)
        runtime = proc.time() - start_time
        runtime_z = as.numeric(runtime[3])
        Z_hat = res$Z_0
        alpha_hat = res$alpha_0
        nll = res$nll
        save(Z_hat, alpha_hat, nll, runtime_z, file = paste("./res/synthetic/separate_z_", mode, file_seq, sep=""))
      }
      # evaluation
      Theta = df$alpha_true %*% t(ones) + ones %*% t(df$alpha_true) + df$Z_true %*%  t(df$Z_true)
      Theta_hat = alpha_hat %*% t(ones) + ones %*% t(alpha_hat) + Z_hat %*%  t(Z_hat)
      e_theta = norm(Theta_hat - Theta, "F") / norm(Theta, "F")
      otg = oppo(Z_hat, df$Z_true) 
      e_Z = otg$error / norm(df$Z_true, "F")
      iter_Z = sum(nll>0)
      
      # second step: estimate eta
      if (file.exists(paste("./res/synthetic/separate_v_", mode, file_seq, sep=""))){
        load(paste("./res/synthetic/separate_v_", mode, file_seq, sep=""))
      } else{
        # initialization
        if (mode == "random_init"){
          set.seed(100)
          v_0 = matrix(rnorm(n * 1, 0, 1), ncol = 1)
        } else if (mode == "near_true_init"){
          v_0 = df$Z_true %*% df$w_true + rep(df$u_true, n)
        } else if (mode == "usvt_init"){
          init_A = abs(df$A_true) * (df$A_true + 1) /2
          init_A[which(df$A_true == 0)] = NA
          Q_hat = fill.USVT(init_A, eta = 0.01)$X
          Q_hat[which(Q_hat>0.9)] = 0.9# project Q into [0.1, 0.9]
          Q_hat[which(Q_hat<0.1)] = 0.1
          Q_hat = (Q_hat + t(Q_hat)) / 2
          Eta_hat = log(Q_hat / (1-Q_hat))
          svd_result = svd(Eta_hat)
          v_0 = as.matrix(svd_result$u[,1] * sqrt(svd_result$d[1]))
        }
        # learning rate
        tau = 1.0
        tau_v = tau / norm(v_0, "F")^2
        # fitting
        start_time = proc.time()
        res = eta_estimation(df$A_true, M, v_0, tau_v, num_iter, rel_tol)
        runtime = proc.time() - start_time
        runtime_v = as.numeric(runtime[3])
        v_hat = res$v_0
        nll = res$nll
        save(v_hat, nll, runtime_v, file = paste("./res/synthetic/separate_v_", mode, file_seq, sep=""))
      }
      # evaluation
      cls = df$Z_true %*% df$w_true + rep(df$u_true, n)
      Eta <- cls %*% t(cls)
      Eta_hat <- v_hat %*% t(v_hat)
      e_eta = norm(Eta_hat - Eta, "F") / norm(Eta, "F")
      otg = oppo(v_hat, cls) 
      e_v = otg$error / norm(cls, "F")
      iter_v = sum(nll>0)
      
      results = rbind(results, c(e_theta, e_eta, e_Z, e_v, iter_Z, iter_v, runtime_z, runtime_v, n, k, seed))
    }
    save(results, file = paste("./res/synthetic/separate_estimation_", mode, "_n_", n, "_k_", k, ".RData", sep=""))
  }
}
save(results, file = paste("./res/synthetic/separate_estimation_", mode, ".RData", sep=""))