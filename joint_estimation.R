setwd("~/SignedLSM")
source("function_generate_data.R")
source("function_joint_PGD.R")
source("util.R")

# setup 
n_list = c(500, 1000, 2000, 4000)
k_list = c(2, 4, 8)
num_iter_coef = 100
seed_list = 1:20
rel_tol = 1e-5
mode = "usvt_init"

results = NULL
for (n in n_list){
  ones = rep(1,n)
  M = matrix(1, n, n)
  num_iter = ceiling(num_iter_coef * log(n))
  for (k in k_list){
    for (seed in seed_list){
      file_seq = paste("_n_", n, "_k_", k, "_seed_", seed, ".RData", sep="")
      print(file_seq)
      # load data
      data_f = paste("./data/synthetic/data", file_seq, sep="")
      if (file.exists(data_f)){
        load(data_f)
      } else{
        df = generate_data(n, k, seed, 0, 0)
        save(df, file = data_f)
      }
      
      if (file.exists(paste("./res/synthetic/joint_", mode, file_seq, sep=""))){
        load(paste("./res/synthetic/joint_", mode, file_seq, sep=""))
      } else{
        ############# initialization ################
        # init step: estimate z
        if (file.exists(paste("./res/synthetic/separate_z_", mode, file_seq, sep=""))){
          load(paste("./res/synthetic/separate_z_", mode, file_seq, sep=""))
        } else{
          print("Please run separate estimation for initialization!")
        }
        
        # init step: estimate v
        if (file.exists(paste("./res/synthetic/separate_v_", mode, file_seq, sep=""))){
          load(paste("./res/synthetic/separate_v_", mode, file_seq, sep=""))
        } else{
          print("Please run separate estimation for initialization!")
        }
        # init step : u and w
        reg = lm(v_hat~1+Z_hat)
        u_hat = coef(reg)[[1]]
        w_hat = as.vector(coef(reg)[2:(k+1)])
        v_hat = Z_hat %*% w_hat + rep(u_hat,n)
        
        Z_0 = Z_hat
        alpha_0 = alpha_hat
        w_0 = w_hat
        u_0 = u_hat
        
        ############## step size #############
        tau = 1.0
        tau_z = tau / max(norm(Z_hat, "2")^2, norm(v_hat, "F")^2)
        tau_alpha = tau / (2*n)
        tau_w = tau / (norm(v_hat, "F")^2 * norm(Z_hat, "2")^2)
        tau_u = tau / (norm(v_hat, "F")^2 * n)
        lambda = 0.5
        
        ############## fitting #############
        res = PGD_wu(df$A_true, M, Z_0, alpha_0, w_0, u_0, 
                     tau_z, tau_alpha, tau_w, tau_u, lambda, num_iter, rel_tol)
        Z_hat = res$Z_0
        alpha_hat = res$alpha_0
        w_hat = res$w_0
        u_hat = res$u_0
        nll = res$nll
        save(Z_hat, alpha_hat, w_hat, u_hat, nll, file = paste("./res/synthetic/joint_", mode, file_seq, sep=""))
      }
      
      # estimation error
      Theta = df$alpha_true %*% t(ones) + ones %*% t(df$alpha_true) + df$Z_true %*%  t(df$Z_true)
      Theta_hat = alpha_hat %*% t(ones) + ones %*% t(alpha_hat) + Z_hat %*%  t(Z_hat)
      e_theta = norm(Theta_hat - Theta, "F") / norm(Theta, "F")
      otg = oppo(Z_hat, df$Z_true) 
      e_Z = otg$error / norm(df$Z_true, "F")
      
      cls = df$Z_true %*% df$w_true + rep(df$u_true, n)
      Eta <- cls %*% t(cls)
      v_hat = Z_hat %*% w_hat + rep(u_hat,n)
      Eta_hat <- v_hat %*% t(v_hat)
      e_eta = norm(Eta_hat - Eta, "F") / norm(Eta, "F")
      otg = oppo(v_hat, cls) 
      e_v = otg$error / norm(cls, "F")
      
      results = rbind(results, c(e_theta, e_eta, e_Z, e_v, n, k, seed))
    }
    save(results, file = paste("./res/synthetic/joint_estimation_", mode, "_n_", n, "_k_", k, ".RData", sep=""))
  }
}
save(results, file = paste("./res/synthetic/joint_estimation_", mode, ".RData", sep=""))