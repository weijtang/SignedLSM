setwd("~/SignedLSM")
source("function_generate_data.R")
source("function_separate_PGD.R")
source("util.R")

# setup 
n_list = c(500, 1000, 2000, 4000)
k_list = c(2, 4, 8)
seed_list = 1:20
rel_tol = 1e-5
mode = "usvt_init" 

results = NULL
for (n in n_list){
  ones = rep(1,n)
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
      
      # init step: estimate z
      if (file.exists(paste("./res/synthetic/separate_z_", mode, file_seq, sep=""))){
        load(paste("./res/synthetic/separate_z_", mode, file_seq, sep=""))
      } else{
        print("Please run separate estimation!")
      }
      
      # init step: estimate v
      if (file.exists(paste("./res/synthetic/separate_v_", mode, file_seq, sep=""))){
        load(paste("./res/synthetic/separate_v_", mode, file_seq, sep=""))
      } else{
        print("Please run separate estimation!")
      }
      
      # one-step joint: estimate w, u through linear regression, then update z
      reg = lm(v_hat~1+Z_hat)
      u_hat = coef(reg)[[1]]
      w_hat = as.vector(coef(reg)[2:(k+1)])
      v_hat = Z_hat %*% w_hat + u_hat * ones
      Eta_hat <- v_hat %*% t(v_hat)
      
      # evaluation of v
      cls = df$Z_true %*% df$w_true + df$u_true * ones
      Eta <- cls %*% t(cls)
      e_eta = norm(Eta_hat - Eta, "F") / norm(Eta, "F")
      otg2 = oppo(v_hat, cls) 
      e_v = otg2$error / norm(cls, "F")
      
      # one-step joint: estimate w, u through linear regression, then update z
      # step size
      tau = 1.0
      tau_z = tau / max(norm(Z_hat, "2")^2, norm(v_hat, "F")^2)
      lambda = 0.5
      Theta_hat = alpha_hat %*% t(ones) + ones %*% t(alpha_hat) + Z_hat %*%  t(Z_hat)
      temp = (abs(df$A_true) * ((df$A_true+1)/2 - sigmoid(Eta_hat)))  %*% v_hat
      Z_hat = Z_hat + 2 * tau_z * ((1-lambda) * (abs(df$A_true) - sigmoid(Theta_hat)) %*% Z_hat
                                   + lambda * temp %*% t(w_hat))
      Z_hat = scale(Z_hat, center = TRUE, scale = FALSE)  # Centering Z
      
      # evaluation of v
      Theta = df$alpha_true %*% t(ones) + ones %*% t(df$alpha_true) + df$Z_true %*%  t(df$Z_true)
      Theta_hat = alpha_hat %*% t(ones) + ones %*% t(alpha_hat) + Z_hat %*%  t(Z_hat)
      e_theta = norm(Theta_hat - Theta, "F") / norm(Theta, "F")
      otg = oppo(Z_hat, df$Z_true) 
      e_Z = otg$error / norm(df$Z_true, "F")
      
      results = rbind(results, c(e_theta, e_eta, e_Z, e_v, n, k, seed))
    }
  }
}
save(results, file = paste("./res/synthetic/one_step_joint_estimation_", mode, ".RData", sep=""))