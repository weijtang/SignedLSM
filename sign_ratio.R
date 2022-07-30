setwd("~/SignedLSM")
source("function_generate_data.R")
source("function_separate_PGD.R")
source("function_joint_PGD.R")
source("util.R")
library(filling)

############## separate estimation ################
# setup
n=2000
k = 4
gamma_list = seq(0,2,by=0.3)
num_iter_coef = 100
ones = rep(1,n)
M = matrix(1, n, n)
num_iter = ceiling(num_iter_coef * log(n))

seed_list = 1:20
rel_tol = 1e-5
mode = "usvt_init"

results = NULL
for (gamma in gamma_list){
  for (seed in seed_list){
    file_seq = paste("_n_", n, "_k_", k, "_gamma_", gamma,  "_seed_", seed, ".RData", sep="")
    print(file_seq)
    # generate data
    data_f = paste("./data/synthetic/data", file_seq, sep="")
    if (file.exists(data_f)){
      load(data_f)
    } else{
      df = generate_data(n, k, seed, 0, gamma)
      save(df, file = data_f)
    }
    pos_sign_ratio = sum(df$A_true>0)/sum(abs(df$A_true))
    print(pos_sign_ratio)

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
      } else if (mode == "near_true_init"){
        Z_0 = df$Z_true
        alpha_0 = df$alpha_true
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

    v_hat = scale(v_hat, center = TRUE, scale = FALSE)  # Centering v
    e_Zw = norm(v_hat - df$Z_true %*% df$w_true %*% otg$adjust, "F") / norm(df$Z_true %*% df$w_true, "F")

    results = rbind(results, c(e_theta, e_eta, e_Z, e_v, iter_Z, iter_v, runtime_z, runtime_v, n, k, seed, gamma, pos_sign_ratio, e_Zw))
  }
  save(results, file = paste("./res/synthetic/separate_estimation_", mode, "_n_", n, "_k_", k, "_gamma_", gamma, ".RData", sep=""))
}

############## one-step-joint estimation ################
# setup
n=2000
k = 4
gamma_list = seq(0,2,by=0.3)
lambda_list = seq(0.1,0.9, by=0.1)
num_iter_coef = 100
ones = rep(1,n)
M = matrix(1, n, n)
num_iter = ceiling(num_iter_coef * log(n))

seed_list = 1:20
rel_tol = 1e-5
mode = "usvt_init"

results = NULL
for (gamma in gamma_list){
  for (seed in seed_list){
    for (lambda in lambda_list){
      file_seq = paste("_n_", n, "_k_", k, "_gamma_", gamma,  "_seed_", seed, ".RData", sep="")
      print(file_seq)
      # generate data
      data_f = paste("./data/synthetic/data", file_seq, sep="")
      if (file.exists(data_f)){
        load(data_f)
      } else{
        df = generate_data(n, k, seed, 0, gamma)
        save(df, file = data_f)
      }
      pos_sign_ratio = sum(df$A_true>0)/sum(abs(df$A_true))

      # init step: estimate z
      if (file.exists(paste("./res/synthetic/separate_z_", mode, file_seq, sep=""))){
        load(paste("./res/synthetic/separate_z_", mode, file_seq, sep=""))
      } else{
        print("Please run separate estimation!")
        next
      }

      # init step: estimate v
      if (file.exists(paste("./res/synthetic/separate_v_", mode, file_seq, sep=""))){
        load(paste("./res/synthetic/separate_v_", mode, file_seq, sep=""))
      } else{
        print("Please run separate estimation!")
        next
      }

      # one-step joint: estimate w, u through linear regression, then update z
      reg = lm(v_hat~1+Z_hat)
      u_hat = coef(reg)[[1]]
      w_hat = as.vector(coef(reg)[2:(k+1)])
      v_hat = Z_hat %*% w_hat + rep(u_hat,n)
      Eta_hat <- v_hat %*% t(v_hat)

      # evaluation of v
      cls = df$Z_true %*% df$w_true + rep(df$u_true,n)
      Eta <- cls %*% t(cls)
      e_eta = norm(Eta_hat - Eta, "F") / norm(Eta, "F")
      otg2 = oppo(v_hat, cls)
      e_v = otg2$error / norm(cls, "F")

      v_hat = scale(v_hat, center = TRUE, scale = FALSE)  # Centering v
      e_Zw = norm(v_hat - df$Z_true %*% df$w_true %*% otg2$adjust, "F") / norm(df$Z_true %*% df$w_true, "F")

      # one-step joint: estimate w, u through linear regression, then update z
      # step size
      tau = 1.0
      tau_z = tau / max(norm(Z_hat, "2")^2, norm(v_hat, "F")^2)
      # lambda = 0.5
      Theta_hat = alpha_hat %*% t(ones) + ones %*% t(alpha_hat) + Z_hat %*%  t(Z_hat)
      temp = (abs(df$A_true) * ((df$A_true+1)/2 - sigmoid(Eta_hat)))  %*% v_hat
      Z_hat = Z_hat + 2 * tau_z * ((1-lambda) * (abs(df$A_true) - sigmoid(Theta_hat)) %*% Z_hat
                                   + lambda * temp %*% t(w_hat))
      Z_hat = scale(Z_hat, center = TRUE, scale = FALSE)  # Centering Z

      # evaluation of z
      Theta = df$alpha_true %*% t(ones) + ones %*% t(df$alpha_true) + df$Z_true %*%  t(df$Z_true)
      Theta_hat = alpha_hat %*% t(ones) + ones %*% t(alpha_hat) + Z_hat %*%  t(Z_hat)
      e_theta = norm(Theta_hat - Theta, "F") / norm(Theta, "F")
      otg = oppo(Z_hat, df$Z_true)
      e_Z = otg$error / norm(df$Z_true, "F")

      results = rbind(results, c(e_theta, e_eta, e_Z, e_v, seed, gamma, pos_sign_ratio, e_Zw, lambda))
    }
  }
  save(results, file = paste("./res/synthetic/one_step_joint_estimation_", mode, "_n_", n, "_k_", k, "_gamma_", gamma, "_tune_lambda", ".RData", sep=""))
}


############## joint estimation ################
# setup 
n=2000
k = 4
gamma_list = 1.2 #seq(0,2,by=0.3)
lambda_list = seq(0.1,0.9, by=0.1)
num_iter_coef = 100
ones = rep(1,n)
M = matrix(1, n, n)
num_iter = ceiling(num_iter_coef * log(n))

seed_list = 1:20
rel_tol = 1e-5
mode = "usvt_init" 

results = NULL
for (gamma in gamma_list){
  for (seed in seed_list){
    for (lambda in lambda_list){
      file_seq = paste("_n_", n, "_k_", k, "_gamma_", gamma,  "_seed_", seed, ".RData", sep="")
      lambda_file_seq = paste("_n_", n, "_k_", k, "_gamma_", gamma,  "_seed_", seed, "_lambda_", lambda,".RData", sep="")
      print(file_seq)
      # generate data
      data_f = paste("./data/synthetic/data", file_seq, sep="")
      if (file.exists(data_f)){
        load(data_f)
      } else{
        df = generate_data(n, k, seed, 0, gamma)
        save(df, file = data_f)
      }
      pos_sign_ratio = sum(df$A_true>0)/sum(abs(df$A_true))
      
      if (file.exists(paste("./res/synthetic/joint_", mode, lambda_file_seq, sep=""))){
        load(paste("./res/synthetic/joint_", mode, lambda_file_seq, sep=""))
      } else{
        # init step: estimate z
        if (file.exists(paste("./res/synthetic/separate_z_", mode, file_seq, sep=""))){
          load(paste("./res/synthetic/separate_z_", mode, file_seq, sep=""))
        } else{
          print("Please run separate estimation!")
          next
        }
        
        # init step: estimate v
        if (file.exists(paste("./res/synthetic/separate_v_", mode, file_seq, sep=""))){
          load(paste("./res/synthetic/separate_v_", mode, file_seq, sep=""))
        } else{
          print("Please run separate estimation!")
          next
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
        # lambda = 0.5
        ############## fitting #############
        res = PGD_wu(df$A_true, M, Z_0, alpha_0, w_0, u_0, 
                     tau_z, tau_alpha, tau_w, tau_u, lambda, num_iter, rel_tol)
        Z_hat = res$Z_0
        alpha_hat = res$alpha_0
        w_hat = res$w_0
        u_hat = res$u_0
        nll = res$nll
        save(Z_hat, alpha_hat, w_hat, u_hat, nll, file = paste("./res/synthetic/joint_", mode, lambda_file_seq, sep=""))
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
      
      v_hat = scale(v_hat, center = TRUE, scale = FALSE)  # Centering v
      e_Zw = norm(v_hat - df$Z_true %*% df$w_true %*% otg$adjust, "F") / norm(df$Z_true %*% df$w_true, "F")
      
      results = rbind(results, c(e_theta, e_eta, e_Z, e_v, seed, gamma, pos_sign_ratio, e_Zw, lambda))
    }
  }
  save(results, file = paste("./res/synthetic/joint_estimation_", mode, "_n_", n, "_k_", k, "_gamma_", gamma, "_tune_lambda", ".RData", sep=""))
}

################## summary #######################
require(scales)
library(dplyr)
library(ggplot2)
library(ggpubr)

n=2000
k = 4
mode = "usvt_init"
# load separate estimation
load(paste("./res/synthetic/separate_estimation_", mode, "_n_", n, "_k_", k, "_gamma_", 1.8, ".RData", sep=""))
result = data.frame(results)
result$X5 = NULL # remove iter_z
result$X6 = NULL # remove iter_v
result$X7 = NULL # remove runtime_z
result$X8 = NULL # remove runtime_v
result$X9 = NULL # remove n
result$X10 = NULL # remove k
colnames(result) = c("Theta", "Eta", "Z", "v", "seed", "gamma", "Positive_ratio_seed", "Zw")
result_sep = result %>% group_by(gamma) %>% summarise(Positive_ratio = mean(Positive_ratio_seed), 
                                                      Z_mean = mean(Z), Z_sd = sd(Z), 
                                                      v_mean = mean(Zw), v_sd = sd(Zw), .groups = 'drop')
result_sep$Algorithm = rep("Separate Estimation", dim(result_sep)[1])

# load one-step-joint estimation 
load(paste("./res/synthetic/one_step_joint_estimation_", mode, "_n_", n, "_k_", k, "_gamma_", 1.8, ".RData", sep=""))
result = data.frame(results)
colnames(result) = c("Theta", "Eta", "Z", "v", "seed", "gamma", "Positive_ratio_seed", "Zw")
result_one_step = result %>% group_by(gamma) %>% summarise(Positive_ratio = mean(Positive_ratio_seed), 
                                                           Z_mean = mean(Z), Z_sd = sd(Z), 
                                                           v_mean = mean(Zw), v_sd = sd(Zw), .groups = 'drop')
result_one_step$Algorithm = rep("One-step-joint Estimation", dim(result_one_step)[1])

# load joint estimation 
load(paste("./res/synthetic/joint_estimation_", mode, "_n_", n, "_k_", k, "_gamma_", 1.8, ".RData", sep=""))
result = data.frame(results)
colnames(result) = c("Theta", "Eta", "Z", "v", "seed", "gamma", "Positive_ratio_seed", "Zw")
result_joint = result %>% group_by(gamma) %>% summarise(Positive_ratio = mean(Positive_ratio_seed), 
                                                        Z_mean = mean(Z), Z_sd = sd(Z), 
                                                        v_mean = mean(Zw), v_sd = sd(Zw), .groups = 'drop')
result_joint$Algorithm = rep("Joint Estimation", dim(result_joint)[1])


result = rbind(result_sep, result_one_step,result_joint)
result$Algorithm = factor(result$Algorithm, 
                          levels = c("Joint Estimation", 
                                     "One-step-joint Estimation", 
                                     "Separate Estimation"))
pp = list()
p_z <- ggplot(result, aes(x=Positive_ratio, y=Z_mean, colour=Algorithm)) + 
  geom_errorbar(aes(ymin=Z_mean-Z_sd, ymax=Z_mean+Z_sd), width=.02) +
  geom_line(aes(linetype = Algorithm), size=0.5) +
  geom_point(aes(shape = Algorithm), size=1.5) +
  xlab("Proportion of Positive Edges") + ylab("Relative Error - Z")+
  theme(text = element_text(size=11)) + ylim(0.09, 0.12)
p_v <- ggplot(result, aes(x=Positive_ratio, y=v_mean, colour=Algorithm)) + 
  geom_errorbar(aes(ymin=v_mean-v_sd, ymax=v_mean+v_sd), width=.02) +
  geom_line(aes(linetype = Algorithm), size=0.5) +
  geom_point(aes(shape = Algorithm), size=1.5) +
  xlab("Proportion of Positive Edges") + ylab(expression(paste("Relative Error - ", v[cen], sep="")))+
  theme(text = element_text(size=11)) + ylim(0.05, 0.24)
pp[[1]] = p_z
pp[[2]] = p_v
p = ggarrange(plotlist=pp, ncol=2, nrow=1, common.legend = TRUE, legend="top")
p
pdf("./plots/sign_ratio_comparison_n_2000_k_4.pdf", width=7, height = 3.8)
print(p)
dev.off()


