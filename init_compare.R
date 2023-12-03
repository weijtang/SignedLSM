setwd("~/SignedLSM")
require(scales)
source("util.R")

# setup 
n_list = c(500, 1000, 2000, 4000)
k = 2
seed_list = 1:20

for (mode in c("usvt_init", "random_init")){
  results = NULL
  for (n in n_list){
    ones = rep(1,n)
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
        print(paste("run separaete estimation: z_", mode, file_seq, sep=""))
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
        print(paste("run separaete estimation: v_", mode, file_seq, sep=""))
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
  }
  save(results, file = paste("./res/synthetic/separate_estimation_", mode, "_k_", k, ".RData", sep=""))
}


setwd("~/Dropbox (Personal)/Research/SignedNetwork/code copy")
require(scales)
source("util.R")

############### separate estimation usvt init ###############
load(paste("./res/synthetic/separate_estimation_", "usvt_init", "_k_", k, ".RData", sep=""))
usvt_result = data.frame(results)
usvt_result$X7 = NULL # remove runtime_z
usvt_result$X8 = NULL # remove runtime_v
colnames(usvt_result) = c("Theta", "Eta", "Z", "v", "Iter_Z", "Iter_v", "n", "k", "seed")
usvt_result$Initialization = rep("USVT", dim(usvt_result)[1])

############### separate estimation usvt init ###############
load(paste("./res/synthetic/separate_estimation_", "random_init", "_k_", k, ".RData", sep=""))
random_result = data.frame(results)
random_result$X7 = NULL # remove runtime_z
random_result$X8 = NULL # remove runtime_v
colnames(random_result) = c("Theta", "Eta", "Z", "v", "Iter_Z", "Iter_v", "n", "k", "seed")
random_result$Initialization = rep("Random", dim(random_result)[1])

results = rbind(usvt_result, random_result)

library(ggplot2)
library(ggpubr)
results$n = as.factor(results$n)
results$k = as.factor(results$k)
results$Initialization = factor(results$Initialization, levels = c("USVT", "Random"))
pp = list()
i = 0
for (var in colnames(results)[c(3,4)]){
  p <- ggplot(results, aes_string(x = "n", y = var, fill="Initialization")) + 
    geom_boxplot(lwd=0.2, outlier.size =0.3) +
    scale_y_continuous(trans = "log2", 
                       breaks = c(0.03,0.06,0.12,0.24)) +
    xlab("Number of Nodes") + ylab(paste("Relative Error - ", var, sep = ""))+
    theme(text = element_text(size=10))
  i = i+1
  pp[[i]] = p
}
p <- ggplot(results, aes_string(x = "n", y = "Iter_Z", fill="Initialization")) + 
  geom_boxplot(lwd=0.2, outlier.size =0.3) +
  # scale_y_continuous(trans = "log2", 
  #                    breaks = c(0.03,0.06,0.12,0.24)) +
  xlab("Number of Nodes") + ylab(paste("Num. of Iterations until Convergence - Z", sep = ""))+
  theme(text = element_text(size=10))
i = i+1
pp[[i]] = p

p <- ggplot(results, aes_string(x = "n", y = "Iter_v", fill="Initialization")) + 
  geom_boxplot(lwd=0.2, outlier.size =0.3) +
  # scale_y_continuous(trans = "log2", 
  #                    breaks = c(0.03,0.06,0.12,0.24)) +
  xlab("Number of Nodes") + ylab(paste("Num. of Iterations until Convergence - v", sep = ""))+
  theme(text = element_text(size=10))
i = i+1
pp[[i]] = p

p = ggarrange(plotlist=pp, ncol=2, nrow=2, common.legend = TRUE, legend="top")
pdf("./plots/separate_estimation_init_comparison_k_2.pdf")
print(p)
dev.off()