setwd("~/SignedLSM")
require(scales)
############### separate estimation ###############
load("./res/synthetic/separate_estimation_usvt_init.RData")
sep_results = data.frame(results)
sep_results$X5 = NULL # remove iter_z
sep_results$X6 = NULL # remove iter_v
sep_results$X7 = NULL # remove runtime_z
sep_results$X8 = NULL # remove runtime_v
colnames(sep_results) = c("Theta", "Eta", "Z", "v", "n", "k", "seed")
sep_results$Algorithm = rep("Separate Estimation", dim(sep_results)[1])

############### one-step-joint estimation ###############
load("./res/synthetic/one_step_joint_estimation_usvt_init.RData")
one_step_results = data.frame(results)
colnames(one_step_results) = c("Theta", "Eta", "Z", "v", "n", "k", "seed")
one_step_results$Algorithm = rep("One-step-joint Estimation", dim(one_step_results)[1])

############### joint estimation ###############
load("./res/synthetic/joint_estimation_usvt_init.RData")
joint_results = data.frame(results)
colnames(joint_results) = c("Theta", "Eta", "Z", "v", "n", "k", "seed")
joint_results$Algorithm = rep("Joint Estimation", dim(joint_results)[1])

results = rbind(sep_results, one_step_results, joint_results)

library(ggplot2)
library(ggpubr)
results$n = as.factor(results$n)
results$k = as.factor(results$k)
results$Algorithm = factor(results$Algorithm, 
                           levels = c("Joint Estimation", 
                                      "One-step-joint Estimation", 
                                      "Separate Estimation"))
######### fix k=2 ###############
pp = list()
i = 0
for (var in colnames(results)[c(3,4)]){
  p <- ggplot(results[results$k==2,], aes_string(x = "n", y = var, fill="Algorithm")) + 
    geom_boxplot(lwd=0.2, outlier.size =0.3) +
    scale_y_continuous(trans = "log2", 
                       breaks = c(0.03,0.06,0.12,0.24)) +
    xlab("Number of Nodes") + ylab(paste("Relative Error - ", var, sep = ""))+
    theme(text = element_text(size=11))
  i = 2* i+1
  pp[[i]] = p
}
p <- ggplot(results[results$k==2,], aes_string(x = "n", y = "Theta", fill="Algorithm")) + 
  geom_boxplot(lwd=0.2, outlier.size =0.3) +
  scale_y_continuous(trans = "log2", 
                     breaks = c(0.03,0.06,0.12,0.24)) +
  xlab("Number of Nodes") + ylab(expression(paste("Relative Error - ", Theta, sep = "")))+
  theme(text = element_text(size=11))
pp[[2]] = p
p <- ggplot(results[results$k==2,], aes_string(x = "n", y = "Eta", fill="Algorithm")) + 
  geom_boxplot(lwd=0.2, outlier.size =0.3) +
  scale_y_continuous(trans = "log2", 
                     breaks = c(0.03,0.06,0.12,0.24)) +
  xlab("Number of Nodes") + ylab(expression(paste("Relative Error - ", eta, sep = "")))+
  theme(text = element_text(size=11))
pp[[4]] = p
p = ggarrange(plotlist=pp, ncol=2, nrow=2, common.legend = TRUE, legend="top")
# pdf("./plots/separate_estimation_comparizon.pdf", width=7, height = 3.8)
pdf("./plots/all_estimation_comparison_k_2.pdf")
print(p)
dev.off()

######### fix n=2000 ###############
pp = list()
p <- ggplot(results[results$n==2000,], aes_string(x = "k", y = "Z", fill="Algorithm")) + 
  geom_boxplot(lwd=0.2, outlier.size =0.3) +
  scale_y_continuous(trans = "log2", 
                     breaks = c(0.04, round(0.04*sqrt(2),2), 0.08, 
                                round(0.08*sqrt(2),2),0.16,
                                round(0.16*sqrt(sqrt(2)),2),
                                round(0.16*sqrt(2),2), 
                                round(0.16*sqrt(2)*sqrt(sqrt(2)),2), 0.32))+#, 0.09, 0.11, 0.13)) +
  xlab("Dimension of the Latent Space") + ylab(paste("Relative Error - ", "Z", sep = ""))+
  theme(text = element_text(size=11))
pp[[1]] = p
p <- ggplot(results[results$n==2000,], aes_string(x = "k", y = "Theta", fill="Algorithm")) + 
  geom_boxplot(lwd=0.2, outlier.size =0.3) +
  scale_y_continuous(trans = "log2", 
                     breaks = c(0.04, round(0.04*sqrt(2),2), 0.08, 
                                round(0.08*sqrt(2),2),0.16,
                                round(0.16*sqrt(sqrt(2)),2),
                                round(0.16*sqrt(2),2), 
                                round(0.16*sqrt(2)*sqrt(sqrt(2)),2), 0.32))+#, 0.09, 0.11, 0.13)) +
  xlab("Dimension of the Latent Space") + ylab(expression(paste("Relative Error - ", Theta, sep = "")))+
  theme(text = element_text(size=11))
pp[[2]] = p
p = ggarrange(plotlist=pp, ncol=2, nrow=1, common.legend = TRUE, legend="top")
pdf("./plots/all_estimation_comparison_n_2000.pdf", width=7, height = 3.8)
print(p)
dev.off()

