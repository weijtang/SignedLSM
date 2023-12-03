# SignedLSM

This repo provides implementation for the latent space approach for signed networks in the following manuscript:

[Population-Level Balance for Signed Networks](https://arxiv.org/abs/2309.00193)

[Weijing Tang](https://sites.google.com/andrew.cmu.edu/weijingtang/home) and [Ji Zhu](http://dept.stat.lsa.umich.edu/~jizhu/). 

## Getting Started
### Prerequisites
The implementation of the proposed method is build on top of the following packages in [R](https://www.r-project.org/) (version 4.0.2).  

- 'pracma'
- 'scales'
- 'ggplot2'
- 'ggpubr'
- 'dplyr'


## Usage
[function_generate_data.R](https://github.com/weijtang/SignedLSM/blob/master/function_generate_data.R): Generate simulation data. 

[function_separate_PGD.R](https://github.com/weijtang/SignedLSM/blob/master/function_separate_PGD.R): Projected gradient descent algorithms for the separate estimation method

[function_joint_PGD.R](https://github.com/weijtang/SignedLSM/blob/master/function_joint_PGD.R): Projected gradient descent algorithm for the joint estimation method

[util.R](https://github.com/weijtang/SignedLSM/blob/master/util.R): Includes the sigmoid function, the function to adjust for the optimal orthonormal transformation when evaluating the estimation error, and the initialization algorithm based on UUniversal Singular Value Thresholding (USVT).

[separate_estimation.R](https://github.com/weijtang/SignedLSM/blob/master/separate_estimation.R): Run the separate estimation method with different initializations. We implement two initialization: (1) usvt_init: universal singular value thresholding in Algorithm S3 and S4 in Supplymental Material; (2) random_init: random initialization. We find that the USVT initialization algorithms and random initialization achieve similar estimation errors for Z and v when the algorithms converge, while the proposed algorithms requires fewer iterations to converge.

[one_step_joint_estimation.R](https://github.com/weijtang/SignedLSM/blob/master/one_step_joint_estimation.R): Run the one-step joint estimation

[joint_estimation.R](https://github.com/weijtang/SignedLSM/blob/master/joint_estimation.R): Run the joint estimation method initialized by the estimates from the separate estimation.

[summary.R](https://github.com/weijtang/SignedLSM/blob/master/summary.R): Compare three estimation methods when varying the network size $n$ and the dimension of the latent space $k$. This script can reproduce Figure 3 in the main text and Figure S3 in the Supplemental Material. 

[network_density.R](https://github.com/weijtang/SignedLSM/blob/master/network_density.R): Compare three estimation methods when varying the network density. This script can reproduce Figure S4 in the Supplemental Material.

[sign_ratio.R](https://github.com/weijtang/SignedLSM/blob/master/sign_ratio.R): Compare three estimation methods when varying the ratio of positive and negative edges. This script can reproduce Figure 4 in the main text.

[init_compare.R](https://github.com/weijtang/SignedLSM/blob/master/init_compare.R): Analyze the impact of initialization on the estimation errors and the number of iterations till convergenc. This script can reproduce Figure S1 in the Supplemental Material.


## Contact

Weijing Tang - weijingt@andrew.cmu.edu
