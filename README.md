# SignedLSM

This repository implements the latent space approach for signed networks.

function_generate_data.R : Generate simulation data. 

function_separate_PGD.R : Projected gradient descent algorithms for the separate estimation method

function_joint_PGD.R : Projected gradient descent algorithm for the joint estimation method

util.R : includes the sigmoid function and the function to adjust for the optimal orthonormal transformation when evaluating the estimation error.

separate_estimation.R : Run the separate estimation method with different initializations. We implement two initialization: (1) usvt_init: universal singular value thresholding in Algorithm S3 and S4 in Supplymental Material; (2) random_init: random initialization. We find that the USVT initialization algorithms and random initialization achieve similar estimation errors for Z and v when the algorithms converge, while the proposed algorithms requires fewer iterations to converge.

one_step_joint_estimation.R : Run the one-step joint estimation

joint_estimation.R : Run the joint estimation method initialized by the estimates from the separate estimation.

network_density.R : Compare three estimation methods when varying the network density

sign_ratio.R : Compare three estimation methods when varying the ratio of positive and negative edges.



