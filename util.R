sigmoid = function(x){
  1/(1+exp(-x))
}

# Orthogonal Procrustes Problem Optimal Value ---------------------------------
# Description: a.k.a. oppo value, Z and X are of the same shape
# -----------------------------------------------------------------------------
# Input:
# Z - numeric matrix
# X - numeric matrix
# -----------------------------------------------------------------------------
# Output:
# numeric - the optimal value of the Orthogonal Procrustes Problem
#           induced by the matrices Z and X
# -----------------------------------------------------------------------------
oppo = function(Z, X){
  svd_result = svd(t(X) %*% Z)
  A = svd_result$u
  B = svd_result$v
  R = A %*% t(B)
  return(list(error = norm(Z - X %*% R ,'F'), 
              adjust = R))
}