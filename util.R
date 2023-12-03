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

# Universal Singular Value Thresholding ---------------------------------------
fill.USVT <- function(A, eta=0.01){
  # Parameter
  #   1. data check, dimension, and missing-1 matrix
  X = check_data(A)
  if (check_bycol(X)==FALSE){   warning("* fill.USVT : there exists at least one column full of missing entries.")}
  if (check_bycol(t(X))==FALSE){warning("* fill.USVT : there exists at least one row full of missing entries.")}
  # 2. eta
  if ((length(eta)>1)||(!is.numeric(eta))||(eta<=0)||(eta>=1)){
    stop("* fill.USVT : a parameter 'eta' should be in (0,1).")
  }
  
  ##############################################################
  # Preparation Step
  # 1. scaling : later the output must be rescaled back : bounded by abs <= 1
  a = min(X[!is.na(X)])
  b = max(X[!is.na(X)])
  # 2. size must be (mxn) with m <= n
  if (nrow(X) > ncol(X)){
    tflag = TRUE   # must be transposed later
    X     = t(X)
  } else {
    tflag = FALSE
  }
  # 3. size parameters
  m = nrow(X)
  n = ncol(X)
  # 4. rescaling of the data
  M = ((X - array((a+b)/2, c(m,n)))/((b-a)/2))
  
  ##############################################################
  # Main Computation as in Paper, under $1.2
  # (1) fill in the proxy Y
  Y = M
  Y[is.na(M)] = 0.0
  # (2) singular value decomposition
  svdY = svd(Y)
  s    = svdY$d
  # (3) proportion of observed values of X
  phat = (sum(!is.na(A))/length(A))
  # (4) thresholding value
  thrval = ((2+eta)*sqrt(n*phat))
  S      = (svdY$d >= phat)
  # (5) Define W
  if (sum(S)==1){
    W = (outer(svdY$u[,], svdY$v[,S])*svdY$d[S])/phat
  } else {
    W = (svdY$u[,S]%*%diag(svdY$d[S])%*%t(svdY$v[,S]))/phat
  }
  # (6) Alter Mhat
  Mhat = W
  Mhat[(W>1)]  = 1.0
  Mhat[(W< -1)] = -1.0
  # (7) re-scale back
  result = ((Mhat*((b-a)/2)) + array((a+b)/2, c(m,n)))
  # (8) control the transpose
  if (tflag==TRUE){
    result = t(result)
  }
  
  
  ##############################################################
  ## RETURN
  output = list()
  output$X = result
  return(output)
}