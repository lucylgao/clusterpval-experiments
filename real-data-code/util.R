norm_vec <- function(x) sqrt(sum(x^2))

#' Principal Orthogonal ComplEment Thresholding (POET) method
#'
#' This is a somewhat more computationally efficient version of
#'the code from the POET package on CRAN.
#'
#' @param Y: p by n matrix of raw data, where p is the dimensionality, n is the sample size
#' @param K: number of factors - recommendation is to count the number of very spiked (much
# larger than others) eigenvalues of the p by p sample covariance matrix of Y.
#' @param C: the positive constant for thresholding, user-specified
#' @param thres: choice of thresholding
#' @param matrix: the option of thresholding either correlation or covariance matrix
#'
#' @return
#' \item{SigmaY}{estimated p by p covariance matrix of y_t}
#' \item{SigmaU}{estimated p by p covariance matrix of u_t}
#' \item{factors}{estimated unobservable factors in a K by T matrix form}
#' \item{loadings}{estimated factor loadings in a p by K matrix form}
#'
#' @references Fan, J., Liao, Y., & Mincheva, M. (2013). Large covariance estimation by thresholding principal orthogonal complements. Journal of the Royal Statistical Society: Series B, 75(4).
POET <- function (Y, K = -Inf, C = -Inf, thres = "soft", matrix = "cor") {
  p = nrow(Y)
  n = ncol(Y)
  Y <- Y - t(t(apply(t(Y), 2, mean))) %*% matrix(1, 1, n)
  if (K == -Inf) {
    K1 = 0.25 * (POETKhat(Y)$K1HL + POETKhat(Y)$K2HL + POETKhat(Y)$K1BN +
                   POETKhat(Y)$K2BN)
    K = floor(K1) + 1
  }
  if (K > 0) {
    V <- eigen(t(Y) %*% Y)$vectors
    V = as.matrix(V)
    Dd <- eigen(t(Y) %*% Y)$values
    Dd = as.vector(Dd)
    W <- sort(diag(Dd), index.return = TRUE)$x
    W = as.matrix(W)
    Id <- sort(diag(Dd), index.return = TRUE)$ix
    Id = as.matrix(Id)
    F <- sqrt(n) * V[, 1:K]
    LamPCA = Y %*% F/n
    uhat = Y - LamPCA %*% t(F)
    Lowrank = LamPCA %*% t(LamPCA)
    rate = 1/sqrt(p) + sqrt((log(p))/n)
  }
  else {
    uhat = Y
    rate = sqrt((log(p))/n)
    Lowrank = matrix(0, p, p)
  }
  SuPCA = uhat %*% t(uhat)/n
  SuDiag = diag(diag(SuPCA))
  if (matrix == "cor") {
    R = solve(SuDiag^(1/2)) %*% SuPCA %*% solve(SuDiag^(1/2))
  }
  if (matrix == "vad") {
    R = SuPCA
  }
  if (C == -Inf) {
    C1 = POETCmin(Y, K, thres, matrix)
    C = C1 + 0.1
  }
  roottheta = array(0, dim = c(p, p))
  lambda = array(0, dim = c(p, p))
  for (i in 1:p) {
    for (j in 1:i) {
      temp = uhat[i, ] * uhat[j, ]
      roottheta[i, j] = sd(temp)
      lambda[i, j] = roottheta[i, j] * rate * C
      lambda[j, i] = lambda[i, j]
    }
  }
  Rthresh = matrix(0, p, p)
  if (thres == "soft") {
    for (i in 1:p) {
      for (j in 1:i) {
        if (abs(R[i, j]) < lambda[i, j] && j < i) {
          Rthresh[i, j] = 0
        }
        else {
          if (j == i) {
            Rthresh[i, j] = R[i, j]
          }
          else {
            Rthresh[i, j] = sign(R[i, j]) * (abs(R[i,
                                                   j]) - lambda[i, j])
          }
        }
        Rthresh[j, i] = Rthresh[i, j]
      }
    }
  }
  if (thres == "hard") {
    for (i in 1:p) {
      for (j in 1:i) {
        if (abs(R[i, j]) < lambda[i, j] && j < i) {
          Rthresh[i, j] = 0
        }
        else {
          Rthresh[i, j] = R[i, j]
        }
        Rthresh[j, i] = Rthresh[i, j]
      }
    }
  }
  if (thres == "scad") {
    for (i in 1:p) {
      for (j in 1:i) {
        if (j == i) {
          Rthresh[i, j] = R[i, j]
        }
        else {
          if (abs(R[i, j]) < lambda[i, j]) {
            Rthresh[i, j] = 0
          }
          else {
            if (abs(R[i, j]) < 2 * lambda[i, j]) {
              Rthresh[i, j] = sign(R[i, j]) * (abs(R[i,
                                                     j]) - lambda[i, j])
            }
            else {
              if (abs(R[i, j]) < 3.7 * lambda[i, j]) {
                Rthresh[i, j] = ((3.7 - 1) * R[i, j] -
                                   sign(R[i, j]) * 3.7 * lambda[i, j])/(3.7 -
                                                                          2)
              }
              else {
                Rthresh[i, j] = R[i, j]
              }
            }
          }
        }
        Rthresh[j, i] = Rthresh[i, j]
      }
    }
  }
  SigmaU = matrix(0, p, p)
  if (matrix == "cor") {
    SigmaU = SuDiag^(1/2) %*% Rthresh * SuDiag^(1/2)
  }
  if (matrix == "vad") {
    SigmaU = Rthresh
  }
  SigmaY = SigmaU + Lowrank
  result <- list(SigmaU = SigmaU, SigmaY = SigmaY, factors = t(F),
                 loadings = LamPCA)
  return(result)
}

wald_test <- function(X, link, K, k1, k2, iso=TRUE, sig=NULL, SigInv = NULL) { 
  n <- nrow(X)
  q <- ncol(X)
  
  hcl <- fastcluster::hclust(dist(X)^2, method=link)
  hcl_at_K <- cutree(hcl, K)
  
  diff_means <- colMeans(X[hcl_at_K == k1, , drop=F]) - colMeans(X[hcl_at_K == k2, , drop=F])
  n1 <- sum(hcl_at_K == k1)  
  n2 <- sum(hcl_at_K == k2) 
  squared_norm_nu <- 1/n1 + 1/n2
  
  if(iso) { 
    if(is.null(sig)) { 
      sig <- sqrt(sum(scale(X, scale=F)^2)/(n*q - q))
    }
    
    stat <- norm_vec(diff_means) 
    scale_factor <- squared_norm_nu*sig^2
  } else { 
    if(is.null(SigInv)) { 
      Sig <- cov(scale(X, scale=F))
      SigInv <- solve(Sig)
    }
    
    stat <- sqrt(as.numeric(t(diff_means)%*%SigInv%*%diff_means))
    scale_factor <- squared_norm_nu
  }

  
  pval <- 1 - pchisq(stat^2/scale_factor, df=q)
  return(list(stat=stat, pval=pval))
}
