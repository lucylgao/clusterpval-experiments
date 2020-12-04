norm_vec <- function(x) sqrt(sum(x^2))

multivariate_Z_test <- function(X, linkage, K, k1, k2, sig) { 
  q <- ncol(X)
  hcl <- fastcluster::hclust(dist(X)^2, method=linkage)
  hcl_at_K <- cutree(hcl, K)
  
  diff_means <- colMeans(X[hcl_at_K == k1, , drop=F]) - colMeans(X[hcl_at_K == k2, , drop=F])
  stat <- norm_vec(diff_means) 
  n1 <- sum(hcl_at_K == k1)  
  n2 <- sum(hcl_at_K == k2) 
  squared_norm_nu <- 1/n1 + 1/n2
  scale_factor <- squared_norm_nu*sig^2
  
  pval <- 1 - pchisq(stat^2/scale_factor, df=q)
  return(list(stat=stat, pval=pval))
}

sample_split_Z_test <- function(X, linkage, K, k1, k2, sig) { 
  n <- nrow(X) 
  q <- ncol(X)
  X.train <- X[1:(n/2), ]
  X.test <-  X[(n/2+1):n, ]
  
  hcl <- fastcluster::hclust(dist(X.train)^2, method=linkage)
  
  new_class <- class::knn(X.train, X.test, cutree(hcl, K), k=K)
  
  if((sum(new_class == k1) > 0 ) & (sum(new_class == k2) > 0 )) { 
    diff_means <- colMeans(X.test[new_class == k1, , drop=F]) - 
      colMeans(X.test[new_class == k2, , drop=F])
    stat <- norm_vec(diff_means) 
    n1 <- sum(new_class == k1)  
    n2 <- sum(new_class == k1)  
    squared_norm_nu <- 1/n1 + 1/n2
    scale_factor <- squared_norm_nu*sig^2
    
    pval <- 1 - pchisq(stat^2/scale_factor, df=q)
    return(list(stat=stat, pval=pval))
  } else { 
    return(list(stat=NA, pval=1))
  }
  
}

make_sample_split <- function(link, K) {
  new_method(name = sprintf("%s-split-Z-test-K-%s", link, K),
             label = sprintf("sample splitting Z test (%s, K=%s)", link, K),
             settings = list(link = link, K=K),
             method = function(model, draw, link, K) {
               X <- draw$data
               
               pair <- sample(1:K, 2)
               k1 <- pair[1] 
               k2 <- pair[2]
               
               results <- sample_split_Z_test(X, link, K, k1, k2, model$sig)
               
               return(list(results=results, effect=(-1), recover=(-1), npair=c(-1, -1)))
             })
}

make_multivariate_Z_test <- function(link, K) {
  new_method(name = sprintf("%s-Z-test-K-%s", link, K),
             label = sprintf("Multivariate Z test (%s, K=%s)", link, K),
             settings = list(link = link, K=K),
             method = function(model, draw, link, K) {
               X <- draw$data
               
               pair <- sample(1:K, 2)
               k1 <- pair[1] 
               k2 <- pair[2]
               
               results <- multivariate_Z_test(X, link, K, k1, k2, model$sig)
               
               hcl <- fastcluster::hclust(dist(X)^2, method=link)
               hcl_at_K <- cutree(hcl, K)
               n1 <- sum(hcl_at_K == k1)  
               n2 <- sum(hcl_at_K == k2) 
               
               true_means <- draw$means[draw$clusters, ]
               effect <- norm_vec(colMeans(true_means[hcl_at_K == k1, , drop=F]) - colMeans(true_means[hcl_at_K == k2, , drop=F]))/model$sig
               crosstab <- table(hcl_at_K, draw$clusters)
               recover <- ifelse(sum(crosstab != 0) == 3, 1, 0)
               npair <- c(n1, n2)
               
               return(list(results=results, effect=effect, recover=recover, npair=npair))
             })
}

make_iso <- function(link, K) {
  new_method(name = sprintf("%s-iso-test-K-%s", link, K),
             label = sprintf("Isotropic (%s, K=%s)", link, K),
             settings = list(link = link, K=K),
             method = function(model, draw, link, K) {
               X <- draw$data
               
               pair <- sample(1:K, 2)
               k1 <- pair[1]
               k2 <- pair[2]
               
               library(fastcluster)
               hcl <- hclust(dist(X)^2, method=link)
               
               if(link == "complete") {
                 results <- test_complete_hier_clusters_approx(X, hcl, K, k1, k2, sig=model$sig)
               } else {
                 results <- test_hier_clusters_exact(X, link, hcl, K, k1, k2, sig=model$sig)
               }
               
               return(list(results=results))
             })
}

multivariate_Z_test_methods <- sapply(c("single", "average", "centroid", "complete"), 
                                      function(link) make_multivariate_Z_test(link, 3))

sample_split_methods <- sapply(c("single", "average", "centroid", "complete"), 
                                      function(link) make_sample_split(link, 2))

selective_methods <- sapply(c("single", "average", "centroid", "complete"),
                            function(link) make_iso(link, 3))