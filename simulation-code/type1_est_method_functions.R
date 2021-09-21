norm_vec <- function(x) sqrt(sum(x^2))

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

make_iso_est <- function(link, K) {
  new_method(name = sprintf("%s-iso-est-test-K-%s", link, K),
             label = sprintf("Isotropic, estimated sigma (%s, K=%s)", link, K),
             settings = list(link = link, K=K),
             method = function(model, draw, link, K) {
               X <- draw$data
               n <- nrow(X)
               q <- ncol(X)
               
               pair <- sample(1:K, 2)
               k1 <- pair[1]
               k2 <- pair[2]
               
               train.id <- sample(1:n, n/2, replace=FALSE)
               test.id <- setdiff(1:n, train.id)
               est.sig <- sqrt(sum(scale(X[test.id, ], scale=FALSE)^2)/((n/2)*q - q))
               X.train <- X[train.id, ]
               means.train <- draw$means[draw$clusters, ][train.id, ]
          
               library(fastcluster)
               hcl <- hclust(dist(X.train)^2, method=link)
               hcl_at_K <- cutree(hcl, K)
               n1 <- sum(hcl_at_K == k1)  
               n2 <- sum(hcl_at_K == k2) 
               effect <- norm_vec(colMeans(means.train[hcl_at_K == k1, , drop=F]) - colMeans(means.train[hcl_at_K == k2, , drop=F]))
               npair <- c(n1, n2)

               if(link == "complete") {
                 results <- test_complete_hier_clusters_approx(X.train, hcl, K, k1, k2, sig=est.sig)
               } else {
                 results <- test_hier_clusters_exact(X.train, link, hcl, K, k1, k2, sig=est.sig)
               }
               
               return(list(results=results, effect=effect, npair=npair, sd=est.sig))
             })
}

selective_methods <- sapply(c("single", "average", "centroid", "complete"),
                            function(link) make_iso(link, 3))

selective_methods_est <- sapply(c("single", "average", "centroid", "complete"),
                            function(link) make_iso_est(link, 3))