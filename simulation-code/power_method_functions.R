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
               hcl_at_K <- cutree(hcl, K)
               n1 <- sum(hcl_at_K == k1)  
               n2 <- sum(hcl_at_K == k2) 
               
               true_means <- draw$means[draw$clusters, ]
               effect <- norm_vec(colMeans(true_means[hcl_at_K == k1, , drop=F]) - colMeans(true_means[hcl_at_K == k2, , drop=F]))/model$sig
               npair <- c(n1, n2)
               crosstab <- table(hcl_at_K, draw$clusters)
               recover <- ifelse(sum(crosstab != 0) == 3, 1, 0)
               
               if(link == "complete") {
                  results <- test_complete_hier_clusters_approx(X, hcl, K, k1, k2, sig=model$sig)
               } else {
                  results <- test_hier_clusters_exact(X, link, hcl, K, k1, k2, sig=model$sig)
               }

               return(list(results=results, effect=effect, recover=recover, npair=npair))
             })
}


selective_test_methods <- sapply(c("single", "average", "centroid", "complete"),
                           function(link) make_iso(link, 3))

