library(ggplot2)
library(patchwork)
library(clusterpval)
library(fastcluster)

norm_vec <- function(x) sqrt(sum(x^2))

set.seed(1)
n <- 30
cl <- c(rep(1, 10), rep(2, 10), rep(3, 10))
mu <- rbind(c(0, 2), c(0, -2), c(sqrt(12), 0))
sig <- 1
q <- 2

X <- matrix(rnorm(n*q, sd=sig), n, q) + mu[cl, ]
hcl <- hclust(dist(X)^2, method="average")
clusters <- cutree(hcl, 3)

p1 <- ggplot(data.frame(X), aes(x=X1, y=X2, col=as.factor(clusters))) + 
  geom_point() + xlab("Feature 1") + ylab("Feature 2") + 
  theme_classic(base_size=18) + theme(legend.position="none") + 
  xlim(c(-4.2, 7)) + ylim(c(-4, 4.5)) + 
  scale_colour_manual(values=c("dodgerblue3", "rosybrown", "orange")) + 
  ggtitle(expression(paste("(a) Original data (", phi, "=4)",sep=""))) 

k1 <- 1
k2 <- 3
prop_k2 <- sum(clusters == k2)/(sum(clusters == k1) + sum(clusters == k2))
diff_means <- colMeans(X[clusters == k1, ]) - colMeans(X[clusters == k2, ])
phi <- sqrt(0)
stat <- sqrt(sum(diff_means^2))
Xphi <- X 
Xphi[clusters == k1, ] <- t(t(X[clusters == k1, ]) + prop_k2*(phi - stat)*diff_means/norm_vec(diff_means))
Xphi[clusters == k2, ] <- t(t(X[clusters == k2, ]) + (prop_k2 - 1)*(phi - stat)*diff_means/norm_vec(diff_means))
hcl_Xphi <- hclust(dist(Xphi)^2, method="average")
clusters_Xphi <- cutree(hcl_Xphi, 3)

p2 <- ggplot(data.frame(Xphi), aes(x=X1, y=X2, col=as.factor(clusters))) + 
  geom_point() + xlab("Feature 1") + ylab("Feature 2") + 
  theme_classic(base_size=18) + theme(legend.position="none") + 
  xlim(c(-4.2, 7)) + ylim(c(-4, 4.5)) + 
  scale_colour_manual(values=c("dodgerblue3", "rosybrown", "orange")) + 
  ggtitle(expression(paste("(b) Perturbed data (", phi, "=0)",sep=""))) 

phi2 <- sqrt(64)
Xphi2 <- X 
Xphi2[clusters == k1, ] <- t(t(X[clusters == k1, ]) + prop_k2*(phi2 - stat)*diff_means/norm_vec(diff_means))
Xphi2[clusters == k2, ] <- t(t(X[clusters == k2, ]) + (prop_k2 - 1)*(phi2 - stat)*diff_means/norm_vec(diff_means))
hcl_Xphi2 <- hclust(dist(Xphi2)^2, method="average")
p3 <- ggplot(data.frame(Xphi2), aes(x=X1, y=X2, col=as.factor(clusters))) + 
  geom_point() + xlab("Feature 1") + ylab("Feature 2") + 
  theme_classic(base_size=18) + theme(legend.position="none") + 
  xlim(c(-4.2, 7)) + ylim(c(-4, 4.5)) + 
  scale_colour_manual(values=c("dodgerblue3", "rosybrown", "orange")) + 
  ggtitle(expression(paste("(c) Perturbed data (", phi, "=8)",sep="")))

ggsave(p1 + p2 + p3, file="../figures/Figure3.pdf", width=14, height=4)

clusterpval::test_hier_clusters_exact(X, "average", 3, k1, k2, sig=sig, hcl=hcl)$trunc # This is S
