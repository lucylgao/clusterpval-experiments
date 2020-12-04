library(ggplot2)
library(patchwork)
library(fastcluster)

set.seed(1)
n <- 100
q <- 2
sig <- 1
X <- data.frame(matrix(rnorm(n*q, sd=sig), n, q)) 
colnames(X) <- c("Feat1", "Feat2")
X$clusters <- as.factor(cutree(hclust(dist(X)^2, method="average"), 3))

centroids <- aggregate(cbind(Feat1, Feat2)~clusters, X, mean)

p1 <- ggplot() + geom_point(aes(x=Feat1, y=Feat2, col=clusters), alpha=0.5, cex=3, data=X) + 
  xlab("Feature 1") + ylab("Feature 2") + 
  geom_point(aes(x=Feat1, y=Feat2, col=clusters), shape=17, cex=6,data=centroids) + 
  scale_colour_discrete(name="", labels=c(expression(bar(x)[hat(C)[1]]), 
                                          expression(bar(x)[hat(C)[2]]), 
                                          expression(bar(x)[hat(C)[3]]))) + 
  ggtitle("(a) Data") + xlim(c(-3, 3)) + ylim(c(-3, 3)) + 
  theme_bw(base_size=17) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_rect(colour = "black", size=1), 
        plot.margin = margin(0, 0.1, 0.1, 0.1, "cm")) + 
  theme(legend.position="left") + guides(color=guide_legend(nrow=3,byrow=TRUE))

load("../simulation-results/naive-type1-n100-q2-sig1.Rdata")

p2 <- ggplot() + stat_qq(aes(sample=ev[ev$Method == "average-Z-test-K-3", ][1:2000, ]$pval), distribution=qunif) + 
  geom_abline(slope=1, intercept=0, col="red") + xlab("Theoretical Quantiles") + 
  ylab("Empirical Quantiles") + ggtitle("(b) Wald test") + 
  xlim(c(0, 1)) + ylim(c(0, 1)) + theme_classic(base_size=17) 
p3 <- ggplot() + stat_qq(aes(sample=ev[ev$Method == "average-iso-test-K-3", ][1:2000, ]$pval), distribution=qunif) + 
  geom_abline(slope=1, intercept=0, col="red") + xlab("Theoretical Quantiles") + 
  ylab("Empirical Quantiles") + ggtitle("(c) Our proposed test") + 
  xlim(c(0, 1)) + ylim(c(0, 1)) + theme_classic(base_size=17)

ggsave(p1+p2+p3, file="../figures/Figure1.pdf", width=12, height=3)

