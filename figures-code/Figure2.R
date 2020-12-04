library(ggplot2)
library(patchwork)
library(class)
library(fastcluster)

set.seed(123)
n <- 100
q <- 2
sig <- 1
X <- matrix(rnorm(n*q, sd=sig), n, q)

dat <- data.frame(X1=X[, 1], X2=X[, 2])

dat.train <- dat[1:50, ]
dat.test <- dat[51:100, ]

blankdata <- ggplot(dat) + geom_point(aes(X1, X2), cex=2) +
  xlab("") + ylab("")  +
  theme_bw() + theme(axis.ticks=element_blank(), axis.text=element_blank(),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.border=element_blank()) +
  xlim(c(-2.8, 2.8)) +  ylim(c(-2.3, 3.3))


shapedata <- ggplot(dat) + geom_point(aes(X1, X2, shape=as.factor(c(rep(1, 50), rep(2, 50)))),
                                      cex=2, alpha=1) +
  xlab("") + ylab("") + theme_bw(base_size=18) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                      legend.position="none") +
  scale_shape_manual(name="Train",
                     values=c("square", "triangle"))  +
  xlim(c(-2.8, 2.8)) +  ylim(c(-2.3, 3.3)) +
  ggtitle("(a) Data") + xlab("Feature 1") + ylab("Feature 2")

hc <- hclust(dist(X[1:50, ])^2, method="average")
dat.train$est <- as.factor(cutree(hc, 2))

traindata <- ggplot(dat.train) +
  geom_point(aes(X1, X2, col=est), cex=2, alpha=1, shape="square") +
  xlab("") + ylab("")  + scale_colour_manual(name="Clusters",
                                             values=c("#66c2a5", "#fc8d62")) +
  theme_bw(base_size=18) +theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position="none") +
  xlim(c(-2.8, 2.8)) +  ylim(c(-2.3, 3.3)) + ggtitle("(b) Training set") + xlab("Feature 1") + ylab("Feature 2")
traindata

dat.test$est <- knn(X[1:50, ], X[51:100, ], dat.train$est, k=3)
testdata <- ggplot(dat.test) + geom_point(aes(X1, X2, col=est), cex=2, alpha=1, shape="triangle") +
  xlab("") + ylab("")  + scale_colour_manual(name="Clusters",
                                             values=c("#66c2a5", "#fc8d62")) +
  theme_bw(base_size=18) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                legend.position="none")  +
  xlim(c(-2.8, 2.8)) +  ylim(c(-2.3, 3.3)) + ggtitle("(c) Test set") + xlab("Feature 1") + ylab("Feature 2")
testdata

load("../simulation-results/naive-type1-n100-q2-sig1.Rdata")
average <- ev[ev$Method == "average-split-Z-test-K-2", ]

qq <- ggplot(average[!is.na(average$stat), ][1:2000, ]) +
  geom_qq(aes(sample=pval, group=Model), size=0.5, distribution=qunif) +
  geom_abline(slope=1, intercept=0, col="black") +
  xlab("Theoretical quantiles") +
  ylab("Empirical quantiles") + theme_bw(base_size = 18) +
  ggtitle("(d) QQ plot") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

shapedata+traindata+testdata + qq + plot_layout(nrow=1)

ggsave(shapedata+traindata+testdata + qq + plot_layout(nrow=1),
       file="../figures/Figure2.pdf", width=12, height=3)

