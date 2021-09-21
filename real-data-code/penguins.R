library(palmerpenguins)
library(ggplot2)
library(clusterpval)
library(fastcluster)
library(patchwork)

source("util.R")

options(ggplot2.discrete.colour=list(RColorBrewer::brewer.pal(6, "Accent")))

# Subset to female penguins in years 2008/2009 & bill/flipper length variables
penguins <- penguins[complete.cases(penguins), ]
dat <- penguins[penguins$sex == "female" & penguins$year != 2009, c(1, 3, 5)]

X <- dat[, -c(1)]
X <- as.matrix(X)

# Split off a part to estimate the variance
X.ss <- penguins[penguins$sex == "female" & penguins$year == 2009, c(3, 5)]
sig <- sqrt(sum(scale(X.ss, scale = FALSE)^2)/(nrow(X.ss) * ncol(X.ss) -ncol(X.ss)))

# Cluster and visualize data; this is Figure 6
K <- 5
hc <- hclust(dist(X)^2, method="average")
hc$labels <- dat$species
table( hc$labels, cutree(hc, K))

p1 <- ggdendro::ggdendrogram(hc, labels=F) + ggtitle("(a) Dendrogram") + 
  geom_hline(yintercept=hc$height[nrow(X)-4] - 20, col="red", linetype="dashed") + 
  ylab("Height") + theme(axis.text.x=element_blank(), 
                         plot.title = element_text(size=22), 
                         axis.text.y=element_text(size=22)) 

p2 <- ggplot(data = penguins[penguins$sex == "female" & penguins$year != 2009, ]) + 
  geom_point(aes(x=bill_length_mm, y = flipper_length_mm, 
                 fill = as.factor(cutree(hc, 5)),
                 shape=as.factor(species)), size = 5, colour="black") + 
  scale_fill_discrete(name="Clusters", 
                      guide=guide_legend(ncol=2, 
                                         override.aes=list(shape=21))) + 
  scale_shape_manual(name="Species", values=c(21, 24, 22), 
                     guide=guide_legend(override.aes=list(fill="black"))) + 
  xlab("Bill length (mm)") + ylab("Flipper length (mm)") + 
  coord_flip() +  
  theme_bw(base_size=22) + ggtitle("(b) Data") + theme(legend.position="right")

p1+p2

# Run tests
test_hier_clusters_exact(as.matrix(X), "average", K=5, k1=1, k2=2, hcl=hc, sig=sig)
test_hier_clusters_exact(as.matrix(X), "average", K=5, k1=1, k2=3, hcl=hc, sig=sig)
test_hier_clusters_exact(as.matrix(X), "average", K=5, k1=1, k2=4, hcl=hc, sig=sig)
test_hier_clusters_exact(as.matrix(X), "average", K=5, k1=2, k2=3, hcl=hc, sig=sig)
test_hier_clusters_exact(as.matrix(X), "average", K=5, k1=2, k2=4, hcl=hc, sig=sig)
test_hier_clusters_exact(as.matrix(X), "average", K=5, k1=3, k2=4, hcl=hc, sig=sig)

wald_test(as.matrix(X), "average", K=5, k1=1, k2=2, sig=sig)
wald_test(as.matrix(X), "average", K=5, k1=1, k2=3, sig=sig)
wald_test(as.matrix(X), "average", K=5, k1=1, k2=4, sig=sig)
wald_test(as.matrix(X), "average", K=5, k1=2, k2=3, sig=sig)
wald_test(as.matrix(X), "average", K=5, k1=2, k2=4, sig=sig)
wald_test(as.matrix(X), "average", K=5, k1=3, k2=4, sig=sig)