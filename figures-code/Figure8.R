library(palmerpenguins)
library(ggplot2)
library(patchwork)
library(fastcluster)

options(ggplot2.discrete.colour=list(RColorBrewer::brewer.pal(6, "Accent")))

# Subset to female penguins & bill/flipper length variables
penguins <- penguins[complete.cases(penguins), ]
dat <- penguins[penguins$sex == "female", c(1, 3, 5)]

X <- dat[, -c(1)]
X <- as.matrix(X)


# Cluster and visualize data
K <- 6
hc <- hclust(dist(X)^2, method="average")
hc$labels <- dat$species
table( hc$labels, cutree(hc, K))

p1 <- ggdendro::ggdendrogram(hc, labels=F) + ggtitle("(a) Dendrogram") + 
  geom_hline(yintercept=hc$height[nrow(X)-5] - 20, col="red", linetype="dashed") + 
  ylab("Height") + theme(axis.text.x=element_blank(), 
                         plot.title = element_text(size=22), 
                         axis.text.y=element_text(size=22)) 

p2 <- ggplot(data = penguins[penguins$sex == "female", ]) + 
  geom_point(aes(x=bill_length_mm, y = flipper_length_mm, 
                 fill = as.factor(cutree(hc, 6)),
                 shape=as.factor(species)), size = 5, colour="black") + 
  scale_fill_discrete(name="Clusters", 
                        guide=guide_legend(ncol=2, 
             override.aes=list(shape=21))) + 
  scale_shape_manual(name="Species", values=c(21, 24, 22), 
                     guide=guide_legend(override.aes=list(fill="black"))) + 
  xlab("Bill length (mm)") + ylab("Flipper length (mm)") + 
  coord_flip() +  
  theme_bw(base_size=22) + ggtitle("(b) Data") + theme(legend.position="right")

ggsave(p1+p2, file="../figures/Figure8.pdf", width=14, height=4)
