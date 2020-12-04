library(ggplot2)
library(patchwork)

ev_cat <- NULL

n <- 150

for(nfeat in c(2, 10, 100)) {
  for(sig in c(1, 2)) { 
    load(paste("../simulation-results/type1-n", n, "-q", nfeat, "-sig", sig, ".Rdata", sep="")) 
    ev_cat <- rbind(ev_cat, ev)
  }
}


average <- ev_cat[ev_cat$Method == "average-iso-test-K-3", ]
centroid <- ev_cat[ev_cat$Method == "centroid-iso-test-K-3", ]
single <- ev_cat[ev_cat$Method == "single-iso-test-K-3", ]
complete <- ev_cat[ev_cat$Method == "complete-iso-test-K-3", ]

#### QQ Plots #### 
p1 <- ggplot(average) + 
  geom_qq(aes(sample=pval, group=Model, colour=Model), size=0.5, distribution=qunif) + 
  geom_abline(slope=1, intercept=0, col="black") + 
  xlab("Theoretical quantiles") + 
  ylab("Empirical quantiles") + theme_bw(base_size = 15) + 
  ggtitle("(a) Average linkage")

p2 <- ggplot(centroid) + 
  geom_qq(aes(sample=pval, group=Model, colour=Model), size=0.5, distribution=qunif) + 
  geom_abline(slope=1, intercept=0, col="black") + 
  xlab("Theoretical quantiles") + 
  ylab("Empirical quantiles") + theme_bw(base_size = 15) + 
  ggtitle("(b) Centroid linkage")

p3 <- ggplot(single) + 
  geom_qq(aes(sample=pval, group=Model, colour=Model), size=0.5, distribution=qunif) + 
  geom_abline(slope=1, intercept=0, col="black") + 
  xlab("Theoretical quantiles") + 
  ylab("Empirical quantiles") + theme_bw(base_size = 15) + 
  ggtitle("(c) Single linkage")

p4 <- ggplot(complete) + 
  geom_qq(aes(sample=pval, group=Model, colour=Model), size=0.5, distribution=qunif) + 
  geom_abline(slope=1, intercept=0, col="black") + 
  xlab("Theoretical quantiles") + 
  ylab("Empirical quantiles") + theme_bw(base_size = 15) + 
  ggtitle("(d) Complete linkage")

combined <- (p1 + p2 + p3 + p4) & theme(legend.position = "bottom") & 
  scale_colour_manual(name="Model parameters", 
                      values = alpha(c("blue", "darkblue", "orange", "darkorange", "darkgreen", "darkolivegreen3"), 
                                     1),
                      labels=c(expression("q=2,"~sigma~"=1 "),
                               expression("q=2,"~sigma~"=10 "), 
                               expression("q=10,"~sigma~"=1 "), 
                               expression("q=10,"~sigma~"=10 "), 
                               expression("q=100,"~sigma~"=1 "), 
                               expression("q=100,"~sigma~"=2"))) & 
  guides(colour=guide_legend(nrow=1,byrow=TRUE, 
                             override.aes = list(size=5))) 

ggsave(combined + plot_layout(nrow=1, guides = "collect"), file="../figures/Figure4.pdf", 
       height=3.5, width=12.5)
