library(ggplot2)
library(patchwork)
library(dplyr)
ev_cat <- NULL

n <- 200
nfeat <- 10
sig <- 1

for(id in 1:4) { 
  name_of_sim <- paste("../simulation-results/type1-est-n", n, "-q", nfeat, "-pt", id, ".Rdata", sep="")
  
  if(file.exists(name_of_sim)) { 
    load(name_of_sim) 
    ev_cat <- rbind(ev_cat, ev)
  } 
}

ev_cat$len <- rep(NA, nrow(ev_cat))
ev_cat$len[grep("len_2", ev_cat$Model)] <- 2
ev_cat$len[grep("len_4", ev_cat$Model)] <- 4
ev_cat$len[grep("len_6", ev_cat$Model)] <- 6
ev_cat$len <- as.factor(ev_cat$len)

average <- ev_cat[ev_cat$Method == "average-iso-est-test-K-3" & ev_cat$effect == 0, ]
centroid <- ev_cat[ev_cat$Method == "centroid-iso-est-test-K-3" & ev_cat$effect == 0, ]
single <- ev_cat[ev_cat$Method == "single-iso-est-test-K-3" & ev_cat$effect == 0, ]
complete <- ev_cat[ev_cat$Method == "complete-iso-est-test-K-3" & ev_cat$effect == 0, ]

average_filter <- as.data.frame(average %>% group_by(len) %>% slice_head(n=500))
centroid_filter <- as.data.frame(centroid %>% group_by(len) %>% slice_head(n=500))
single_filter <- as.data.frame(single %>% group_by(len) %>% slice_head(n=500))
complete_filter <- as.data.frame(complete %>% group_by(len) %>% slice_head(n=500))

#### QQ Plots #### 
p1 <- ggplot(average_filter) + 
  geom_qq(aes(sample=pval, group=len, colour=len), size=0.5, distribution=qunif) + 
  geom_abline(slope=1, intercept=0, col="black") + 
  xlab("Uniform(0, 1) quantiles") + 
  ylab("Empirical quantiles") + theme_bw(base_size = 15) + 
  ggtitle("(a) Average linkage")

p2 <- ggplot(centroid_filter) + 
  geom_qq(aes(sample=pval, group=len, colour=len), size=0.5, distribution=qunif) + 
  geom_abline(slope=1, intercept=0, col="black") + 
  xlab("Uniform(0, 1) quantiles") + 
  ylab("Empirical quantiles") + theme_bw(base_size = 15) + 
  ggtitle("(b) Centroid linkage")

p3 <- ggplot(single_filter) + 
  geom_qq(aes(sample=pval, group=len, colour=len), size=0.5, distribution=qunif) + 
  geom_abline(slope=1, intercept=0, col="black") + 
  xlab("Uniform(0, 1) quantiles") + 
  ylab("Empirical quantiles") + theme_bw(base_size = 15) + 
  ggtitle("(c) Single linkage")

p4 <- ggplot(complete_filter) + 
  geom_qq(aes(sample=pval, group=len, colour=len), size=0.5, distribution=qunif) + 
  geom_abline(slope=1, intercept=0, col="black") + 
  xlab("Uniform(0, 1) quantiles") + 
  ylab("Empirical quantiles") + theme_bw(base_size = 15) + 
  ggtitle("(d) Complete linkage")

combined <- (p1 + p2 + p3 + p4) & theme(legend.position = "bottom") & 
  scale_colour_manual(name=expression("Distance between clusters ("~delta~")"), 
                      values = c("#9ecae1", "#4292c6", "#08519c")) & 
  guides(colour=guide_legend(nrow=1,byrow=TRUE, 
                             override.aes = list(size=5))) 

combined

ggsave(combined + plot_layout(nrow=1, guides = "collect"), file="../figures/FigureS1.pdf", 
              height=3.5, width=12.5)
