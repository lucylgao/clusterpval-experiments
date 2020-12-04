library(ggplot2)
library(patchwork)

ev_cat <- NULL

n <- 150
for(q in c(10)) { 
for(effect in rev(seq(3, 7, by=0.5))) { 
    name_of_sim <- paste("../simulation-results/power-n", n, "-q", q, "-effect-", effect, ".Rdata", sep="")
    if(file.exists(name_of_sim)) { 
      load(name_of_sim) 
      ev_cat <- rbind(ev_cat, ev)
    } 
}
}

average <- ev_cat[ev_cat$Method == "average-iso-test-K-3", ]
centroid <- ev_cat[ev_cat$Method == "centroid-iso-test-K-3", ]
single <- ev_cat[ev_cat$Method == "single-iso-test-K-3", ]
complete <- ev_cat[ev_cat$Method == "complete-iso-test-K-3", ]

threshold <- 10
g1 <- ggplot(average[average$nmin >= threshold, ])+ 
  geom_smooth(aes(x=effect, y=as.numeric(rejects)), 
              method="gam", method.args = list(family = "binomial")) + ylim(c(0, 1)) +
  ylab(expression("Power at"~alpha~"= 0.05")) + 
  xlab(expression(paste("Effect size (", Delta, ")", sep=""))) + 
  ggtitle("(a) Average linkage") 
  
g2 <- ggplot(centroid[centroid$nmin >= threshold, ])+ 
  geom_smooth(aes(x=effect, y=as.numeric(rejects)), 
              method="gam", method.args = list(family = "binomial")) + ylim(c(0, 1)) + 
  ylab(expression("Power at"~alpha~"= 0.05")) + 
  xlab(expression(paste("Effect size (", Delta, ")", sep=""))) + 
  ggtitle("(b) Centroid linkage")

g3 <- ggplot(single[single$nmin >= threshold, ])+ 
  geom_smooth(aes(x=effect, y=as.numeric(rejects)), 
              method="gam", method.args = list(family = "binomial")) + ylim(c(0, 1)) + 
  ylab(expression("Power at"~alpha~"= 0.05")) + 
  xlab(expression(paste("Effect size (", Delta, ")", sep=""))) + 
  ggtitle("(c) Single linkage")

g4 <- ggplot(complete[complete$nmin >= threshold, ])+ 
  geom_smooth(aes(x=effect, y=as.numeric(rejects)), 
              method="gam", method.args = list(family = "binomial")) + ylim(c(0, 1)) + 
  ylab(expression("Power at"~alpha~"= 0.05")) + 
  xlab(expression(paste("Effect size (", Delta, ")", sep=""))) + 
  ggtitle("(d) Complete linkage")

ggsave(((g1 + g2 + g3 + g4) + plot_layout(nrow=1)) & theme_bw(base_size = 13), 
       file="../figures/Figure6abcd.pdf", 
       width=12, height=2.5)

p1 <- ggplot(average[average$nmin < threshold, ])+ 
  geom_smooth(aes(x=effect, y=as.numeric(rejects)), 
              method="gam", method.args = list(family = "binomial")) + ylim(c(0, 1)) + 
  ylab(expression("Power at"~alpha~"= 0.05")) + 
  xlab(expression(paste("Effect size (", Delta, ")", sep=""))) + 
  ggtitle("(e) Average linkage") 

p2 <- ggplot(centroid[centroid$nmin <  threshold, ])+ 
  geom_smooth(aes(x=effect, y=as.numeric(rejects)), 
              method="gam", method.args = list(family = "binomial")) + ylim(c(0, 1)) + 
  ylab(expression("Power at"~alpha~"= 0.05")) + 
  xlab(expression(paste("Effect size (", Delta, ")", sep=""))) + 
  ggtitle("(f) Centroid linkage") 

p3 <- ggplot(single[single$nmin < threshold, ])+ 
  geom_smooth(aes(x=effect, y=as.numeric(rejects)), 
              method="gam", method.args = list(family = "binomial")) + ylim(c(0, 1)) + 
  ylab(expression("Power at"~alpha~"= 0.05")) + 
  xlab(expression(paste("Effect size (", Delta, ")", sep=""))) + 
  ggtitle("(g) Single linkage") 

p4 <- ggplot(complete[complete$nmin < threshold, ])+ 
  geom_smooth(aes(x=effect, y=as.numeric(rejects)), 
              method="gam", method.args = list(family = "binomial")) + ylim(c(0, 1)) + 
  ylab(expression("Power at"~alpha~"= 0.05")) + 
  xlab(expression(paste("Effect size (", Delta, ")", sep=""))) + 
  ggtitle("(h) Complete linkage")

ggsave(((p1 + p2 + p3 + p4) + plot_layout(nrow=1)) & theme_bw(base_size = 13), 
       file="../figures/Figure6efgh.pdf", width=12, height=2.5)


