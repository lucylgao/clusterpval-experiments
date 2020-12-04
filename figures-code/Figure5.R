library(ggplot2)
library(patchwork)

ev_cat <- NULL

for(id in 1:3) { 
  name_of_sim <- paste("../simulation-results/cond-n", 30, "-q", 10, "-pt", id, ".Rdata", sep="")
  if(file.exists(name_of_sim)) { 
      load(name_of_sim) 
      ev_cat <- rbind(ev_cat, ev)
  } 
}

collate <- NULL
for(effect in seq(4, 7, by=0.5)) { 
  ev <- ev_cat[ev_cat$Model %in% unique(ev_cat$Model)[grep(paste("len_", effect, "/", sep=""), 
                                      unique(ev_cat$Model))], ]
  
  df <- data.frame(effect=rep(effect, 4), 
                   Method=c("average-test-K-3", "centroid-test-K-3", "complete-test-K-3", 
                            "single-test-K-3"), 
                   rejects=tapply(ev[ev$recover == 1, ]$rejects, ev[ev$recover == 1, ]$Method, mean), 
                   nsim=tapply(ev[ev$recover == 1, ]$rejects, ev[ev$recover == 1, ]$Method, length), 
                   total=rep(nrow(ev)/4, 4))
  collate <- rbind(collate, df)
}

collate$detect <- collate$nsim/collate$total

plot1 <- ggplot(collate) + geom_line(aes(x=effect, y=rejects, group=Method, colour=Method)) + 
  geom_errorbar(aes(group=Method, colour=Method, 
                    x=effect,
                    y=rejects, ymin=rejects - qnorm(0.975)*sqrt(rejects*(1 - rejects)/nsim),
                    ymax=rejects + qnorm(0.975)*sqrt(rejects*(1 - rejects)/nsim), width=0.05))+ 
  theme_bw(base_size=20) + scale_colour_discrete(name="Linkage", 
                                                 labels=c("Average", "Centroid", "Complete", "Single")) + 
  xlab(expression(paste("Distance between clusters (", delta, ")", sep=""))) + 
  ylab("Conditional power") + ggtitle("(a) Conditional power") 

plot2 <- ggplot(collate) + geom_line(aes(x=effect, y=detect, group=Method, colour=Method)) + 
  geom_errorbar(aes(group=Method, colour=Method, 
                    x=effect,
                    y=detect, ymin=detect - qnorm(0.975)*sqrt(detect*(1 - detect)/total),
                    ymax=detect + qnorm(0.975)*sqrt(detect*(1 - detect)/total), width=0.05)) + 
  theme_bw(base_size=20) + scale_colour_discrete(name="Linkage", 
                                                 labels=c("Average", "Centroid", "Complete", "Single")) + 
xlab(expression(paste("Distance between clusters (", delta, ")", sep=""))) + 
  ylab("Detection probability") + ggtitle("(b) Detection probability") + ylim(c(0, 1))

combined <- (plot1 + plot2) & 
  theme(legend.position = "right") & guides(colour=guide_legend(nrow=4, byrow=TRUE) )

ggsave(combined + plot_layout(guides = 'collect'), file="../figures/Figure5.pdf", width=14, height=4)
