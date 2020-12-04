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

threshold <- 10
p1 <- ggplot(average[average$effect == 5 & average$nmin >= threshold, ]) + 
  geom_smooth(aes(x=boundary, y=log10(pval))) + 
  ylab(expression(log[10]~"(p-value)")) + 
  xlab("Distance to left boundary") + ggtitle(expression(paste("(a) min {|", C[k], "|, |", C["k'"], "|} \u2265 10")))

p2 <- ggplot(average[average$effect == 5 & average$nmin < threshold, ]) + 
  geom_smooth(aes(x=boundary, y=log10(pval))) + 
  ylab(expression(log[10]~"(p-value)")) + 
  xlab("Distance to left boundary") + ggtitle(expression(paste("(b) min {|", C[k], "|, |", C["k'"], "|} < 10")))

ggsave((p1+p2) & theme_bw(base_size=25), 
       file="../figures/Figure7.pdf", width=16, height=4.5, device=cairo_pdf)
