## @knitr init
library(simulator) 

if(Sys.getenv("SGE_TASK_ID") != "") { 
  this.sim.id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
} else { 
  this.sim.id <- 1
}

n <- 200
nfeat <- 10
sig <- 1

source("model_functions.R")
source("type1_est_method_functions.R")
source("eval_functions.R")

name_of_simulation <- paste("type1-est-n", n, "-q", nfeat, "-pt", this.sim.id, sep="")
sim <- new_simulation(name=name_of_simulation, label=name_of_simulation)

sim <- sim %>% 
  generate_model(make_two_clusters_mod_with_id, n=n, q=nfeat, 
                 len=as.list(seq(2, 6, by=2)), id=this.sim.id, 
                 sig=sig, seed=this.sim.id, vary_along = "len")
sim <- sim %>% simulate_from_model(nsim=5000, index=1:10) 
sim <- sim %>% run_method(c(selective_methods_est), 
                          parallel=list(socket_names=10, libraries=c("clusterpval"))) 
sim <- sim %>%  evaluate(c(pval, stat, rejects, effect, n1, n2, nmin, nmax, sd))
ev <- sim %>% evals %>% as.data.frame

save(ev, file=paste("../simulation-results/", name_of_simulation, ".Rdata", sep=""))
