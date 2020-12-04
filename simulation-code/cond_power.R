## @knitr init
library(simulator) 

if(Sys.getenv("SGE_TASK_ID") != "") { 
  this.sim.id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
} else { 
  this.sim.id <- 1
}

source("model_functions.R")
source("cond_power_method_functions.R")
source("eval_functions.R")

n <- 30
nfeat <- 10
sig <- 1

name_of_simulation <- paste("cond-n", n, "-q", nfeat, "-pt", this.sim.id, sep="")
sim <- new_simulation(name=name_of_simulation, label=name_of_simulation)

sim <- sim %>% 
  generate_model(make_three_clusters_mod_with_id, n=n, q=nfeat, sig=sig, id=this.sim.id, 
                 seed=this.sim.id, len = as.list(seq(4, 7, by=0.5)), vary_along = "len")
sim <- sim %>% simulate_from_model(nsim=10000, index=1:10) 
sim <- sim %>% run_method(selective_test_methods,
                          parallel=list(socket_names=10, libraries=c("clusterpval"))) 
sim <- sim %>% evaluate(c(pval, stat, rejects, effect, recover, n1, n2, nmin))
ev <- sim %>% evals %>% as.data.frame

save(ev, file=paste("../simulation-results/", name_of_simulation, ".Rdata", sep=""))
