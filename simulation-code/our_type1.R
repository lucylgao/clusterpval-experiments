## @knitr init
library(simulator) 

if(Sys.getenv("SGE_TASK_ID") != "") {
  this.sim.id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
} else {
  this.sim.id <- 1
}

all.settings <- expand.grid(c(c(150)), c(2, 10, 100), c(1, 2))

this.setting <- all.settings[this.sim.id, ]
n <- as.numeric(this.setting[1])
nfeat <- as.numeric(this.setting[2])
sig <- as.numeric(this.setting[3])

source("model_functions.R")
source("type1_method_functions.R")
source("eval_functions.R")

name_of_simulation <- paste("type1-n", n, "-q", nfeat, "-sig", sig,  sep="")
sim <- new_simulation(name=name_of_simulation, label=name_of_simulation)

sim <- sim %>% 
  generate_model(make_null_mod, n=n, q=nfeat, sig=sig, seed=this.sim.id)
sim <- sim %>% simulate_from_model(nsim=200, index=1:10) 
sim <- sim %>% run_method(c(selective_methods), 
                          parallel=list(socket_names=10, libraries=c("clusterpval"))) 
sim <- sim %>% evaluate(c(stat, pval, rejects))
ev <- sim %>% evals %>% as.data.frame

save(ev, file=paste("../simulation-results/", name_of_simulation, ".Rdata", sep=""))
