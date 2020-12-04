## @knitr init
library(simulator) 

if(Sys.getenv("SGE_TASK_ID") != "") {
  this.sim.id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
} else {
  this.sim.id <- 1
}

source("model_functions.R")
source("power_method_functions.R")
source("eval_functions.R")

all.settings <- expand.grid(c(150), c(10), seq(3, 7, by=0.5), c(1))
all.settings <- all.settings[order(all.settings$Var2), ]

this.setting <- all.settings[this.sim.id, ]
n <- as.numeric(this.setting[1])
nfeat <- as.numeric(this.setting[2])
len <- as.numeric(this.setting[3])
sig <- as.numeric(this.setting[4])

name_of_simulation <- paste("power-n", n, "-q", nfeat, "-effect-", len/sig, sep="")
sim <- new_simulation(name=name_of_simulation, label=name_of_simulation)

sim <- sim %>% 
  generate_model(make_three_clusters_mod, n=n, q=nfeat, len=len, sig=sig, seed=this.sim.id)
sim <- sim %>% simulate_from_model(nsim=1000, index=1:10) 
sim <- sim %>% run_method(selective_test_methods,
                          parallel=list(socket_names=10, libraries=c("clusterpval"))) 
sim <- sim %>% evaluate(c(pval, stat, rejects, effect, recover, n1, n2, nmin, boundary))
ev <- sim %>% evals %>% as.data.frame

save(ev, file=paste("../simulation-results/", name_of_simulation, ".Rdata", sep=""))
