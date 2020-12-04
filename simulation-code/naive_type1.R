## @knitr init
library(simulator) 

n <- 100
nfeat <- 2
sig <- 1

source("model_functions.R")
source("type1_method_functions.R")
source("eval_functions.R")

name_of_simulation <- paste("naive-type1-n", n, "-q", nfeat, "-sig", sig,  sep="")
sim <- new_simulation(name=name_of_simulation, label=name_of_simulation)

sim <- sim %>% 
  generate_model(make_null_mod, n=n, q=nfeat, sig=sig)
sim <- sim %>% simulate_from_model(nsim=400, index=1:10) 
sim <- sim %>% run_method(c(multivariate_Z_test_methods, 
                            sample_split_methods, 
                            selective_methods), 
                          parallel=list(socket_names=10, libraries=c("clusterpval"))) 
sim <- sim %>% evaluate(c(stat, pval, rejects))
ev <- sim %>% evals %>% as.data.frame

save(ev, file=paste("../simulation-results/", name_of_simulation, ".Rdata", sep=""))
