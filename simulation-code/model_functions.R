generate_null_data <- function(n, q, sig) { 
  dat <- matrix(rnorm(n*q, sd=sig), n, q) 
  
  return(list(data=dat, means=matrix(0, 1, q), clusters=rep(1, n)))
}


generate_three_clusters_data <- function(n, q, len, sig, id) { 
  cl <- rep_len(1:3, length.out=n)
  
  a <- len/2
  mu <- rbind(c(-a, rep(0, q-1)), c(rep(0, q-1), sqrt(3)*a), c(a, rep(0, q-1)))
  
  dat <- matrix(rnorm(n*q, sd=sig), n, q) + mu[cl, ]
  
  return(list(data=dat, means=mu, clusters=cl))
}

make_null_mod <- function(n, q, sig) { 
  new_model(name="null-mod", 
            label=sprintf("Global null data (n=%s, q=%s, sig=%s)", n, q,  sig), 
            params=list(n=n, q=q,  sig=sig), 
            simulate=function(nsim) { 
              lapply(1:nsim, function(x) generate_null_data(n, q, sig))
            })  
}

make_three_clusters_mod <- function(n, q, len, sig) { 
  new_model(name="three-clusters-mod", 
            label=sprintf("Three clusters data (n=%s, q=%s, len=%s, sig=%s)", n, q, len, sig), 
            params=list(n=n, sig=sig), 
            simulate=function(nsim) { 
              lapply(1:nsim, function(x) generate_three_clusters_data(n, q, len, sig))
            })  
}

make_three_clusters_mod_with_id <- function(n, q, len, sig, id) { 
  new_model(name="three-clusters-mod-id", 
            label=sprintf("Three clusters data (n=%s, q=%s, len=%s, sig=%s, id=%s)", n, q, len, sig, id), 
            params=list(n=n, sig=sig, id=id), 
            simulate=function(nsim) { 
              lapply(1:nsim, function(x) generate_three_clusters_data(n, q, len, sig, id))
            })  
}

