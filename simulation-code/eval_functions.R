## @knitr metrics

stat <- new_metric(name="stat", label="stat", metric=function(model, out) out$results$stat)
pval <- new_metric(name="pval", label="pval", metric=function(model, out) out$results$pval)
rejects <- new_metric(name="rejects", label="rejects", metric=function(model, out) out$results$pval <= 0.05)
effect <- new_metric(name="effect", label="effect", metric=function(model, out) out$effect)
recover <- new_metric(name="recover", label="recover", metric=function(model, out) out$recover)
n1 <- new_metric(name="n1", label="n1", metric=function(model, out) out$npair[1])
n2 <- new_metric(name="n2", label="n2", metric=function(model, out) out$npair[2])
nmin <- new_metric(name="nmin", label="nmin", metric=function(model, out) min(out$npair))
boundary <- new_metric(name="boundary", label="boundary", 
                       metric=function(model, out) { 
                         
                         if(!is.null(out$results$trunc)) { 
                           row <- which(apply(out$results$trunc, 1, 
                                              function(x) (out$results$stat >= x[1] & out$results$stat <= x[2])))
                           out$results$stat - out$results$trunc[row, 1]
                         } else { 
                           NA 
                         }
                       } 
)