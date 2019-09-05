#! /usr/bin/env Rscript
d<-scan("stdin", quiet=TRUE)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
cat("Median ",median(d),"Geometric Mean ",gm_mean(d), sep="\n")
