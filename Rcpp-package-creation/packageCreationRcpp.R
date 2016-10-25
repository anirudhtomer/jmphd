library(devtools)
library(Rcpp)
library(RcppArmadillo)

setwd(choose.dir())

RcppArmadillo.package.skeleton("JMbayesCpp")

compileAttributes(pkgdir = ".", verbose = getOption("verbose"))

setwd("..")
install("JMbayesCpp")
