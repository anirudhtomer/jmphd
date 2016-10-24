library(devtools)
library(Rcpp)
library(RcppArmadillo)

setwd("C:/Users/838035/Dropbox/PhD/src/testing/")
setwd("D:/Dropbox/PhD/src/testing/")

RcppArmadillo.package.skeleton("JMbayesCpp")

setwd("C:/Users/838035/Dropbox/PhD/src/testing/JMbayesCpp")
setwd("D:/Dropbox/PhD/src/testing/JMbayesCpp")

compileAttributes(pkgdir = ".", verbose = getOption("verbose"))

setwd("..")
install("JMbayesCpp")
