as.ggsList <- function(jmfit){
  if(inherits(jmfit, "JMbayes")|inherits(jmfit, "mvJMbayes")){
    lapply(jmfit$mcmc[[1]], function(resMatrix){
      ggs(as.mcmc(resMatrix))
    })
  }else{
    stop("The object you passed is not a JMbayes object")
  }
}