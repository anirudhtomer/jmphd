library(ggplot2)

numObs = 100
len = 100
cutoff = 30
origEstimator = vector("numeric", len)
newEstimator = vector("numeric", len)
for(i in 1:len){
  sample = rexp(n = numObs, rate = 0.05)
  origEstimator[i] = mean(sample)
  
  truncSample = sample[sample<cutoff]
  r = length(truncSample)
  newEstimator[i] = (sum(truncSample) + (numObs-r)*cutoff)/r
}

ggplot(data = data.frame(mean=c(origEstimator, newEstimator), type=c(rep("Orig", len), rep("New", len)))) + 
  geom_density(aes(x=mean, color=type))
  