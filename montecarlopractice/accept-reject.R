sampsize = 100000

Y = runif(sampsize, 0 , 1)
U = runif(sampsize, 0, 1)
fY = Y*(1-y)^4
fY = Y*(1-Y)^4
fY
cutoff = fY/0.2
U<cutoff
sum(U<cutoff)
simulated  = Y[U<cutoff]

mydf  = data.frame(val = c(simulated, sampled = rbeta(n = length(simulated), shape1 = 2, shape2 = 5)), 
                   isSim=c(rep("yes", length(simulated)), rep("no", length(simulated))))

ggplot(mydf, aes(x=val)) + geom_density(aes(color=isSim)) + 
 ggtitle(paste("A",length(simulated)))
