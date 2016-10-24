#target is beta(2,5)

#uniform(0,1) is instrumental

sampsize = 1000000

#theta = runif(sampsize, 0 , 1)
theta = rbeta(sampsize, 5 , 10)
fx = theta*(1-theta)^4
gx = dbeta(theta, 5, 10)

oldweight = fx/gx
weights = oldweight/sum(oldweight)

multinomres = c(rmultinom(1, sampsize/4, weights))

samples = rep(theta, multinomres)
plot(density(samples))