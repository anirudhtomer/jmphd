a=2.7; b=6.3; c=2.669 # initial values

Nsim=5000

X=rep(runif(1),Nsim) # initialize the chain

for (i in 2:Nsim){
  Y=runif(1)
  rho=dbeta(Y,a,b)/dbeta(X[i-1],a,b)
  X[i]=X[i-1] + (Y-X[i-1])*(runif(1)<rho)
}