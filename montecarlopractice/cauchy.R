xm=rcauchy(500)

f=function(y){
  -sum(log(1+(x-y)^2))
}

for (i in 1:500){
  x=xm[1:i]
  mi[i]=optimise(f,interval=c(-10,10),maximum = TRUE)$max
}
