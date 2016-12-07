ND <- pbc2[as.numeric(pbc2$id) < 2 , ]

start = Sys.time()

sfit.pbc15 = survfitJM(jointFit.pbc1, newdata = ND, survTimes = seq(11, 15, 1), simulate = TRUE, M = 500)

end = Sys.time()

print(end-start)