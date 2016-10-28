ND <- pbc2[as.numeric(pbc2$id) < 100 , ]

start = Sys.time()

sfit.pbc15 = survfitJM(jointFit.pbc15, newdata = ND, survTimes = seq(11, 15, 1), simulate = TRUE, M = 500)

end = Sys.time()

print(end-start)