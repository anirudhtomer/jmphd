ND <- pbc2[as.numeric(pbc2$id) ==2, ]

start = Sys.time()

sfit.pbc15 = survfitJM(jointFit.pbc3, newdata = ND, survTimes = seq(10, 15, 1), simulate = TRUE)

end = Sys.time()

print(end-start)

