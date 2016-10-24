ND <- pbc2[as.numeric(pbc2$id) == 2, ]


sfit.pbc15 = survfitJM(jointFit.pbc15, newdata = ND, survTimes = c(20:25))
