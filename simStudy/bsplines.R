data = read.table("MontrealTemp1.txt", header = T)


############
# x + x^2 + x^3 + (x-knot1)^3 + (x-knot2)^3 + (x-knotk)^3
#where if x-knot is negative then the cubic term is assumed zero. i.e. truncated values are used
############
D = 3
K = 5
knots = nrow(data) * (1:K ) / (K+1)
X1 = outer(data$x ,1:D , "^")
X2 = outer(data$x , knots , ">") * outer ( data$x , knots , "-")^D
X = cbind(X1, X2)

fitlm = lm(y~X, data = data)

data$yfit1 = fitlm$fitted.values

ggplot(data = data) + geom_line(aes(x=data$x, y=data$yfit1)) + 
  geom_point(aes(x=data$x, y=data$y, colour="red"))

######### B splines ######
X = bs (data$x , knots = knots,
         degree=D , intercept = TRUE )

fitlm_bs = lm(y~X, data = data)

data$yfit_bs = fitlm_bs$fitted.values

ggplot(data = data) + geom_line(aes(x=data$x, y=data$yfit_bs)) + 
  geom_point(aes(x=data$x, y=data$y, colour="red"))
