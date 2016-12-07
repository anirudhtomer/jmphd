knots = seq(1,100, by = 10)
splinemat = data.frame(bs(1:100, knots = knots))
splinemat$X = 1:100

splinemat_long=reshape(splinemat, direction='long', idvar='X',
                   varying=list(c(1:4)),
                   v.names=c('XSpline'))
splinemat_long$time = as.factor(splinemat_long$time)

ggplot(data=splinemat_long) + geom_line(aes(x = X, y=XSpline, group=time, color=time)) + 
    scale_x_continuous(breaks = round(seq(1, 100, by = 10),0)) + xlab("X (knots as labels)")
