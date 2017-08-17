# Stuff to add

#Daily return vs SPY
head(rets)
plot(rets[,"XBP"], xlab = "Time",
     main = "Daily Return Comparison", major.ticks= "years",
     minor.ticks = FALSE, col = "black")
lines(rets[,'SPY'], col='blue')
legend(x = 'topright', legend = c('BTC','SPY'),
       lty = 1, col = c('black','blue'))


#Risk return scatter
chart.RiskReturnScatter(rets)

# Show kurtosis
par(mfrow=c(2,1))
chart.Histogram(rets[,1], 
                main = "Bitcoin Return Distribution with Normal Overlay", 
                breaks = 200, show.outliers = FALSE, 
                methods = c("add.density", "add.normal", "add.rug"))


chart.Histogram(rets[,2], 
                main = "SPY Return Distribution with Normal Overlay", 
                breaks = 200, show.outliers = FALSE, 
                methods = c("add.density", "add.normal", "add.rug"))

chart.Correlation(rets[,c(1,2)],histogram = TRUE, pch="+")







# MAYBE USE THIS PLOT IF YOU CAN GET IT TO LOOK NICE
chart.Scatter(eff1$Std.Dev,y = eff1$Exp.Return,xlim = c(-.15, .15), ylim = c(0, .5), type = 'l', lty = 7, 
              xlab = 'risk(sd)', ylab = 'return', main = 'Efficient Frontier', col = 'navyblue')






# Random
chart.RollingCorrelation(rets[,1], rets[,2:6], width = 24*7, colorset = rich10equal, legend.loc = 'bottomleft', main = '6 Month Rolling Correlations')

ggplot() + 
  geom_line(data = eff1, aes(x=Std.Dev, y=Exp.Return), color = 'blue')+
  geom_line(data = eff2, aes(x=Std.Dev, y=Exp.Return), color = 'red') + 
  coord_cartesian(xlim = c(0,0.2), ylim=c(0, 0.3))




plot(eff1[c('Std.Dev','Exp.Return')],  xlim = c(0, 1), ylim = c(-.03, 1),type = 'l', lty = 7, 
     xlab = 'Standard Deviation', ylab = 'Expected Return', main = 'Efficient Frontier', col = 'navyblue', lwd=2)
grid(lwd = 1.5)
lines(eff2[c('Std.Dev','Exp.Return')], type = 'l', lty = 7, col = 'red', lwd=2)
lines(seq(0,1,length.out = nrow(eff1)),max(eff1['sharpe'])*seq(0,1,length.out = nrow(eff1)), type = 'l', lty=2)
lines(seq(0,1,length.out = nrow(eff1)),max(eff2['sharpe'])*seq(0,1,length.out = nrow(eff1)), type = 'l', lty=2)
points(sd.annualized(rets), colMeans(rets*252))
text(sd.annualized(rets), colMeans(rets*252), c('Bitcoin',securities), cex = .7, adj = -.3)
