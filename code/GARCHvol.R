# Need to show that bitcoin volatility is decreasing:
library(rugarch)

btc <- read.csv('~/googleDrive/Research/BTC_diversification/Data/coindesk-bpi-USD-close_data-2010-07-17_2016-04-16.csv')
btc <- btc[-c(nrow(btc),nrow(btc)-1),]
btc <- as.xts(btc[,2],order.by = as.POSIXct(btc[,1],tz = 'GMT'))
colnames(btc) = c('BTC')

rets <- na.omit(ROC(btc))

# Define daily GARCH spec:
spec_d <- ugarchspec(mean.model = list(armaOrder = c(1,1)),
                     variance.model = list(model = 'sGARCH',
                                           garchOrder = c(1,1)),
                     distribution.model = 'std')


fit <- ugarchfit(spec = spec_d,data = clean.boudt(rets)[[1]])

garchVol <- xts(fit@fit$sigma, order.by = index(rets))
plot(garchVol, main='Conditional GARCH(1,1) Volatility')
lines(rollmean(garchVol,k = 252), col='red', lwd=2)






volume <- function(sec = securities){
  volume <- read.csv(
    'https://blockchain.info/charts/trade-volume?showDataPoints=false&timespan=all&show_header=true&daysAverageString=1&scale=0&format=csv&address=')
  volume <- xts(volume[,2],order.by = as.POSIXct(volume[,1], format='%d/%m/%Y', tz='GMT'))['2011-01-03/2016-05-20']
  btc <- read.csv('~/googleDrive/Research/BTC_diversification/Data/coindesk-bpi-USD-close_data-2010-07-17_2016-04-16.csv')
  btc <- btc[-c(nrow(btc),nrow(btc)-1),]
  btc <- as.xts(btc[,2],order.by = as.POSIXct(btc[,1],tz = 'GMT'))['2011-01-03/2016-5-20']
  volume <- volume*btc
  colnames(volume) = c('BTC')
  
  
  
  getSymbols(Symbols = securities,
             src = 'yahoo',
             from = '2010-01-03', to = '2016-05-20')
  df <- volume
  for (sym in securities){
    df <- cbind(df,Vo(get(sym))*Ad(get(sym)))
    df <- cbind(df,Vo(get(sym)))
  }
  colnames(df)[2:length(colnames(df))] <- securities
  
  df <- na.omit(df)
  df
}

volume <- volume()

plot(rollmean(volume[,1],k = 10))


mean(volume[,1]['2011-01-01/2012-01-01'])
mean(volume[,1]['2012-01-01/2013-01-01'])
mean(volume[,1]['2013-01-01/2014-01-01'])
mean(volume[,1]['2014-01-01/2015-01-01'])
mean(volume[,1]['2015-01-01/2016-05-20'])
mean(volume[,1]['2016-01-01/2016-05-20'])


1356898902/3353027
