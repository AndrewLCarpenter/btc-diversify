# Size penalty
securities <- c('SPY',       # S&P 500
                'IWM',       # Russell 2000 (small cap)
                'EFA',       # Foreign developed market equities
                'VNQ',       # Real Estate
                'GSG',       # Commodities
                'BND',       # Vangaurd Bond Index
                'HYG')       # High Yield Bond

importData <- function(sec = securities, returns = TRUE){
  # Import Bitcoin price data
  btc <- read.csv('~/googleDrive/Research/BTC_diversification/Data/coindesk-bpi-USD-close_data-2010-07-17_2016-04-16.csv')
  btc <- btc[-c(nrow(btc),nrow(btc)-1),]
  btc <- as.xts(btc[,2],order.by = as.POSIXct(btc[,1],tz = 'GMT'))['2011-01-03::']
  colnames(btc) = c('BTC')
  
  getSymbols(Symbols = securities,
             src = 'yahoo',
             from = '2011-01-03', to = '2016-04-16')
  df <- btc
  for (sym in securities){
    df <- cbind(df,Ad(get(sym)))
  }
  colnames(df)[2:length(colnames(df))] <- securities
  
  if(returns==TRUE){
    df <- na.omit(ROC(na.omit(df)))
  }
  head(df)
}

rets <- importData()




volume <- function(sec = securities){
  volume <- read.csv(
    'https://blockchain.info/charts/trade-volume?showDataPoints=false&timespan=all&show_header=true&daysAverageString=1&scale=0&format=csv&address=')
  volume <- xts(volume[,2],order.by = as.POSIXct(volume[,1], format='%d/%m/%Y', tz='GMT'))['2011-01-03/2016-05-20']
  ##btc <- read.csv('~/googleDrive/Research/BTC_diversification/Data/coindesk-bpi-USD-close_data-2010-07-17_2016-04-16.csv')
  #btc <- btc[-c(nrow(btc),nrow(btc)-1),]
  #btc <- as.xts(btc[,2],order.by = as.POSIXct(btc[,1],tz = 'GMT'))['2011-01-03/2016-5-20']
  #volume <- volume*btc
  colnames(volume) = c('BTC')
  
  
  
  getSymbols(Symbols = securities,
             src = 'yahoo',
             from = '2011-01-03', to = '2016-05-20')
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

plot(volume[,1])
mean(volume[,2])/mean(volume[,3])

colMeans(rets)[1]/100


plot(volume[,3])
lines(volume[,1], col='red')
