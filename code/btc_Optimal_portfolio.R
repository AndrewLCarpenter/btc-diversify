library(PerformanceAnalytics)
library(quadprog)
library(quantmod)
library(ggplot2)


securities <- c('SPY',       # S&P 500
                'IWM',       # Russel 2000 (small cap)
                'EFA',       # Foreign developed market equities
                'EEM',       # Foreign emerging market equities
                'BND')       # Vangaurd Bond)

importData <- function(dir = '~/googleDrive/Research/BTC_diversification/Data/coindesk-bpi-USD-close_data-2010-07-17_2016-03-01.csv', sec = securities){
  # Import Bitcoin price data
  btc <- read.csv(dir)
  btc <- btc[-c(nrow(btc),nrow(btc)-1),]
  btc <- as.xts(btc[,2],order.by = as.POSIXct(btc[,1],tz = 'UTC'))['2011-01-03::']
  colnames(btc) = c('XBP')
  getSymbols(Symbols = securities,
             src = 'yahoo',
             from = '2011-01-03', to = Sys.Date())
  df <- btc
  for (sym in securities){
    df <- cbind(df,Ad(get(sym)))
  }
  df
}

df <- importData()

colnames(df)[2:length(colnames(df))] <- securities
rets <- na.omit(ROC(na.omit(df)))

head(rets)

# Correlation Table
corr <- cor(rets, method = 'pearson')
corr[upper.tri(corr)] <- ' '
corr

rets2 <- rets[,-1]
head(rets2)


efront <- function(returns, nport = 100, max.weight = NULL){
  # solving the QP
  Dmat <- cov(returns)
  n <- nrow(Dmat)
  
  # Constraints (short selling prohibited)
  Amat <- cbind(1, diag(n))
  bvec <- c(1, rep(0, n))
  meq <- 1
  port <- 1
  
  # If max weight is set:
  if(!is.null(max.weight)){
    if(max.weight > 1 | max.weight <0){
      stop("max.weight must be greater than 0 and less than 1")
    }
    if(max.weight * n < 1){
      stop("max.weights must add to 1")
    }
    Amat <- cbind(Amat, -1*diag(n))
    bvec <- c(bvec, rep(-max.weight, n))
  }
  
  
  # Initialize matrix for allocation and statistics
  eff <- matrix(nrow=nport, ncol=n+3)
  colnames(eff) <- c(colnames(returns), "Std.Dev", "Exp.Return", "sharpe")
  
  # Fill matrix
  for (i in seq(from = 0,to = .6,length.out = nport)){
    dvec <- colMeans(returns)*i
    qp <- solve.QP(Dmat, dvec, Amat, bvec, meq)
    qp$solution[abs(qp$solution) <= 1e-7] <- 0
    
    
    # load results into data frame
    eff[port,1:n] <- qp$solution
    eff[port,"Std.Dev"] <- sqrt(qp$solution %*% Dmat %*% qp$solution * 252)   # weights^T * covarMatrix * weights
    eff[port,"Exp.Return"] <- as.numeric(colMeans(returns) %*% qp$solution * 252)            # returns * weights
    eff[port,"sharpe"] <- eff[port,"Exp.Return"] / eff[port,"Std.Dev"]
    
    port <- port+1
    
  }
  as.data.frame(eff)
}

eff1 <- efront(returns=rets)
eff2 <- efront(returns = rets2)
plot(eff1[c('Std.Dev','Exp.Return')], xlim = c(0, .20), ylim = c(0, .25), type = 'l', lty = 7, 
     xlab = 'risk(sd)', ylab = 'return', main = 'Efficient Frontier', col = 'navyblue', lwd=2)
lines(eff2[c('Std.Dev','Exp.Return')], type = 'l', lty = 7, col = 'red', lwd=2)
lines(seq(0,1,length.out = nrow(eff1)),max(eff1['sharpe'])*seq(0,1,length.out = nrow(eff1)), type = 'l', lty=2)
lines(seq(0,1,length.out = nrow(eff1)),max(eff2['sharpe'])*seq(0,1,length.out = nrow(eff1)), type = 'l', lty=2)
grid(lwd = 1.5)
points(sd.annualized(rets, scale = 252), colMeans(rets*252))
text(sd.annualized(rets, scale = 252), colMeans(rets*252), c('Bitcoin',securities), cex = .7, adj = -.3)
