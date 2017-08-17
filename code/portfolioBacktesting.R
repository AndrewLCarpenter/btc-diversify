# Backtesting time
library(quantmod)
library(PerformanceAnalytics)
library(quadprog)

# DATA
########
securities <- c('SPY',       # S&P 500
                'IWM',       # Russell 2000 (small cap)
                'EFA',       # Foreign developed market equities
                'VNQ',       # Real Estate
                'GSG',       # Commodities
                'BND')       # Vangaurd Bond Index

importData <- function(sec = securities, returns = TRUE){
  # Import Bitcoin price data
  btc <- read.csv('~/googleDrive/Research/BTC_diversification/Data/coindesk-bpi-USD-close_data-2010-07-17_2016-04-16.csv')
  btc <- btc[-c(nrow(btc),nrow(btc)-1),]
  btc <- as.xts(btc[,2],order.by = as.POSIXct(btc[,1],tz = 'UTC'))['2011-01-03::']
  colnames(btc) = c('BTC')
  getSymbols(Symbols = securities,
             src = 'yahoo',
             from = '2011-01-03', to = Sys.Date())
  df <- btc
  for (sym in securities){
    df <- cbind(df,Ad(get(sym)))
  }
  colnames(df)[2:length(colnames(df))] <- securities
  if(returns==TRUE){
    df <- na.omit(ROC(na.omit(df)))
  }
  df
}
rets <- importData()


# Calculate Returns and remove NAs

################################################################################
# Now write a backtesting function
#####################################
# Efficient frontier function that returns optimal Portfolios
efront <- function(returns, nport = 100, max_concentration = .20, BTC_percentage = 0){
  # If bitcoin percentage is specified:
  if(BTC_percentage > 0){
    returns <- returns[,-1]
    Dmat <- cov(returns)
    n <- nrow(Dmat)
    if (1/n > max_concentration)
      stop("Maximum asset concentration is too small.",call. = FALSE)
    
    # Constraints
    Amat <- cbind(1,  
                  diag(nrow(Dmat)),
                  -1*diag(nrow(Dmat))) 
    bvec <- c(1-BTC_percentage, 
              rep(0, nrow(Dmat)), 
              rep(-1*max_concentration, nrow(Dmat)))
    meq <- 1
  }
  
  else{
    Dmat <- cov(returns)
    n <- nrow(Dmat)
    if (1/n > max_concentration)
      stop("Maximum asset concentration is too small.",call. = FALSE)
    
    # Constraints: sum(w_i) = 1 & w_i >= 0 & w_i <= max_concentration
    Amat <- cbind(1, 
                  diag(nrow(Dmat)), 
                  -1*diag(nrow(Dmat)))
    bvec <- c(1, 
              rep(0, nrow(Dmat)), 
              rep(-1*max_concentration, nrow(Dmat)))
    meq <- 1
  }
  
  # Initialize matrix for allocation and statistics
  eff <- matrix(nrow=nport, ncol=n+3)
  colnames(eff) <- c(colnames(returns), "Std.Dev", "Exp.Return", "Sharpe")
  port <- 1 #(counter variable)
  
  # Fill matrix
  for (i in seq(from = -1,to = 1,length.out = nport)){
    dvec <- colMeans(returns)*i
    qp <- solve.QP(Dmat, dvec, Amat, bvec, meq)
    qp$solution[abs(qp$solution) <= 1e-7] <- 0
    
    # load results into data frame
    eff[port,1:n] <- qp$solution
    eff[port,"Std.Dev"] <- sqrt(qp$solution %*% Dmat %*% qp$solution * 252)   # weights^T * covarMatrix * weights
    eff[port,"Exp.Return"] <- as.numeric(colMeans(returns) %*% qp$solution * 252)            # returns * weights
    eff[port,"Sharpe"] <- eff[port,"Exp.Return"] / eff[port,"Std.Dev"]
    
    qp$solution %*% Dmat 
    port <- port+1
  }
  as.data.frame(eff)
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

backtest <- function(rets, type = "equal", max_weight = .3, BTC_weight = .05, plot=FALSE){
  # If equal then do it
  if(type == "equal"){
    n <- ncol(rets)
    
    #Create weight matrix
    weights <- xts(rep.row(rep(1/n, ncol(rets)),nrow(rets)), order.by = index(rets))
    
    #Performance
    bt <- Return.portfolio(rets, weights <- weights, rebalance_on = 'monthly')
    if(plot==TRUE){
      out <- xts(exp(cumsum(bt)),order.by = index(bt))
      return(out)
      }
    }
  
  if(type == "markowitz"){
    # Find optimal portfolio using mean-variance
    optimal.portfolios <- efront(rets, max_concentration = max_weight)
    maxSharpe <- optimal.portfolios[optimal.portfolios$Sharpe == max(optimal.portfolios$Sharpe),]
    if(nrow(maxSharpe) > 1){maxSharpe <- maxSharpe[1,]}
    
    # Create weight matrix
    weights <- xts(rep.row(as.numeric(maxSharpe[,1:ncol(rets)]),nrow(rets)), order.by = index(rets))
    
    # Now we have portfoio weights, so evolve it forward in time
    bt <- Return.portfolio(rets, weights = weights, rebalance_on = 'monthly')
    if(plot==TRUE){
      out <- xts(exp(cumsum(bt)),order.by = index(bt))
      return(out)
    }
  }
  
  if(type == 'BTC'){
    optimal.portfolios <- efront(rets, max_concentration = max_weight,BTC_percentage = BTC_weight)
    maxSharpe <- optimal.portfolios[optimal.portfolios$Sharpe == max(optimal.portfolios$Sharpe),]
    if(nrow(maxSharpe) > 1){maxSharpe <- maxSharpe[1,]}
    last <- ncol(optimal.portfolios)-3
    maxSharpe <- cbind(.05, maxSharpe[,1:last])
    
    # Create weight matrix
    weights <- xts(rep.row(as.numeric(maxSharpe[,1:ncol(rets)]),nrow(rets)), order.by = index(rets))
    
    # Now we have portfoio weights, so evolve it forward in time
    bt <- Return.portfolio(rets, weights = weights, rebalance_on = 'monthly')
    if(plot==TRUE){
      out <- xts(exp(cumsum(bt)),order.by = index(bt))
      return(out)
    }
    }
  return(bt)
}

plotFront <- function(rets, max_weight, BTC_weight, xlim = c(0,.3), ylim = c(-.3,.3), main = 'Efficient Frontier'){
  q <- efront(rets,nport = 10000, max_concentration = max_weight, BTC_percentage = BTC_weight)
  q2 <- efront(rets[,-1],nport = 10000, max_concentration = max_weight, BTC_percentage = BTC_weight)
  plot(q[c('Std.Dev','Exp.Return')], xlim = xlim, ylim = ylim, type = 'l', lty = 7, 
       xlab = 'Standard Deviation', ylab = 'Expected Return', main = main, col = 'navyblue', lwd=2)
  lines(q2[c('Std.Dev','Exp.Return')], type = 'l', lty = 7, col = 'red', lwd=2)
  grid(lwd = 1.5)
  points(sd.annualized(rets, scale = 252), colMeans(rets*252))
  text(sd.annualized(rets, scale = 252), colMeans(rets*252), c('Bitcoin',securities), cex = .7, adj = -.3)
  lines(seq(0,1,length.out = nrow(q)),max(q['Sharpe'])*seq(0,1,length.out = nrow(q)), type = 'l', lty=2)
  lines(seq(0,1,length.out = nrow(q2)),max(q2['Sharpe'])*seq(0,1,length.out = nrow(q2)), type = 'l', lty=2)
}
sp <- rets[,2]

long1 <- backtest(rets, type = 'markowitz', max_weight = 1)
long2 <- backtest(rets[,-1], type = 'markowitz', max_weight = 1)
q <- table.AnnualizedReturns(merge(sp,long2,long1))
#plotFront(rets, max_weight = 1, BTC_weight = 0)

long_noHeavy1 <- backtest(rets, type = 'markowitz', max_weight = .3)
long_noHeavy2 <- backtest(rets[,-1], type = 'markowitz', max_weight = .3)
r <- table.AnnualizedReturns(merge(sp,long_noHeavy2,long_noHeavy1))
#plotFront(rets, max_weight = .3, BTC_weight = 0, ylim = c(-.25,.25))

eqWeighted1 <- backtest(rets, type = 'equal')
eqWeighted2 <- backtest(rets[,-1], type = 'equal')
s <- table.AnnualizedReturns(merge(sp,eqWeighted2,eqWeighted1))

btcFive1 <- backtest(rets, type = 'BTC', max_weight = .3, BTC_weight = .05)
btcFive2 <- backtest(rets[,-1], type = 'BTC', max_weight = .3, BTC_weight = .05)
t <- table.AnnualizedReturns(merge(sp,btcFive2, btcFive1))

annRets <- cbind(c(0,0,0),t(q[1,]),t(r[1,]),t(s[1,]),t(t[1,]))
annRets[,1] <- c('S&P 500', 'Market Portfolio', 'Market Portfolio with Bitcoin')
colnames(annRets) <- c('Annualized Return','Portfolio 1','Portfolio 2','Portfolio 3','Portfolio 4')
rownames(annRets) <- NULL
annRets















rets2 <- rets[,-1]
q <- efront(rets,nport = 100,max_concentration = 1)
q2 <- efront(rets2,nport = 10000,max_concentration = 1)
plot(q[c('Std.Dev','Exp.Return')], xlim = c(0, .3), ylim = c(-.3, .3), type = 'l', lty = 7, 
     xlab = 'risk(sd)', ylab = 'return', main = 'Efficient Frontier', col = 'navyblue', lwd=2)
lines(q2[c('Std.Dev','Exp.Return')], type = 'l', lty = 7, col = 'red', lwd=2)
grid(lwd = 1.5)
points(sd.annualized(rets, scale = 252), colMeans(rets*252))
text(sd.annualized(rets, scale = 252), colMeans(rets*252), c('Bitcoin',securities), cex = .7, adj = -.3)
lines(seq(0,1,length.out = nrow(q)),max(q['Sharpe'])*seq(0,1,length.out = nrow(q)), type = 'l', lty=2)
lines(seq(0,1,length.out = nrow(q2)),max(q2['Sharpe'])*seq(0,1,length.out = nrow(q2)), type = 'l', lty=2)
max(q['Sharpe'])
max(q2['Sharpe'])





















rows <- make.names(as.character(round(q$Sharpe, digits = 4)), unique = TRUE)
for (i in 1:length(rows)){
  rownames(x)[i] <- substring(rows[i],3,6)

}

chart.StackedBar(x,
                 colorset =  rich8equal,
                 space = 0)
axis(3, at=1:100, labels=letters[1:100])

library(reshape2)
library(ggplot2)
x <- t(q[,c(1:7)])
colnames(x) <- round(q$Exp.Return,3)
head(x)
barplot(x, col = rainbow(ncol(rets)),legend.text = names(q[,1:7]))

legend("topright", legend = legendtext, bty = "n", 
       cex = 0.7, fill = col)
)

