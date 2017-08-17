#Paper Codessss
library(PerformanceAnalytics)
library(quadprog)
library(quantmod)
library(xtable)
library(Matrix)
getData <- function(startDate='2012-01-01', path='~/googleDrive/Research/BTC_diversification/Data/logReturns.csv'){
  rets <- read.csv(path)
  rets <- xts(rets[,-1], order.by = as.POSIXct(rets[,1]))[paste(startDate,'/2016-05-20',sep='')]
  rets
}
rets <- getData()
btc <- read.csv('~/googleDrive/Research/BTC_diversification/Data/coindesk-bpi-USD-close_data-2010-07-17_2016-05-23.csv')
btc <- btc[-c(nrow(btc),nrow(btc)-1),]
btc <- as.xts(btc[,2],order.by = as.POSIXct(btc[,1],tz = 'UTC'))['2011-01-01/2016-05-20']
colnames(btc) = c('BTC/USD')
# CAPM stuff
# Assume risk-free rate of 0 for now
mean_returns <- colMeans(rets)
sd_returns <- apply(rets, 2, sd)
rbind(mean_returns, sd_returns)


CAPM <- function(rets, market){
  mkt_return <- mean(market)
  mkt_variance <- var(market)
  df <- matrix(nrow = 3, ncol = ncol(rets))
  
  for (i in 1:ncol(rets)){
    reg <- lm(rets[,i] ~ rets$SPY)
    beta <- reg$coefficients[2]
    exp_return <- beta * mkt_return 
    variance <- beta * as.numeric(mkt_variance) + var(reg$residuals)
    df[1,i] <- exp_return
    df[2,i] <- variance
    df[3,i] <- beta
  }
  df <- as.data.frame(df)
  colnames(df) <- colnames(rets)
  
  covar <- matrix(nrow=ncol(rets), ncol=ncol(rets))
  # Now make covar matrix
  for (i in 1:ncol(rets)){
    for (j in 1:ncol(rets)){
      covar[i,j] <- ifelse(i==j,yes =  df[2,i], df[3,i]*df[3,j]*mkt_variance)
    }
  }
  colnames(covar) <- colnames(df)
  rownames(covar) <- colnames(df)
  
  out <- list(ExpRet = df[1,], 
              covar = covar)
  return(out)
  
}


CAPM_efront <- function(returns, nport = 100, max_concentration = .35){
  # Do the regressions
  capm <- CAPM(returns, returns$SPY)
  Dmat <- as.matrix(nearPD(capm$covar)$mat)
  n <- nrow(Dmat)
  
  if (1/n > max_concentration)
    stop("Maximum asset concentration is too small.",call. = FALSE)
  
  # Constraints: sum(w_i) = 1 & w_i >= 0 & w_i <= max_concentration
  Amat <- cbind(1, 
                diag(n), 
                -1*diag(n))
  bvec <- c(1, 
            rep(0, n), 
            rep(-1*max_concentration, n))
  meq <- 1
  
  # Initialize matrix for allocation and statistics
  eff <- matrix(nrow=nport, ncol=n+3)
  colnames(eff) <- c(colnames(returns), "Std.Dev", "Exp.Return", "Sharpe")
  port <- 1 #(counter variable)
  
  # Fill matrix
  for (i in seq(from = -.5,to = 1,length.out = nport)){
    dvec <- capm$ExpRet*i
    qp <- solve.QP(Dmat, dvec, Amat, bvec, meq)
    qp$solution[abs(qp$solution) <= 1e-7] <- 0
    
    # load results into data frame
    eff[port,1:n] <- qp$solution
    eff[port,"Std.Dev"] <- sqrt(qp$solution %*% Dmat %*% qp$solution * 252)   # weights^T * covarMatrix * weights
    eff[port,"Exp.Return"] <- as.numeric(as.matrix(capm$ExpRet) %*% qp$solution * 252)            # returns * weights
    eff[port,"Sharpe"] <- eff[port,"Exp.Return"] / eff[port,"Std.Dev"]
    port <- port+1
  }
  as.data.frame(eff)
}

q <- CAPM_efront(rets, max_concentration = .35, nport =1000)
q2 <- CAPM_efront(rets[,-1], max_concentration = .35, nport =1000)
plot(q[c('Std.Dev','Exp.Return')], xlim = c(0,.2), ylim = c(-.1,.25), type = 'l', lty = 7, 
     xlab = 'Standard Deviation', ylab = 'Expected Return', main = 'Efficient Frontier', col = 'navyblue', lwd=2)
lines(q2[c('Std.Dev','Exp.Return')], type = 'l', lty = 3, col = 'red', lwd=2)
grid(lwd = 1.5)

max(q$Sharpe)
max(q2$Sharpe)


hybrid_efront <- function(returns, nport = 100, max_concentration = .35, BTC=TRUE){
  # Do the regressions
  if(BTC==FALSE){
    returns <- returns[,-1]
    capm <- CAPM(returns, returns$SPY)
    Dmat <- as.matrix(nearPD(capm$covar)$mat)
    Er <- capm$ExpRet
  }else{
    capm <- CAPM(returns, returns$SPY)
    Dmat <- cov(returns)
    Er <- capm$ExpRet
    Er[1,1] <- mean(returns[,1])
  }
  
  n <- nrow(Dmat)
  if (1/n > max_concentration)
    stop("Maximum asset concentration is too small.",call. = FALSE)
  
  # Constraints: sum(w_i) = 1 & w_i >= 0 & w_i <= max_concentration
  Amat <- cbind(1, 
                diag(n), 
                -1*diag(n))
  bvec <- c(1, 
            rep(0, n), 
            rep(-1*max_concentration, n))
  meq <- 1
  
  # Initialize matrix for allocation and statistics
  eff <- matrix(nrow=nport, ncol=n+3)
  colnames(eff) <- c(colnames(returns), "Std.Dev", "Exp.Return", "Sharpe")
  port <- 1 #(counter variable)
  
  # Fill matrix
  for (i in seq(from = -.25,to = 1,length.out = nport)){
    dvec <- Er*i
    qp <- solve.QP(Dmat, dvec, Amat, bvec, meq)
    qp$solution[abs(qp$solution) <= 1e-7] <- 0
    
    # load results into data frame
    eff[port,1:n] <- qp$solution
    eff[port,"Std.Dev"] <- sqrt(qp$solution %*% Dmat %*% qp$solution * 252)   # weights^T * covarMatrix * weights
    eff[port,"Exp.Return"] <- as.numeric(as.matrix(Er) %*% qp$solution * 252)            # returns * weights
    eff[port,"Sharpe"] <- eff[port,"Exp.Return"] / eff[port,"Std.Dev"]
    port <- port+1
  }
  as.data.frame(eff)
}


q <- hybrid_efront(rets, max_concentration = .75, nport =1000)
q2 <- hybrid_efront(rets, max_concentration = .75, nport =1000, BTC=FALSE)
plot(q[c('Std.Dev','Exp.Return')], xlim = c(0,.2), ylim = c(-.1,.25), type = 'l', lty = 7, 
     xlab = 'Standard Deviation', ylab = 'Expected Return', main = 'Efficient Frontier', col = 'navyblue', lwd=2)
lines(q2[c('Std.Dev','Exp.Return')], type = 'l', lty = 3, col = 'red', lwd=2)

z <- CAPM_efront(rets, max_concentration = 1, nport = 1000)
z2 <- CAPM_efront(rets[,-1], max_concentration = 1, nport = 1000)
lines(z[c('Std.Dev','Exp.Return')], type = 'l', lty = 3, col = 'blue', lwd=2)
lines(z2[c('Std.Dev','Exp.Return')], type = 'l', lty = 4, col = 'red', lwd=2)

grid(lwd = 1.5)

max(q$Sharpe)
max(q2$Sharpe)
head(q)
q2


summary(lm(rets[,1] ~ rets[,2]))

















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


library(quantmod)
source("http://www.quantmod.com/examples/chartSeries3d/chartSeries3d.alpha.R")
getSymbols("DGS1MO;DGS3MO;DGS6MO;DGS1;DGS2;DGS3;DGS5;DGS7;DGS10;DGS20;DGS30",
           src="FRED")
TR <- merge(DGS1MO,DGS3MO,DGS6MO,DGS1,DGS2,DGS3,DGS5,
            DGS7,DGS10,DGS20,DGS30, all=FALSE)
colnames(TR) <- c("1mo","3mo","6mo","1yr","2yr","3yr","5yr",
                  "7yr","10yr","20yr","30yr")
TR <- na.locf(TR)

chartSeries3d0(TR['2011'])
head(TR)






CAPM_plotFront <- function(rets, max_weight, BTC_weight, xlim = c(0,.3), ylim = c(-.3,.3), main = 'Efficient Frontier'){
  q <- CAPM_efront(rets,nport = 10000, max_concentration = max_weight, BTC_percentage = BTC_weight)
  q2 <- CAPM_efront(rets[,-1],nport = 10000, max_concentration = max_weight, BTC_percentage = BTC_weight)
  plot(q[c('Std.Dev','Exp.Return')], xlim = xlim, ylim = ylim, type = 'l', lty = 7, 
       xlab = 'Standard Deviation', ylab = 'Expected Return', main = main, col = 'navyblue', lwd=2)
  lines(q2[c('Std.Dev','Exp.Return')], type = 'l', lty = 7, col = 'red', lwd=2)
  grid(lwd = 1.5)
  points(sqrt(abs(diag(CAPM(rets, rets[,2])[[2]]))*252), CAPM(rets, rets[,2])[[1]] * 252)
  text(sqrt(abs(diag(CAPM(rets, rets[,2])[[2]]))*252), CAPM(rets, rets[,2])[[1]] * 252, c('Bitcoin',securities), cex = .7, adj = -.3)
  lines(seq(0,1,length.out = nrow(q)),max(q['Sharpe'])*seq(0,1,length.out = nrow(q)), type = 'l', lty=2)
  lines(seq(0,1,length.out = nrow(q2)),max(q2['Sharpe'])*seq(0,1,length.out = nrow(q2)), type = 'l', lty=2)
}

CAPM_plotFront(rets['2011::'], max_weight = .35, BTC_weight = 0, xlim = c(0,1))

CAPM_plotFront(rets['2012::'], max_weight = .35, BTC_weight = 0, xlim = c(0,1))

CAPM_plotFront(rets['2013::'], max_weight = .35, BTC_weight = 0, xlim = c(0,1))

CAPM_plotFront(rets['2014::'], max_weight = .35, BTC_weight = 0, xlim = c(0,1))

CAPM_plotFront(rets['2015::'], max_weight = .35, BTC_weight = 0, xlim = c(0,1))

CAPM_plotFront(rets['2016::'], max_weight = .35, BTC_weight = 0, xlim = c(0,1))









long_noHeavy1 <- backtest(rets, type = 'markowitz', max_weight = .3)
long_noHeavy2 <- backtest(rets[,-1], type = 'markowitz', max_weight = .3)
r <- table.AnnualizedReturns(merge(sp,long_noHeavy2,long_noHeavy1))


eqWeighted1 <- backtest(rets, type = 'equal')
eqWeighted2 <- backtest(rets[,-1], type = 'equal')
s <- table.AnnualizedReturns(merge(sp,eqWeighted2,eqWeighted1))















