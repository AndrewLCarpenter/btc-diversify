#DEoptim
library(PortfolioAnalytics)
library(DEoptim)
library(fGarch)
library(robustbase)

# Load some canned R data
data(indexes)
indexes <- indexes[,1:4]

## Set the objective function
##############################

# Step 1: Create portfolio specification object
Wcons <- portfolio.spec(assets = colnames(indexes))

# Add box constraints...
Wcons <- add.constraint(portfolio = Wcons, type='box', min = 0, max = 1)

# Full investment constraint makes it so that sum of all asset weights must equal 1
Wcons <- add.constraint(portfolio = Wcons, type = 'full_investment')

# If we violate the constraint
# Good
constrained_objective(w = rep(1/4,4) , R = indexes, portfolio = Wcons)

#This violates the constraints that we have set, so we can turn on normalize, with means that if the sum of the weights exceeds our max constraint, the
# weight vecotr is normalized by multiplying it by sum(weights)/max_sum
constrained_objective(w = rep(1/3,4) , R = indexes, portfolio = Wcons,
                      normalize = FALSE)



# Minimum C-VaR Objective Function
##################################
#Suppose now we want to find the portfolio that minimizes the 95% portfolio CVaR subject to the weight constraints listed above.
ObjSpec = add.objective(portfolio = Wcons , type="risk", name="CVaR", 
                        arguments=list(p=0.95, clean = 'boudt'), 
                        enabled=TRUE)

constrained_objective(w = rep(1/4,4) , R = indexes, portfolio = ObjSpec)
# Same as : ES(indexes[,1:4],weights = rep(1/4,4),p=0.95,clean="boudt", portfolio_method="component")


# Minimum C-VaR Objective Function WITH CONSTANT CONDITIONAL CORRELATION MODEL OF BOLLERSLEV (*1990)
####################################################################################################
ObjSpec = add.objective(portfolio = Wcons , type="risk", name="CVaR", 
                        arguments=list(p=0.95, clean = 'boudt'), 
                        enabled=TRUE, garch=TRUE)
constrained_objective(w = rep(1/4,4) , R = indexes, portfolio = ObjSpec)


# Minimum C-VaR Concentration:
##############################
ObjSpec = add.objective(portfolio = Wcons , type="risk_budget_objective", 
                        name="CVaR", arguments=list(p=0.95, clean = 'boudt'), 
                        min_concentration=TRUE, enabled=TRUE)

constrained_objective(w = rep(1/4,4) , R = indexes, portfolio = ObjSpec,trace = TRUE)
ES(indexes[,1:4],weights = rep(1/4,4),p=0.95,clean="boudt", portfolio_method="component")



# Minimum CVarR portfolio under upper 40% CVaR allocation Constraint
# The portfolio object and functions needed to obtain the minimum CVaR portfolio under an upper
# 40% CVaR allocation objective are the following:


library(PortfolioAnalytics)
library(robustbase)
# Create the portfolio specification object
ObjSpec <- portfolio.spec(assets=colnames(rets))

# Add box constraints
ObjSpec <- add.constraint(portfolio=ObjSpec, type='box', min = 0, max=1)

# Add the full investment constraint that specifies the weights must sum to 1.
ObjSpec <- add.constraint(portfolio=ObjSpec, type="weight_sum",
                              min_sum=0.99, max_sum=1.01)
# Add objective to minimize CVaR
ObjSpec <- add.objective(portfolio=ObjSpec, type="risk", name="CVaR",
                              arguments=list(p=0.95, clean="boudt"))
# Add objective for an upper 40% CVaR allocation
ObjSpec <- add.objective(portfolio=ObjSpec, type="risk_budget_objective",
                              name="CVaR", max_prisk=0.4,
                              arguments=list(p=0.95, clean="boudt"))
set.seed(1234)
out <- optimize.portfolio(R=rets, portfolio=ObjSpec,
                             optimize_method="DEoptim", search_size=2000,
                             traceDE=5, itermax=50, trace=TRUE)
print(out)
out$DEoptim_objective_results[[length(out$DEoptim_objective_results)]]

# Extract stats from the out object into a matrix
xtract <- extractStats(out)
head(xtract)

plot.new()
chart.Weights(out)
plot.new()
chart.RiskBudget(out, risk.type="pct_contrib", col="blue", pch=18)


# Min-CVaR Concentration:

# Create the portfolio specification object
ObjSpec <- portfolio.spec(assets=colnames(rets))

# Add box constraints
ObjSpec <- add.constraint(portfolio=ObjSpec, type='box', min = 0, max=1)

# Add the full investment constraint that specifies the weights must sum to 1.
ObjSpec <- add.constraint(portfolio=ObjSpec, type="weight_sum",
                             min_sum=0.99, max_sum=1.01)
# Add objective for min CVaR concentration
ObjSpec <- add.objective(portfolio=ObjSpec, type="risk_budget_objective",
                             name="CVaR", arguments=list(p=0.95, clean="boudt"),
                             min_concentration=TRUE)
set.seed(1234)
out <- optimize.portfolio(R=rets, portfolio=ObjSpec,
                             optimize_method="DEoptim", search_size=5000,
                             itermax=50, traceDE=5, trace=TRUE)

# Dynamic optimization with rebalancing::
set.seed(1234)
out <- optimize.portfolio.rebalancing(R=rets, portfolio=ObjSpec,
                                         optimize_method="DEoptim", search_size=5000,
                                         rebalance_on="months",
                                         training_period=200,
                                         traceDE=0)



library(ROI)
install.packages('ROI.plugin.quadprog')
require(ROI.plugin.glpk)
require(ROI.plugin.quadprog)
data(edhec)
R <- edhec[, 1:6]
colnames(R) <- c("CA", "CTAG", "DS", "EM", "EQMN", "ED")
funds <- colnames(R)
# Create an initial portfolio object with leverage and box constraints
init <- portfolio.spec(assets=funds)
init <- add.constraint(portfolio=init, type="leverage",
                          min_sum=0.99, max_sum=1.01)
init <- add.constraint(portfolio=init, type="box", min=0.05, max=0.65)
#Maximize mean return with ROI
#Add an objective to maximize mean return.
maxret <- add.objective(portfolio=init, type="return", name="mean")
#Run the optimization.
opt_maxret <- optimize.portfolio(R=R, portfolio=maxret,
                                    optimize_method="ROI",
                                    trace=TRUE)
print(opt_maxret)
par(mfrow=c(1,1))
plot(opt_maxret, risk.col="StdDev", return.col="mean",
     main="Maximum Return Optimization", chart.assets=TRUE,
     xlim=c(0, 0.05), ylim=c(0,0.0085))


minvar <- add.objective(portfolio=init, type="risk", name="var")
#Run the optimization. Note that although ’var’ is the risk metric, ’StdDev’ is returned as an objective measure.
opt_minvar <- optimize.portfolio(R=R, portfolio=minvar,
                                  optimize_method="ROI", trace=TRUE)
print(opt_minvar)
plot(opt_minvar, risk.col="StdDev", return.col="mean",
     main="Minimum Variance Optimization", chart.assets=TRUE,
     xlim=c(0, 0.05), ylim=c(0,0.0085))





#ok Do it NOW
meanETL <- portfolio.spec(assets=colnames(rets))
meanETL <- add.constraint(portfolio=meanETL, type='box', min = 0, max=1)
meanETL <- add.constraint(portfolio=meanETL, type="weight_sum",
                          min_sum=0.99, max_sum=1.01)
meanETL <- add.objective(portfolio = meanETL, type = 'return',name = 'mean')
meanETL <- add.objective(portfolio = meanETL, type = 'risk', name = "CVaR",
                         arguments=list(p=0.95))
meanetl.efficient.frontier(portfolio = meanETL,R = rets,n.portfolios = 100)

opt_meanETL <- optimize.portfolio(R=rets, portfolio = meanETL, trace=TRUE, search_size = 5000)

meanETL2 <- portfolio.spec(assets=colnames(rets[,-1]))
meanETL2 <- add.constraint(portfolio=meanETL2, type='box', min = 0, max=1)
meanETL2 <- add.constraint(portfolio=meanETL2, type="weight_sum",
                          min_sum=0.99, max_sum=1.01)
meanETL2 <- add.objective(portfolio = meanETL2, type = 'return',name = 'mean')
meanETL2 <- add.objective(portfolio = meanETL2, type = 'risk', name = "CVaR",
                         arguments=list(p=0.95))
meanetl.efficient.frontier(portfolio = meanETL,R = rets,n.portfolios = 100)


ef2 <- meanetl.efficient.frontier(portfolio = meanETL, R = rets,n.portfolios = 50)
plot(ef[,2], ef[,1]*sqrt(252), type='l')

ef2 <- meanetl.efficient.frontier(portfolio = meanETL2, R = rets[,-1],n.portfolios = 25)
lines(ef2[,2], ef2[,1]*sqrt(252), col='red')


########################
meanETL <- portfolio.spec(assets=colnames(rets))
meanETL <- add.constraint(portfolio=meanETL, type='box', min = 0, max=1)
meanETL <- add.constraint(portfolio=meanETL, type="weight_sum",
                          min_sum=0.99, max_sum=1.01)
meanETL <- add.objective(portfolio = meanETL, type = 'return',name = 'mean')
meanETL <- add.objective(portfolio = meanETL, type = 'risk', name = "CVaR",
                         arguments=list(p=0.95))
meanETL <- add.objective(portfolio = meanETL, type = 'risk_budget', name = "CVaR",
                         max_prisk=.3,
                         arguments=list(p=0.95))

meanetl.efficient.frontier(portfolio = meanETL,R = rets,n.portfolios = 100)

opt_meanETL <- optimize.portfolio(R=rets, portfolio = meanETL,
                                  optimize_method = "DEoptim",
                                  trace=TRUE, search_size = 5000)
print(opt_meanETL)

plot(opt_meanETL, risk.col="CVaR", return.col="mean",
     main="Risk Budget mean-CVaR Optimization",
     xlim=c(0,0.12), ylim=c(0.005,0.009))
opt_meanETL
