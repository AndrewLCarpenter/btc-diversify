library(DEoptim)
library(ROI)
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



meanETL <- add.objective(portfolio=init, type="return", name="mean")
meanETL <- add.objective(portfolio=meanETL, type="risk", name="ETL",
                            arguments=list(p=0.95))
#Run the optimization. The default random portfolio method is ’sample’.
opt_meanETL <- optimize.portfolio(R=R, portfolio=meanETL,
                                     optimize_method="random",
                                     trace=TRUE, search_size=2000)
print(opt_meanETL)
stats_meanETL <- extractStats(opt_meanETL)
dim(stats_meanETL)
plot(opt_meanETL, risk.col="ETL", return.col="mean",
      main="mean-ETL Optimization", neighbors=20)

#Calculate and plot the portfolio component ETL contribution.
pct_contrib <- ES(R=R, p=0.95, portfolio_method="component",
                     weights=extractWeights(opt_meanETL))
barplot(pct_contrib$pct_contrib_MES, cex.names=0.8, las=3, col="lightblue")


# change the box constraints to long only
init$constraints[[2]]$min <- rep(0, 6)
init$constraints[[2]]$max <- rep(1, 6)
rb_meanETL <- add.objective(portfolio=init, type="return", name="mean")
rb_meanETL <- add.objective(portfolio=rb_meanETL, type="risk", name="ETL",
                               arguments=list(p=0.95))
rb_meanETL <- add.objective(portfolio=rb_meanETL, type="risk_budget",
                               name="ETL", max_prisk=0.4, arguments=list(p=0.95))
#Run the optimization. Set traceDE=5 so that every fifth iteration is printed. The default is to
#print every iteration.
opt_rb_meanETL <- optimize.portfolio(R=R, portfolio=rb_meanETL,
                                        optimize_method="DEoptim",
                                        search_size=2000,
                                        trace=TRUE, traceDE=5)
print(opt_rb_meanETL)


plot(opt_rb_meanETL, risk.col="ETL", return.col="mean",
      main="Risk Budget mean-ETL Optimization",
      xlim=c(0,0.12), ylim=c(0.005,0.009))
plot.new()
chart.RiskBudget(opt_rb_meanETL, risk.type="percentage", neighbors=25)


#Maximize mean return per unit ETL with ETL equal contribution
#to risk

eq_meanETL <- add.objective(portfolio=init, type="return", name="mean")
eq_meanETL <- add.objective(portfolio=eq_meanETL, type="risk", name="ETL",
                               arguments=list(p=0.95))
eq_meanETL <- add.objective(portfolio=eq_meanETL, type="risk_budget",
                               name="ETL", min_concentration=TRUE,
                               arguments=list(p=0.95))
#Run the optimization. Set traceDE=5 so that every fifth iteration is printed. The default is to
#print every iteration.
opt_eq_meanETL <- optimize.portfolio(R=R, portfolio=eq_meanETL,
                                        optimize_method="DEoptim",
                                        search_size=2000,
                                        trace=TRUE, traceDE=5)
print(opt_eq_meanETL)


plot.new()
plot(opt_eq_meanETL, risk.col="ETL", return.col="mean",
       main="Risk Budget mean-ETL Optimization",
       xlim=c(0,0.12), ylim=c(0.005,0.009))
#Chart the contribution to risk in percentage terms. It is clear in this chart that the optimization
#results in a near equal risk contribution portfolio.
plot.new()
chart.RiskBudget(opt_eq_meanETL, risk.type="percentage", neighbors=25)



#Combine the optimizations for easy comparison.
opt_combine <- combine.optimizations(list(meanETL=opt_meanETL,
                                             rbmeanETL=opt_rb_meanETL,
                                             eqmeanETL=opt_eq_meanETL))
# View the weights and objective measures of each optimization
extractWeights(opt_combine)

obj_combine <- extractObjectiveMeasures(opt_combine)
chart.Weights(opt_combine, plot.type="bar", legend.loc="topleft", ylim=c(0, 1))

plot.new()
chart.RiskReward(opt_combine, risk.col="ETL", return.col="mean",
                    main="ETL Optimization Comparison", xlim=c(0.018, 0.024),
                    ylim=c(0.005, 0.008))

#Calculate the STARR of each optimization
STARR <- obj_combine[, "mean"] / obj_combine[, "ETL"]
barplot(STARR, col="blue", cex.names=0.8, cex.axis=0.8,
           las=3, main="STARR", ylim=c(0,1))
plot.new()
chart.RiskBudget(opt_combine, match.col="ETL", risk.type="percent",
                    ylim=c(0,1), legend.loc="topright")
