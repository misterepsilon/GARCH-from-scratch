rm(list = ls())

library(PerformanceAnalytics)
library(GAS)
library(here)

# Import the functions
source(here('src','functions.R'))

# Set and Create directories
data_dir <- here::here('data')
plots_dir <- here::here('plots')
backtest_dir <- here::here('backtest')

# Load the data 
load(here(data_dir,'indices.rda'))

# Take data from January 2005
prices <- prices['2005-01-01/']

# Get the logreturns of the 2 series
logret <-  PerformanceAnalytics::CalculateReturns(prices, 'log')[-1]

# Plot the logreturns
png(file = here('plots', "logret.png"))
par(mfrow = c(2, 1))
plot(logret[, 1], 
     type = 'l', 
     col = 'black', 
     xlab = 'Time', 
     ylab = 'Log-returns',
     main = 'S&P500')
plot(logret[, 2],
     type = 'l', 
     col = 'deeppink3', 
     xlab = 'Time', 
     ylab = 'Log-returns',
     main = 'FTSE100')
mtext("S&P 500 / FTSE100 Log-returns", side = 3, line = -2, outer = TRUE)
dev.off()

# Plot the histogram of the logreturns

png(file = here('plots', "hist_logret_S&P500.png"))
hist(logret[, 1], 
     col = 'blue', 
     xlab = 'Log-returns',
     main = 'Histogram S&P500 log returns', 
     prob = TRUE,
     nclass = round(10*log(length(logret[, 1]))),
     xlim = c(min(logret[, 1]), 
              max(logret[, 1]))
     )
dev.off()

png(file = here('plots', "hist_logret_FTSE.png"))
hist(logret[, 2],
     col = 'deeppink3', 
     xlab = 'Log-returns',
     main = 'Histogram FTSE100',
     prob = TRUE,
     nclass = round(10*log(length(logret[, 2]))),
     xlim = c(min(logret[, 2]), 
              max(logret[, 2]))
     )
dev.off()

##### STATIC ESTIMATIONS : estimation of VaR at T+1 at 95% risk level for both indices

P = 1000
Risk_Level = 0.95

## s-GARCH

forcast_SP500 <- 
      f_forecast_var(as.numeric(logret[1:P, 1]), Risk_Level)
forcast_FTSE <- 
      f_forecast_var(as.numeric(logret[1:P, 2]), Risk_Level)

# predicted VaR at T+1 

VaR_SP500 <- 
      forcast_SP500$VaR_Forecast
VaR_FTSE100 <- 
      forcast_FTSE$VaR_Forecast

# MLE s-GARCH parameters

sGARCH_param_SP500 <- 
      forcast_SP500$GARCH_param
sGARCH_param_FTSE <- 
      forcast_FTSE$GARCH_param



## GJR-GARCH

forcast_SP500_gjr <- 
       f_forecast_var_gjr(as.numeric(logret[1:P, 1]), Risk_Level)
forcast_FTSE_gjr <- 
       f_forecast_var_gjr(as.numeric(logret[1:P, 2]), Risk_Level)

# predicted VaR at T+1 

VaR_SP500_gjr <- 
        forcast_SP500_gjr$VaR_Forecast
VaR_FTSE100_gjr <- 
        forcast_FTSE_gjr$VaR_Forecast


# MLE gjr-GARCH parameters

gjr_GARCH_param_SP500 <- 
        forcast_SP500_gjr$GARCH_param
gjr_GARCH_param_FTSE <- 
        forcast_FTSE_gjr$GARCH_param

##### DYNAMIC ESTIMATION
# Rolling window to get all the VaR for the 1000 obs after the first 1000

## sGARCH
roll_garch_SP500 <- 
  f_rolling_forecast_var(y = logret[1:(2*P), 1], 
                         level = Risk_Level, model = 'sGARCH', 
                         window_size = P, nb_forecast = P)

roll_garch_FTSE100 <- 
  f_rolling_forecast_var(y = logret[1:(2*P), 2], 
                         level = Risk_Level, model = 'sGARCH', 
                         window_size = P, nb_forecast = P)
## gjr-GARCH
roll_gjrgarch_SP500 <- 
  f_rolling_forecast_var(y = as.numeric(logret[1:(2*P), 1]), 
                         level = Risk_Level, model = 'gjrGARCH', 
                         window_size = P, nb_forecast = P)

roll_gjrgarch_FTSE100 <- 
  f_rolling_forecast_var(y = as.numeric(logret[1:(2*P), 2]), 
                         level = Risk_Level, model = 'gjrGARCH', 
                         window_size = P, nb_forecast = P)


# Create a dataframe with the VaR and conditional variances for the sGARCH and gjr-GARCH
VaR_CV <- xts(
  data.frame( var_SP500 = roll_garch_SP500$roll_VaR_Forecast, 
              var_FTSE100 = roll_garch_FTSE100$roll_VaR_Forecast,
              cv_SP500 = roll_garch_SP500$roll_ConditionalVariances,
              cv_FTSE100 = roll_garch_FTSE100$roll_ConditionalVariances,
              gjr_var_SP500 = roll_gjrgarch_SP500$roll_VaR_Forecast,
              gjr_var_FTSE100 = roll_gjrgarch_FTSE100$roll_VaR_Forecast,
              gjr_cv_SP500 = roll_gjrgarch_SP500$roll_ConditionalVariances,
              gjr_cv_FTSE100 = roll_gjrgarch_FTSE100$roll_ConditionalVariances),
           order.by = index(logret)[(P + 1):(2*P)])

# Create a dataframe with the VaR, returns and conditional variances for the sGARCH and gjr-GARCH
VaR_Rets_CV <-  xts(
  data.frame(logret[(P + 1):(2*P), ], 
             VaR_CV),
           order.by = index(logret)[(P + 1):(2*P)])

# Plot the sGARCH VaR series for both indices and save in a file
png(file = here('plots', "sGARCH_VaR.png"))
par(mfrow = c(3, 1))
plot(VaR_Rets_CV$var_SP500, 
     type = 'l', 
     col = 'black', 
     xlab = 'Time', 
     ylab = 'VaR',
     main = 'Predicted VaR S&P500')
plot(VaR_Rets_CV$var_FTSE100, 
     type = 'l', 
     col = 'deeppink3', 
     xlab = 'Time', 
     ylab = 'VaR',
     main = 'Predicted VaR FTSE100')
plot(VaR_Rets_CV[, 3:4], 
     xlab = 'Time', 
     ylab = 'VaR', 
     main = 'Predicted VaR S&P500 & FTSE100')
mtext("GARCH - VaR", side = 3, line = -2, outer = TRUE)
dev.off()

# Plot the gjr-GARCH VaR series for both indices and save in a file
png(file = here('plots', "gjrGARCH_VaR.png"))
par(mfrow = c(3, 1))
plot(VaR_Rets_CV$gjr_var_SP500, 
     type = 'l', 
     col = 'black', 
     xlab = 'Time', 
     ylab = 'VaR',
     main = 'VaR S&P500')
plot(VaR_Rets_CV$gjr_var_FTSE100, 
     type = 'l', 
     col = 'deeppink3', 
     xlab = 'Time', 
     ylab = 'VaR', 
     main = 'VaR FTSE100')
plot(VaR_Rets_CV[, 7:8], 
     xlab = 'Time', 
     ylab = 'VaR', 
     main = 'VaR S&P500 & FTSE100')
mtext("GJR-GARCH - VaR", side = 3, line = -2, outer = TRUE)
dev.off()

# Plot the GARCH-VAR and GJR-GARCH-VAR series for the S&P500 and save in a file
png(file = here('plots', "comp_sGARCH_gjrGARCH_VAR_SP500.png"))
par(mfrow = c(3, 1))
plot(VaR_Rets_CV$var_SP500, 
     type = 'l', 
     col = 'black', 
     xlab = 'Time', 
     ylab = 'VaR', 
     main = 'S&P500 - Predicted VaR sGARCH')
plot(VaR_Rets_CV$gjr_var_SP500, 
     type = 'l', 
     col = 'deeppink3', 
     xlab = 'Time', 
     ylab = 'VaR',
     main = 'S&P500 - Predicted VaR gjr-GARCH')
plot(VaR_Rets_CV[, c(3, 7)], 
     xlab = 'Time', 
     ylab = 'VaR', 
     main = 'S&P500 - sGARCH & gjrGARCH prediction')
mtext("S&P500 - VaR", side = 3, line = -2, outer = TRUE)
dev.off()

# Plot the GARCH-VAR and GJR-GARCH-VAR series for the FTSE100 and save in a file
png(file = here('plots', "comp_sGARCH_gjrGARCH_VAR_FTSE.png"))
par(mfrow = c(3, 1))
plot(VaR_Rets_CV$var_FTSE100, 
     type = 'l', 
     col = 'black', 
     xlab = 'Time', 
     ylab = 'VaR', 
     main = 'FTSE100 - Predicted VaR sGARCH')
plot(VaR_Rets_CV$gjr_var_FTSE100, 
     type = 'l', 
     col = 'deeppink3', 
     xlab = 'Time', 
     ylab = 'VaR', 
     main = 'FTSE100 - Predicted VaR gjr-GARCH')
plot(VaR_Rets_CV[, c(4, 8)], 
     xlab = 'Time', 
     ylab = 'VaR', 
     main = 'FTSE100 - sGARCH & gjrGARCH estimation')
mtext("FTSE100 - VaR", side = 3, line = -2, outer = TRUE)
dev.off()



# Plot the VaR and the returns for the FTSE100 in the same plot 
png(file = here('plots', "Rolling_wind_VaR_FTSE100.png"))
plot(index(VaR_Rets_CV), 
     as.numeric(VaR_Rets_CV$FTSE100), 
     type = 'l', 
     col = 'black', 
     xlab = 'Time', 
     ylab = 'Returns', 
     main = 'FTSE100 - Rolling Window VaR & Returns')
lines(index(VaR_Rets_CV), 
      as.numeric(VaR_Rets_CV$var_FTSE100), 
      col = 'blue')
lines(index(VaR_Rets_CV), 
      as.numeric(VaR_Rets_CV$gjr_var_FTSE100), 
      col = 'deeppink3')
legend("topright", 
       legend = c("Returns", "Predicted sVaR", "Predicted gjrVaR"),
       col = c("black", "blue", "deeppink3"), 
       lty = 1:1, cex = 0.8)
dev.off()


# Plot the VaR and the returns for the S&P500 in the same plot
png(file = here('plots', "Rolling_wind_VaR_SP500.png"))
plot(index(VaR_Rets_CV),
     as.numeric(VaR_Rets_CV$SP500), 
     type = 'l', 
     col = 'black', 
     xlab = 'Time', 
     ylab = 'Returns', 
     main = 'S&P500 - Rolling Window VaR & Returns')
lines(index(VaR_Rets_CV), 
      as.numeric(VaR_Rets_CV$var_SP500), 
      col = 'blue')
lines(index(VaR_Rets_CV), 
      as.numeric(VaR_Rets_CV$gjr_var_SP500), 
      col = 'deeppink3')

legend("topright", 
       legend = c("Returns", "Predicted sVaR", "Predicted gjrVaR"),
       col = c("black", "blue", "deeppink3"), 
       lty = 1:1, cex = 0.8)
dev.off()

# Coverage
# S&P500
sum(VaR_Rets_CV$SP500 < VaR_Rets_CV$var_SP500)
sum(VaR_Rets_CV$SP500 < VaR_Rets_CV$gjr_var_SP500)

# FTSE100
sum(VaR_Rets_CV$FTSE100 < VaR_Rets_CV$var_FTSE100)
sum(VaR_Rets_CV$FTSE100 < VaR_Rets_CV$gjr_var_FTSE100)


#### BONUS QUESTION 3 - BACKTESTING
# Backtesting the VaR with the GAS package
GASS_backtest_SP500 <- 
  GAS::BacktestVaR(data = as.numeric(VaR_Rets_CV$SP500),
                   VaR = as.numeric(VaR_Rets_CV$SP500),
                   alpha = .05)

GASS_backtest_FTSE100 <- 
  GAS::BacktestVaR(data = as.numeric(VaR_Rets_CV$FTSE100),
                    VaR = as.numeric(VaR_Rets_CV$var_FTSE100),
                    alpha = .05)

GASS_backtest_gjr_SP500 <- 
  GAS::BacktestVaR(data = as.numeric(VaR_Rets_CV$SP500),
                    VaR = as.numeric(VaR_Rets_CV$gjr_var_SP500),
                    alpha = .05)

GASS_backtest_gjr_FTSE100 <- 
  GAS::BacktestVaR(data = as.numeric(VaR_Rets_CV$FTSE100),
                    VaR = as.numeric(VaR_Rets_CV$gjr_var_FTSE100),
                    alpha = .05)

# Save the backtest results
save(VaR_Rets_CV, 
     GASS_backtest_SP500, 
     GASS_backtest_FTSE100,
     GASS_backtest_gjr_SP500,
     GASS_backtest_gjr_FTSE100,
     file = here('backtest', "backtest_VaR_sGARCH_gjrGARCH.rda"))

# load(file = here('backtest', "backtest_VaR_sGARCH_gjrGARCH.rda"))

#### BONUS QUESTION 4 - EQW PTF

## Compute the sGARCH innovations for both indices

# sGARCH
std_SP500 <- VaR_Rets_CV$SP500/sqrt(VaR_Rets_CV$cv_SP500)
std_FTSE100 <- VaR_Rets_CV$FTSE100/sqrt(VaR_Rets_CV$cv_FTSE100)

#Plot the sGARCH inovations for both indices

par(mfrow = c(2, 1))
plot(xts(std_SP500, order.by = index(VaR_Rets_CV)), 
     type = 'l', 
     col = 'black', 
     xlab = 'Time', 
     ylab = 'Standardised Returns', 
     main = 'Standardised Returns S&P500')
plot(xts(std_FTSE100, order.by = index(VaR_Rets_CV)), 
     type = 'l', 
     col = 'deeppink3', 
     xlab = 'Time', 
     ylab = 'Standardised Returns', 
     main = 'Standardised Returns FTSE100')
mtext("sGARCH", side = 3, line = -2, outer = TRUE)


## Compute the gjr-GARCH innovations for both indices
# gjrGARCH
std_gjrSP500 <- VaR_Rets_CV$SP500/sqrt(VaR_Rets_CV$gjr_cv_SP500)
std_gjrFTSE100 <- VaR_Rets_CV$FTSE100/sqrt(VaR_Rets_CV$gjr_cv_FTSE100)

#Plot the gjr-GARCH inovations for both indices

par(mfrow = c(2, 1))
plot(xts(std_gjrSP500, order.by = index(VaR_Rets_CV)),
     type = 'l', col = 'black',
     xlab = 'Time',
     ylab = 'Standardised Returns',
     main = 'Standardised Returns S&P500')
plot(xts(std_gjrFTSE100, order.by = index(VaR_Rets_CV)),
     type = 'l',
     col = 'deeppink3',
     xlab = 'Time',
     ylab = 'Standardised Returns',
     main = 'Standardised Returns FTSE100')
mtext("gjrGARCH", side = 3, line = -2, outer = TRUE)


# Simulate the bivariate gaussian distribution for sGARCH and gjr-GARCH
garch_sim <- f_sim_biv_gauss('sGARCH')
gjr_sim <- f_sim_biv_gauss('gjrGARCH')

# Generate returns from the simulated bivariate gaussian distribution for sGARCH and gjr-GARCH
VaR_Rets_CV$sGARCH_sim_SP500 <- garch_sim[, 1]*sqrt(VaR_Rets_CV$cv_SP500)
VaR_Rets_CV$sGARCH_sim_FTSE100 <- garch_sim[, 2]*sqrt(VaR_Rets_CV$cv_FTSE100)

VaR_Rets_CV$gjrGARCH_sim_SP500 <- gjr_sim[, 1]*sqrt(VaR_Rets_CV$gjr_cv_SP500)
VaR_Rets_CV$gjrGARCH_sim_FTSE100 <- gjr_sim[, 2]*sqrt(VaR_Rets_CV$gjr_cv_FTSE100)

# Scatter plot of the sGARCH simulated returns

plot(as.numeric(VaR_Rets_CV$sGARCH_sim_SP500), 
     as.numeric(VaR_Rets_CV$sGARCH_sim_FTSE100), 
     main = 'sGARCH Simulated Returns', 
     xlab = 'S&P500', 
     ylab = 'FTSE100')


# Scatter plot of the gjr-GARCH simulated returns

plot(as.numeric(VaR_Rets_CV$gjrGARCH_sim_SP500), 
     as.numeric(VaR_Rets_CV$gjrGARCH_sim_FTSE100), 
     main = 'gjr-GARCH Simulated Returns', 
     xlab = 'S&P500', 
     ylab = 'FTSE100')


# Create an equal weight portfolio based on sGARCH simulated returns
VaR_Rets_CV$sGARCH_eqw_portfolio_rets <- 
  (VaR_Rets_CV$sGARCH_sim_SP500 + VaR_Rets_CV$sGARCH_sim_FTSE100)/2

# Create an equal weight portfolio based on gjr-GARCH simulated returns
VaR_Rets_CV$gjrGARCH_eqw_portfolio_rets <- 
  (VaR_Rets_CV$gjrGARCH_sim_SP500 + VaR_Rets_CV$gjrGARCH_sim_FTSE100)/2


# Get 95*% VaR for the equal weight portfolio based on sGARCH simulated returns
sGARCH_eqw_portfolio_VaR <- f_forecast_var(y = VaR_Rets_CV$sGARCH_eqw_portfolio_rets,
                                           level = Risk_Level)$VaR_Forecast

# Get 95*% VaR for the equal weight portfolio based on gjr-GARCH simulated returns
gjrGARCH_eqw_portfolio_VaR <- f_forecast_var_gjr(y = VaR_Rets_CV$gjrGARCH_eqw_portfolio_rets,
                                             level = Risk_Level)$VaR_Forecast


# Get Static VaR for each 

var_SP500_ <- f_forecast_var(y = VaR_Rets_CV$SP500,
                                    level = Risk_Level)$VaR_Forecast
var_FTSE100_ <- f_forecast_var(y = VaR_Rets_CV$FTSE100,
                 level = Risk_Level)$VaR_Forecast

var_gjr_SP500_ <- f_forecast_var_gjr(y = VaR_Rets_CV$SP500,
                             level = Risk_Level)$VaR_Forecast
var_gjr_FTSE100_ <- f_forecast_var_gjr(y = VaR_Rets_CV$FTSE100,
                               level = Risk_Level)$VaR_Forecast

# Plot the histogram of the sGARCH equal weight portfolio returns
png(here('plots', 'sGARCH_hist_eqw_portfolio_rets_&_VaR.png'))
hist(as.numeric(VaR_Rets_CV$sGARCH_eqw_portfolio_rets), 
     col = 'blue', 
     xlab = 'Returns',
     main = 'sGARCH Equal Weight Portfolio',
     prob = TRUE,
     nclass = round(10*log(length(VaR_Rets_CV))))
abline(v = sGARCH_eqw_portfolio_VaR, col = 'deeppink3')
text(
  x = (sGARCH_eqw_portfolio_VaR - 0.012),
  y = 40,
  labels = "VaR(95%)",
  col = "deeppink3",
  srt = 360)
dev.off()   

# Plot the histogram of the gjr-GARCH equal weight portfolio returns
png(here('plots', 'gjrGARCH_hist_eqw_portfolio_rets_&_VaR.png'))
hist(as.numeric(VaR_Rets_CV$gjrGARCH_eqw_portfolio_rets), 
     col = 'blue', 
     xlab = 'Returns',
     main = 'gjr-GARCH Equal Weight Portfolio',
     prob = TRUE,
     nclass = round(10*log(length(VaR_Rets_CV))))
abline(v = gjrGARCH_eqw_portfolio_VaR, col = 'deeppink3')
text(
  x = (gjrGARCH_eqw_portfolio_VaR - 0.011),
  y = 50,
  labels = "VaR(95%)",
  col = "deeppink3",
  srt = 360)
dev.off()  






