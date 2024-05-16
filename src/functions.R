f_forecast_var <- function(y, level) {
  
  ### Compute the VaR forecast of a GARCH(1,1) 
  ### model with Normal errors at the desired risk level
  #  INPUTS
  #   y     : [vector] (T x 1) of observations (log-returns)
  #   level : [scalar] risk level (e.g. 0.95 for a VaR at the 95# risk level)
  #  OUTPUTS
  #   VaR   : [scalar] VaR forecast 
  #   sig2  : [vector] (T+1 x 1) conditional variances
  #   theta : [vector] GARCH parameters
  #  NOTE
  #   o the estimation is done by maximum likelihood
  
  # Fit a GARCH(1,1) model with Normal errors
  
  # Starting values and bounds
  theta0 <- c(0.1* var(y), 0.1, 0.8)
  
  # Get the lower bound
  LB     <- c(0,0,0)
  
  # Stationarity conditions
  
  A      <- matrix(c(0, -1, -1), nrow = 1) # Coefficients for alpha and beta in the constraint
  b      <- -1  # Right-hand side of the inequality
  
  
  # Run the optimization using constrOptim
  fit <- constrOptim(theta0,
                     f  = f_nll,
                     grad = NULL,
                     ui = matrix(c(1,0,0, 0,1,0 , 0,0,1, A), 
                                 nrow = 4, 
                                 ncol = 3, 
                                 byrow = TRUE),
                     ci = matrix(c(LB, b), nrow = 4, ncol = 1, byrow = TRUE),
                     y  = y)
  
  # Extract the estimated parameters
  theta <- fit$par
  
  # Recompute the conditional variance
  sig2 <- f_ht(theta, y)
  
  # Compute the next-day ahead VaR for the Normal model
  VaR <- -qnorm(level) * tail(sqrt(sig2), 1)
  
  out <- list(VaR_Forecast = VaR, 
              ConditionalVariances = sig2, 
              GARCH_param = theta)
  
  out
}

f_nll <- function(theta, y) {
  ### Fonction which computes the negative log likelihood value 
  ### of a GARCH model with Normal errors
  #  INPUTS
  #   theta  : [vector] of parameters
  #   y      : [vector] (T x 1) of observations
  #  OUTPUTS
  #   nll    : [scalar] negative log likelihood value
  
  T <- length(y)
  
  # Compute the conditional variance of a GARCH(1,1) model
  sig2 <- f_ht(theta, y)
  
  # Consider the T values
  sig2 <- sig2[1:T]
  
  # Compute the loglikelihood
  ll <- sum(dnorm(x = y, 
                  mean = 0, 
                  sd = sqrt(sig2), 
                  log = TRUE))
  
  # Output the negative value
  nll <- -ll
  
  nll
}

f_ht <- function(theta, y)  {
  ### Function which computes the vector of conditional variance
  #  INPUTS
  #   theta : [vector] (3 x 1)
  #   y     : [vector] (T x 1) log-returns
  #  OUTPUTS 
  #   sig2  : [vector] (T+1 x 1) conditional variances
  
  # Extract the parameters
  a0 <- theta[1]
  a1 <- theta[2]
  b1 <- theta[3]
  
  T <- length(y)
  
  # Initialize the conditional variances
  sig2 <- rep(NA, T+1)
  
  # Start with unconditional variances
  sig2[1] <- a0 / (1 - a1 - b1)
  
  # Compute conditional variance at each step
  for (t in 2:(T+1)) {
    sig2[t] <- a0 + a1  * y[t-1]^2 + b1 * sig2[t-1]
  }
  
  
  sig2
}

f_forecast_var_gjr <- function(y, level) {
  ### Compute the VaR forecast of a gjrGARCH(1,1) 
  ### model with Normal errors at the desired risk level
  #  INPUTS
  #   y     : [vector] (T x 1) of observations (log-returns)
  #   level : [scalar] risk level (e.g. 0.95 for a VaR at the 95# risk level)
  #  OUTPUTS
  #   VaR   : [scalar] VaR forecast 
  #   sig2  : [vector] (T+1 x 1) conditional variances
  #   theta : [vector] GARCH parameters
  #  NOTE
  #   o the estimation is done by maximum likelihood
  
  # Fit a gjrGARCH(1,1) model with Normal errors
  # Starting values and bounds
  theta0 <- c(0.1 * var(y), 0.1, 0.8, 0.1)
  LB     <- c(0, 0, 0, 0) # Ax >= b au format ui %*% theta - ci >= 0 -> Ax >= b
  # Stationarity condition
  A      <-  c(0,-1,-1, -0.5) # Ax <= b au format ui %*% theta - ci >= 0 -> -Ax >= -b
  b      <-  -1
  
  # Run the optimization
  temp <- constrOptim(theta0,
                      f  = f_nll_gjr,
                      grad = NULL,
                      ui = matrix(c(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1, A),
                                  nrow = 5, ncol = 4, byrow = TRUE),
                      ci = matrix(c(LB, b), nrow = 5, ncol = 1, byrow = TRUE),
                      y  = y)
  theta <- temp$par
  # Recompute the conditional variance
  sig2 <- f_ht_gjr(theta, y)
  
  # Compute the next-day ahead VaR for the Normal model
  VaR <- -qnorm(level) * sqrt(tail(sig2,1))
  
  out <- list(VaR_Forecast = VaR, 
              ConditionalVariances = sig2, 
              GARCH_param = theta)
  
  out
}

f_nll_gjr <- function(theta, y) {
  ### Fonction which computes the negative log likelihood value 
  ### of a gjrGARCH model with Normal errors
  #  INPUTS
  #   theta  : [vector] of parameters
  #   y      : [vector] (T x 1) of observations
  #  OUTPUTS
  #   nll    : [scalar] negative log likelihood value
  
  T <- length(y)
  
  # Compute the conditional variance of a gjrGARCH(1,1) model
  sig2 <- f_ht_gjr(theta, y)
  
  # Consider the T values
  sig2 <- sig2[1:T]
  
  # Compute the loglikelihood
  ll <- sum(dnorm(x = y,
                  mean = 0,
                  sd = sqrt(sig2),
                  log = TRUE))
  
  # Output the negative value
  nll <- -ll
  
  nll
}

f_ht_gjr <- function(theta, y)  {
  ### Function which computes the vector of conditional variance
  #  INPUTS
  #   x0 : [vector] (3 x 1)
  #   y     : [vector] (T x 1) log-returns
  #  OUTPUTS 
  #   sig2  : [vector] (T+1 x 1) conditional variances
  
  # Extract the parameters
  a0 <- theta[1]
  a1 <- theta[2]
  b1 <- theta[3]
  g1 <- theta[4]
  
  T <- length(y)
  I <- 1*(y < 0) # I[t] = 1 if y[t] < 0, 0 otherwise
  
  # Initialize the conditional variances
  sig2 <- rep(NA, T + 1)
  
  # Start with unconditional variances
  sig2[1] <- a0 / (1 - a1 - b1 - g1/2)
  
  # Compute conditional variance at each step
  
  for (t in 2:(T + 1)) { 
    sig2[t] <- 
      a0 + a1 * y[t - 1]^2 + b1 * sig2[t - 1] + g1 * y[t - 1]^2 * I[t - 1] 
    
    
  }
  sig2
}

f_rolling_forecast_var <- function(y, level, model = 'sGARCH', 
                                   window_size, nb_forecast) {
  #### Function to get the rolling forecast of the VaR
  # INPUTS
  #   y     : [vector] (T x 1) log-returns
  #   level : [scalar] (1 x 1) confidence level
  #   model : [string] (1 x 1) 'sGARCH' or 'gjrGARCH'
  #   window_size : [scalar] (1 x 1) window size
  #   nb_forecast : [scalar] (1 x 1) number of forecasts
  # OUTPUTS
  #   forecast : [list] with the VaR forecast and the conditional variances
  
  y <- as.numeric(y)
  var_forecast <- rep(NA, nb_forecast)
  cv_forecast <- rep(NA, nb_forecast)
  
  for ( t in 1:(2*nb_forecast - window_size) ) {
    
    if (model == 'sGARCH') {
      
      
      forecast <- f_forecast_var(y = y[t:(t + window_size - 1)], 
                                 level = level)
    } else if (model == 'gjrGARCH'){
      
      forecast <- f_forecast_var_gjr(y = y[t:(t + window_size - 1)], 
                                     level = level)
    }
    else {
      stop('model should be sGARCH or gjrGARCH')
    }
    
    var_forecast[t] <-  forecast$VaR_Forecast
    
    cv_forecast[t] <- tail(forecast$ConditionalVariances, 1)
  }
  
  out <- list(roll_VaR_Forecast = var_forecast, 
              roll_ConditionalVariances = cv_forecast)
  
  return(out)
}

f_sim_biv_gauss <- function(model) {
  #### Function that simulate the bivariate Gaussian distribution
  # INPUTS:
  # model : [string] (1 x 1) 'sGARCH' or 'gjrGARCH'
  # OUTPUTS:
  # simulated : [data.frame] (T x 2) data frame with the simulated returns of the two assets
  
  
  # Standardisation of the residuals
  if (model == 'sGARCH') {
    
    std_SP500 <- 
      VaR_Rets_CV$SP500/sqrt(VaR_Rets_CV$cv_SP500)
    std_FTSE100 <- 
      VaR_Rets_CV$FTSE100/sqrt(VaR_Rets_CV$cv_FTSE100)
    
    # Compute the correlation between the standardised returns
    correlation <- cor(std_SP500, std_FTSE100)
    print(
      paste('The correlation between the sGARCH standardised returns is:',round(correlation, 3)))
    
  } else if (model == 'gjrGARCH') {
    
    
    std_gjrSP500 <- 
      VaR_Rets_CV$SP500/sqrt(VaR_Rets_CV$gjr_cv_SP500)
    std_gjrFTSE100 <- 
      VaR_Rets_CV$FTSE100/sqrt(VaR_Rets_CV$gjr_cv_FTSE100)
    
    correlation <- cor(std_gjrSP500, std_gjrFTSE100)
    print(
      paste('The correlation between the gjr-GARCH standardised returns is:', round(correlation, 3)))
    
    
  } else {
    stop('The model is not valid')
  }
  
  # Simulate a bivariate gaussian distribution with the estimated correlation
  set.seed(123)
  n <- 1000
  mu <- c(0, 0)
  sigma <- matrix(c(1, correlation, correlation, 1), nrow = 2)
  simulated <- MASS::mvrnorm(n = n, mu = mu, Sigma = sigma)
  plot(simulated[, 1], simulated[, 2], 
       main = 'Simulated Bivariate Gaussian Distribution')

  
  return(simulated)
}


