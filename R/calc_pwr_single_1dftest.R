
calc_pwr_single_1dftest <- function(K,            # Number of clusters in each arm
                                    m,            # Individuals per cluster
                                    alpha = 0.05, # Significance level
                                    beta1,        # Effect for outcome 1
                                    beta2,        # Effect for outcome 2
                                    varY1,        # Variance for outcome 1
                                    varY2,        # Variance for outcome 2
                                    rho01,        # ICC for outcome 1
                                    rho02,        # ICC for outcome 2
                                    rho1,         # Inter-subject between-endpoint ICC
                                    rho2          # Intra-subject between-endpoint ICC
                                    ){

  # Calculate critical value
  cv <- qchisq(p = alpha, df = 1, lower.tail = FALSE)

  # Calculate test statistic for first and second outcome
  Z1.sq <- (beta1^2)/(((2*varY1)/(K*m))*(1 + (m - 1)*rho01)) # Z1^2
  Z2.sq <- (beta2^2)/(((2*varY2)/(K*m))*(1 + (m - 1)*rho02)) # Z2^2

  # Calculate correlation between test statistics
  CorrZ1Z2 <- (rho2 + (m - 1)*rho1)/sqrt((1 + (m - 1)*rho01)*(1 + (m - 1)*rho02))

  # Calculate power
  lambda <- ((sqrt(Z1.sq) + sqrt(Z2.sq))^2)/(2*(1 + CorrZ1Z2))
  power <- round(1 - pchisq(cv, ncp = lambda, df = 1, lower.tail = TRUE), 4)

  return(power)
} # End calc_pwr_single_1dftest()

calc_K_single_1dftest <- function(power,        # Desired statistical power
                                  m,            # Individuals per cluster
                                  alpha = 0.05, # Significance level
                                  beta1,        # Effect for outcome 1
                                  beta2,        # Effect for outcome 2
                                  varY1,        # Variance for outcome 1
                                  varY2,        # Variance for outcome 2
                                  rho01,        # ICC for outcome 1
                                  rho02,        # ICC for outcome 2
                                  rho1,         # Inter-subject between-endpoint ICC
                                  rho2          # Intra-subject between-endpoint ICC
                                  ){

  # Calculate correlation between test statistics
  CorrZ1Z2 <- (rho2 + (m - 1)*rho1)/sqrt((1 + (m - 1)*rho01)*(1 + (m - 1)*rho02))

  # Non-centrality parameter for given power and alpha level
  ncp <- calc_ncp_chi2(alpha, power, df = 1)

  # Calculate K
  K <- ceiling((2*ncp*(1 + CorrZ1Z2))/
                 (sqrt((beta1^2)/((2*varY1/m)*(1+(m-1)*rho01))) +
                    sqrt((beta2^2)/((2*varY2/m)*(1+(m-1)*rho02))))^2)

  return(K)
} # End calc_K_single_1dftest()

calc_m_single_1dftest <- function(power,        # Desired statistical power
                                  K,            # Number of clusters in each arm
                                  alpha = 0.05, # Significance level
                                  beta1,        # Effect for outcome 1
                                  beta2,        # Effect for outcome 2
                                  varY1,        # Variance for outcome 1
                                  varY2,        # Variance for outcome 2
                                  rho01,        # ICC for outcome 1
                                  rho02,        # ICC for outcome 2
                                  rho1,         # Inter-subject between-endpoint ICC
                                  rho2          # Intra-subject between-endpoint ICC
                                  ){
  # Non-centrality parameter for given power and alpha level
  ncp <- calc_ncp_chi2(alpha, power, df = 1)

  # Function to solve for m
  power.eq.m <- function(m, para){ # Function for equation solve for m
    # Calculate test statistic for first and second outcome
    Z1.sq <- (para[1]^2)/((2*para[4]/(m*para[3]))*(1 + (m - 1)*para[6])) # Z1^2
    Z2.sq <- (para[2]^2)/((2*para[5]/(m*para[3]))*(1 + (m - 1)*para[7])) # Z1^2
    # Calculate Corr(Z1, Z2)
    CorrZ1Z2 <- (para[9] + (m - 1)*para[8])/
      sqrt((1+(m-1)*para[6])*(1+(m- 1)*para[7]))
    # Non-centrality parameter
    lambda <- ((sqrt(Z1.sq) + sqrt(Z2.sq))^2) / (2*(1 + CorrZ1Z2)) - para[10]
  }

  # Find m using multiroot function
  find.para.m <- multiroot(power.eq.m, start = 10,
                           para = c(beta1, beta2, K, varY1, varY2,
                                    rho01, rho02, rho1, rho2, ncp),
                           positive = TRUE)
  m <- ceiling(find.para.m$root)

  return(m)
} # End calc_m_single_1dftest()
