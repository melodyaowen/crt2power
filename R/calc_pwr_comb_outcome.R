
calc_pwr_comb_outcome <- function(K,            # Number of clusters in each arm
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

  # Calculate combined outcome effect size
  betaC <- beta1 + beta2

  # Calculate variance for combined outcome
  varYC <- round(varY1 + varY2 + 2*rho2*sqrt(varY1)*sqrt(varY2), 2)

  # Calculate ICC for combined outcome
  rho0C  <- (rho01*varY1 + rho02*varY2 + 2*rho1*sqrt(varY1*varY2))/
    (varY1 + varY2 + 2*rho2*sqrt(varY1*varY2))

  # Find critical value
  cv <- qchisq(p = alpha, df = 1, lower.tail = FALSE)

  # Power for Method 2
  lambda <- (betaC^2)/(2*(varYC/(K*m))*(1 + (m - 1)*rho0C))
  power <- round(1 - pchisq(cv, 1, ncp = lambda, lower.tail = TRUE), 4)

  return(power)
} # End calc_pwr_comb_outcome()

calc_K_comb_outcome <- function(power,        # Desired statistical power
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

  # Calculate combined outcome effect size
  betaC <- beta1 + beta2

  # Calculate variance for combined outcome
  varYC <- round(varY1 + varY2 + 2*rho2*sqrt(varY1)*sqrt(varY2), 2)

  # Calculate ICC for combined outcome
  rho0C  <- (rho01*varY1 + rho02*varY2 + 2*rho1*sqrt(varY1*varY2))/
    (varY1 + varY2 + 2*rho2*sqrt(varY1*varY2))

  # Find non-centrality parameter corresponding to desired power
  ncp <- calc_ncp_chi2(alpha, power, df = 1)

  # K for Method 2
  K <- ceiling((2*ncp*varYC*(1 + (m - 1)*rho0C))/(m*(betaC^2)))

  return(K)
} # End calc_K_comb_outcome()

calc_m_comb_outcome <- function(power,        # Desired statistical power
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

  # Calculate combined outcome effect size
  betaC <- beta1 + beta2

  # Calculate variance for combined outcome
  varYC <- round(varY1 + varY2 + 2*rho2*sqrt(varY1)*sqrt(varY2), 2)

  # Calculate ICC for combined outcome
  rho0C  <- (rho01*varY1 + rho02*varY2 + 2*rho1*sqrt(varY1*varY2))/
    (varY1 + varY2 + 2*rho2*sqrt(varY1*varY2))

  # Find non-centrality parameter corresponding to desired power
  ncp <- calc_ncp_chi2(alpha, power, df = 1)

  # m for Method 2
  m <- ceiling((2*ncp*varYC*(1 - rho0C))/((betaC^2)*K - (2*ncp*varYC*rho0C)))

  return(m)
} # End calc_m_comb_outcome()
