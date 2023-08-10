calc_pwr_disj_2dftest <- function(K,            # Number of clusters in each arm
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

  # Critical value, df = 2
  cv <- qchisq(1 - alpha, df = 2, ncp = 0, lower.tail = TRUE, log.p = FALSE)

  # Calculate VIFs
  VIF1 <- 1 + (m - 1)*rho01
  VIF2 <- 1 + (m - 1)*rho02
  VIF12 <- rho2 + (m - 1)*rho1

  # Power calculation
  lambda <- K*m*((beta1^2)*varY2*VIF2 + (beta2^2)*varY1*VIF1 -
                   2*beta1*beta2*sqrt(varY1)*sqrt(varY2)*VIF12)/
    (2*varY1*varY2*(VIF1*VIF2-VIF12^2))

  power <- round(1 - pchisq(cv, df = 2, ncp = lambda, lower.tail = TRUE), 4)

  return(power)
} # End calc_pwr_disj_2dftest

calc_K_disj_2dftest <- function(power,        # Desired statistical power
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

  # Non-centrality parameter for given power and alpha level
  ncp <- calc_ncp_chi2(alpha, power, df = 2)

  # Calculate VIFs
  VIF1 <- 1 + (m - 1)*rho01
  VIF2 <- 1 + (m - 1)*rho02
  VIF12 <- rho2 + (m - 1)*rho1

  # Calculate K
  K <- ceiling((ncp*2*varY1*varY2*(VIF1*VIF2-VIF12^2))/
                 (m*((beta1^2)*varY2*VIF2 + (beta2^2)*varY1*VIF1 -
                       2*beta1*beta2*sqrt(varY1)*sqrt(varY2)*VIF12)))

  return(K)
} # End calc_K_disj_2dftest

calc_m_disj_2dftest <- function(power,        # Desired statistical power
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
  ncp <- calc_ncp_chi2(alpha, power, df = 2)

  # m for Method 4
  # Temporary equation to use multiroot() on to solve for m
  find.para.m <- function(m, para){
    VIF1 <- 1 + (m - 1)*para[6]
    VIF2 <- 1 + (m - 1)*para[7]
    VIF12 <- para[9] + (m - 1)*para[8]
    lambda <- para[3]*m*(para[1]^2*para[5]*VIF2 + para[2]^2*para[4]*VIF1 -
                           2*para[1]*para[2]*sqrt(para[4]*para[5])*VIF12)/
      (2*para[4]*para[5]*(VIF1*VIF2 - VIF12^2)) - para[10]
  }

  # Calculate m
  find.para.m <- multiroot(find.para.m, start = 3,
                           para = c(beta1, beta2, K, varY1, varY2,
                                    rho01, rho02, rho1, rho2, ncp),
                           positive = TRUE)
  m <- ceiling(find.para.m$root) # Cluster size required

  return(m)
} # End calc_m_disj_2dftest

