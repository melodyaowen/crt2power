#' Calculate statistical power for a cluster randomized trial with co-primary endpoints using a disjunctive 2-DF test approach.
#'
#' @description
#' Allows user to calculate the statistical power of a hybrid type 2 cluster randomized trial given a set of study design input values, including the number of clusters in each trial arm, and cluster size. Uses the disjunctive 2-DF test approach.
#'
#' @param K Number of clusters in each arm; numeric.
#' @param m Individuals per cluster; numeric.
#' @param alpha Type I error rate; numeric.
#' @param beta1 Effect size for the first outcome; numeric.
#' @param beta2 Effect size for the second outcome; numeric.
#' @param varY1 Total variance for the first outcome; numeric.
#' @param varY2 Total variance for the second outcome; numeric.
#' @param rho01 Correlation of the first outcome for two different individuals in the same cluster; numeric.
#' @param rho02 Correlation of the second outcome for two different individuals in the same cluster; numeric.
#' @param rho1 Correlation between the first and second outcomes for two individuals in the same cluster; numeric.
#' @param rho2 Correlation between the first and second outcomes for the same individual; numeric.
#' @returns A data frame of numerical values.
#' @examples
#' calc_pwr_disj_2dftest(K = 15, m = 300, alpha = 0.05, beta1 = 0.1, beta2 = 0.1, varY1 = 0.23, varY2 = 0.25, rho01 = 0.025, rho02 = 0.025, rho1 = 0.01, rho2  = 0.05)
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
} # End calc_pwr_disj_2dftest()








#' Calculate required number of clusters per treatment group for a cluster randomized trial with co-primary endpoints using a disjunctive 2-DF test approach.
#'
#' @description
#' Allows user to calculate the number of clusters per treatment arm of a hybrid type 2 cluster randomized trial given a set of study design input values, including the statistical power, and cluster size. Uses the disjunctive 2-DF test approach.
#'
#' @param power Desired statistical power in decimal form; numeric.
#' @param m Individuals per cluster; numeric.
#' @param alpha Type I error rate; numeric.
#' @param beta1 Effect size for the first outcome; numeric.
#' @param beta2 Effect size for the second outcome; numeric.
#' @param varY1 Total variance for the first outcome; numeric.
#' @param varY2 Total variance for the second outcome; numeric.
#' @param rho01 Correlation of the first outcome for two different individuals in the same cluster; numeric.
#' @param rho02 Correlation of the second outcome for two different individuals in the same cluster; numeric.
#' @param rho1 Correlation between the first and second outcomes for two individuals in the same cluster; numeric.
#' @param rho2 Correlation between the first and second outcomes for the same individual; numeric.
#' @returns A data frame of numerical values.
#' @examples
#' calc_K_disj_2dftest(power = 0.8, m = 300, alpha = 0.05, beta1 = 0.1, beta2 = 0.1, varY1 = 0.23, varY2 = 0.25, rho01 = 0.025, rho02 = 0.025, rho1 = 0.01, rho2  = 0.05)
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
} # End calc_K_disj_2dftest()




#' Calculate cluster size for a cluster randomized trial with co-primary endpoints using a disjunctive 2-DF test approach.
#'
#' @description
#' Allows user to calculate the cluster size of a hybrid type 2 cluster randomized trial given a set of study design input values, including the number of clusters in each trial arm, and statistical power. Uses the disjunctive 2-DF test approach.
#'
#' @param power Desired statistical power in decimal form; numeric.
#' @param K Number of clusters in each arm; numeric.
#' @param alpha Type I error rate; numeric.
#' @param beta1 Effect size for the first outcome; numeric.
#' @param beta2 Effect size for the second outcome; numeric.
#' @param varY1 Total variance for the first outcome; numeric.
#' @param varY2 Total variance for the second outcome; numeric.
#' @param rho01 Correlation of the first outcome for two different individuals in the same cluster; numeric.
#' @param rho02 Correlation of the second outcome for two different individuals in the same cluster; numeric.
#' @param rho1 Correlation between the first and second outcomes for two individuals in the same cluster; numeric.
#' @param rho2 Correlation between the first and second outcomes for the same individual; numeric.
#' @returns A data frame of numerical values.
#' @examples
#' calc_m_disj_2dftest(power = 0.8, K = 15, alpha = 0.05, beta1 = 0.1, beta2 = 0.1, varY1 = 0.23, varY2 = 0.25, rho01 = 0.025, rho02 = 0.025, rho1 = 0.01, rho2  = 0.05)
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
} # End calc_m_disj_2dftest()

