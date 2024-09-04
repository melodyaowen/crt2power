#' Calculate statistical power for a cluster-randomized trial with co-primary endpoints using a combined outcomes approach.
#'
#' @import devtools
#' @import knitr
#' @import rootSolve
#' @import tidyverse
#' @import tableone
#' @import foreach
#' @import mvtnorm
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @importFrom stats uniroot dchisq pchisq qchisq rchisq df pf qf rf dt pt qt rt dnorm pnorm qnorm rnorm
#' @description
#' Allows user to calculate the statistical power of a cluster-randomized trial with two co-primary outcomes given a set of study design input values, including the number of clusters in each trial arm, and cluster size. Uses a combined outcomes approach where the two outcome effects are summed together.
#'
#' @param dist Specification of which distribution to base calculation on, either 'Chi2' for Chi-Squared or 'F' for F-Distribution.
#' @param K Number of clusters in treatment arm, and control arm under equal allocation; numeric.
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
#' @param r Treatment allocation ratio - K2 = rK1 where K1 is number of clusters in experimental group; numeric.
#' @returns A numerical value.
#' @examples
#' calc_pwr_comb_outcome(K = 15, m = 300, alpha = 0.05,
#' beta1 = 0.1, beta2 = 0.1, varY1 = 0.23, varY2 = 0.25,
#' rho01 = 0.025, rho02 = 0.025, rho1 = 0.01, rho2  = 0.05)
#' @export
calc_pwr_comb_outcome <- function(dist = "Chi2",# Distribution to base calculation from
                                  K,            # Number of clusters in treatment arm
                                  m,            # Individuals per cluster
                                  alpha = 0.05, # Significance level
                                  beta1,        # Effect for outcome 1
                                  beta2,        # Effect for outcome 2
                                  varY1,        # Variance for outcome 1
                                  varY2,        # Variance for outcome 2
                                  rho01,        # ICC for outcome 1
                                  rho02,        # ICC for outcome 2
                                  rho1,         # Inter-subject between-endpoint ICC
                                  rho2,         # Intra-subject between-endpoint ICC
                                  r = 1         # Treatment allocation ratio
                                  ){
  # Check that input values are valid
  if(!is.numeric(c(K, m, alpha, beta1, beta2, varY1, varY2, rho01, rho02, rho1, rho2, r))){
    stop("All input parameters must be numeric values.")
  }
  if(r <= 0){
    stop("Treatment allocation ratio should be a number greater than 0.")
  }
  if(K < 1 | K != round(K)){
    stop("'K' must be a positive whole number.")
  }
  if(m < 1 | m != round(m)){
    stop("'m' must be a positive whole number.")
  }

  # Defining necessary parameters based on input values
  r_alt <- 1/(r + 1)
  Q <- 2 # Number of outcomes, could extend this to more than 2 in the future
  K_total <- ceiling(K/r_alt) # Total number of clusters

  # Calculate combined outcome effect size
  betaC <- beta1 + beta2

  # Calculate variance for combined outcome
  varYC <- round(varY1 + varY2 + 2*rho2*sqrt(varY1)*sqrt(varY2), 4)

  # Calculate ICC for combined outcome
  rho0C  <- (rho01*varY1 + rho02*varY2 + 2*rho1*sqrt(varY1*varY2))/
    (varY1 + varY2 + 2*rho2*sqrt(varY1*varY2))

  # Lambda value for the combined outcome
  lambda <- (betaC^2)/((1 + 1/r)*(varYC/(K*m))*(1 + (m - 1)*rho0C))

  if(dist == "Chi2"){
    # Find critical value
    cv <- qchisq(p = alpha, df = 1, lower.tail = FALSE)
    # Calculate power
    power <- round(1 - pchisq(cv, 1, ncp = lambda, lower.tail = TRUE), 4)
  } else if(dist == "F"){
    Fscore <- qf(1 - alpha, df1 = 1, df2 = K_total - 2*Q, ncp = 0,
                 lower.tail = TRUE, log.p = FALSE)
    power <- round(1 - pf(Fscore, ncp = lambda,
                          df1 = 1, df2 = K_total - 2*Q,
                          lower.tail = TRUE, log.p = FALSE), 4)
  } else{
    stop("Please choose a valid input parameter for 'dist', either 'Chi2' for Chi-Square or 'F' for F-distribution.")
  }

  return(power)
} # End calc_pwr_comb_outcome()








#' Calculate required number of clusters per treatment group for a cluster-randomized trial with co-primary endpoints using a combined outcomes approach.
#'
#' @description
#' Allows user to calculate the number of clusters per treatment arm of a cluster-randomized trial with two co-primary outcomes given a set of study design input values, including the number of clusters in each trial arm, and cluster size. Uses a combined outcomes approach where the two outcome effects are summed together.
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
#' @param r Treatment allocation ratio - K2 = rK1 where K1 is number of clusters in experimental group; numeric.
#' @returns A data frame of numerical values.
#' @examples
#' calc_K_comb_outcome(power = 0.8, m = 300, alpha = 0.05,
#' beta1 = 0.1, beta2 = 0.1, varY1 = 0.23, varY2 = 0.25,
#' rho01 = 0.025, rho02 = 0.025, rho1 = 0.01, rho2  = 0.05)
#' @export
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
                                rho2,         # Intra-subject between-endpoint ICC
                                r = 1         # Treatment allocation ratio
                                ){

  # Check that input values are valid
  if(!is.numeric(c(power, m, alpha, beta1, beta2, varY1, varY2, rho01, rho02, rho1, rho2, r))){
    stop("All input parameters must be numeric values.")
  }
  if(r <= 0){
    stop("Treatment allocation ratio should be a number greater than 0.")
  }
  if(power > 1 | power < 0){
    stop("'power' must be a number between 0 and 1.")
  }
  if(m < 1 | m != round(m)){
    stop("'m' must be a positive whole number.")
  }

  # Calculate combined outcome effect size
  betaC <- beta1 + beta2

  # Calculate variance for combined outcome
  varYC <- round(varY1 + varY2 + 2*rho2*sqrt(varY1)*sqrt(varY2), 4)

  # Calculate ICC for combined outcome
  rho0C  <- (rho01*varY1 + rho02*varY2 + 2*rho1*sqrt(varY1*varY2))/
    (varY1 + varY2 + 2*rho2*sqrt(varY1*varY2))

  # Find non-centrality parameter corresponding to desired power
  ncp <- calc_ncp_chi2(alpha, power, df = 1)

  # K for Method 2
  if(r == 1){ # When treatment allocation is even
    K1 <- ceiling((2*ncp*varYC*(1 + (m - 1)*rho0C))/(m*(betaC^2)))
    K <- tibble(`Treatment (K)` = K1,
                `Control (K)` = K1)
  } else{ # Unequal treatment allocation
    K1 <- ceiling(((1 + 1/r)*ncp*varYC*(1 + (m - 1)*rho0C))/(m*(betaC^2)))
    #K2 <- ceiling(((1 + r)*ncp*varYC*(1 + (m - 1)*rho0C))/(m*(betaC^2)))
    K2 <- ceiling(r*K1)
    K <- tibble(`Treatment (K1)` = K1,
                `Control (K2)` = K2)
  }

  return(K)
} # End calc_K_comb_outcome()









#' Calculate cluster size for a cluster-randomized trial with co-primary endpoints using a combined outcomes approach.
#'
#' @description
#' Allows user to calculate the cluster size of a cluster-randomized trial with two co-primary endpoints given a set of study design input values, including the number of clusters in each trial arm, and statistical power. Uses a combined outcomes approach where the two outcome effects are summed together.
#'
#' @param power Desired statistical power in decimal form; numeric.
#' @param K Number of clusters in treatment arm, and control arm under equal allocation; numeric.
#' @param alpha Type I error rate; numeric.
#' @param beta1 Effect size for the first outcome; numeric.
#' @param beta2 Effect size for the second outcome; numeric.
#' @param varY1 Total variance for the first outcome; numeric.
#' @param varY2 Total variance for the second outcome; numeric.
#' @param rho01 Correlation of the first outcome for two different individuals in the same cluster; numeric.
#' @param rho02 Correlation of the second outcome for two different individuals in the same cluster; numeric.
#' @param rho1 Correlation between the first and second outcomes for two individuals in the same cluster; numeric.
#' @param rho2 Correlation between the first and second outcomes for the same individual; numeric.
#' @param r Treatment allocation ratio - K2 = rK1 where K1 is number of clusters in experimental group; numeric.
#' @returns A numerical value.
#' @examples
#' calc_m_comb_outcome(power = 0.8, K = 15, alpha = 0.05,
#' beta1 = 0.1, beta2 = 0.1, varY1 = 0.23, varY2 = 0.25,
#' rho01 = 0.025, rho02 = 0.025, rho1 = 0.01, rho2  = 0.05)
#' @export
calc_m_comb_outcome <- function(power,        # Desired statistical power
                                K,            # Number of clusters in treatment arm
                                alpha = 0.05, # Significance level
                                beta1,        # Effect for outcome 1
                                beta2,        # Effect for outcome 2
                                varY1,        # Variance for outcome 1
                                varY2,        # Variance for outcome 2
                                rho01,        # ICC for outcome 1
                                rho02,        # ICC for outcome 2
                                rho1,         # Inter-subject between-endpoint ICC
                                rho2,         # Intra-subject between-endpoint ICC
                                r = 1         # Treatment allocation ratio
                                ){

  # Check that input values are valid
  if(!is.numeric(c(power, K, alpha, beta1, beta2, varY1, varY2, rho01, rho02, rho1, rho2, r))){
    stop("All input parameters must be numeric values.")
  }
  if(r <= 0){
    stop("Treatment allocation ratio should be a number greater than 0.")
  }
  if(power > 1 | power < 0){
    stop("'power' must be a number between 0 and 1.")
  }
  if(K < 1 | K != round(K)){
    stop("'K' must be a positive whole number.")
  }

  # Calculate combined outcome effect size
  betaC <- beta1 + beta2

  # Calculate variance for combined outcome
  varYC <- round(varY1 + varY2 + 2*rho2*sqrt(varY1)*sqrt(varY2), 4)

  # Calculate ICC for combined outcome
  rho0C  <- (rho01*varY1 + rho02*varY2 + 2*rho1*sqrt(varY1*varY2))/
    (varY1 + varY2 + 2*rho2*sqrt(varY1*varY2))

  # Find non-centrality parameter corresponding to desired power
  ncp <- calc_ncp_chi2(alpha, power, df = 1)

  # m for Method 2
  m <- ceiling(((1 + 1/r)*ncp*varYC*(1 - rho0C))/((betaC^2)*K - ((1 + 1/r)*ncp*varYC*rho0C)))

  return(m)
} # End calc_m_comb_outcome()
