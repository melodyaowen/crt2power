#' Calculate statistical power for a cluster randomized trial with co-primary endpoints using the conjunctive intersection-union test approach.
#'
#' @import devtools
#' @import kableExtra
#' @import knitr
#' @import MASS
#' @import pracma
#' @import rootSolve
#' @import tidyverse
#' @import tableone
#' @import foreach
#' @import mvtnorm
#' @import tibble
#' @import dplyr
#' @import tidyr
#'
#' @description
#' Allows user to calculate the statistical power of a hybrid type 2 cluster randomized trial given a set of study design input values, including the number of clusters in each trial arm, and cluster size. Uses the conjunctive intersection-union test approach. Code is adapted from "calPower_ttestIU()" from https://github.com/siyunyang/coprimary_CRT written by Siyun Yang.
#'
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
#' @param cv Cluster variation parameter, set to 0 if assuming all cluster sizes are equal; numeric.
#' @param deltas Vector of non-inferiority margins, set to delta_1 = delta_2 = 0; numeric vector.
#' @param dist Specification of which distribution to base calculation on, either 'T' for T-Distribution or 'MVN' for Multivariate Normal Distribution. Default is T-Distribution.
#' @returns A numerical value.
#' @examples
#' calc_pwr_conj_test(K = 15, m = 300, alpha = 0.05, beta1 = 0.1, beta2 = 0.1, varY1 = 0.23, varY2 = 0.25, rho01 = 0.025, rho02 = 0.025, rho1 = 0.01, rho2  = 0.05)
#' @export
calc_pwr_conj_test <- function(K,            # Number of clusters in treatment arm
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
                               r = 1,        # Treatment allocation ratio
                               cv = 0,       # If equal cluster size, cv=0
                               deltas = c(0,0),
                               dist = "T"    # Distribution to be used
                               ){

  # Check that input values are valid
  if(!is.numeric(c(K, m, alpha, beta1, beta2, varY1, varY2, rho01, rho02, rho1, rho2, r, cv))){
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

  # Helper functions requires ratio be defined as K1/K rather than K2/K1,
  # so define new ratio variable based on the one that was inputted by user
  r_alt <- 1/(r + 1)
  K_total <- ceiling(K/r_alt)
  betas = c(beta1, beta2)
  Q = 2
  message("Using ", K_total, " as the total number of clusters for IU test.")

  # Variance of trt assignment
  sigmaz.square <- r_alt*(1 - r_alt)

  # Helper Function 1: Construct covariance matrix Sigma_E for Y_k -------------
  constrRiE <- function(rho01, rho2, Q, vars){
    rho0q <- diag(rho01)
    SigmaE_Matrix <- diag((1-rho0q)*vars)
    for(row in 1:Q){
      for(col in 1:Q){
        if(row != col){
          SigmaE_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*(rho2[row,col]-rho01[row,col])
        }
      }
    }
    # Check for matrix positive definite
    if(min(eigen(SigmaE_Matrix)$values) <= 1e-08){
      print("Warning: the resulting covariance matrix Sigma_E is not positive definite. Check the inputs for the correlation values.")
    }
  return(SigmaE_Matrix)
  }

  # Helper Function 2: Construct covariance matrix Sigma_phi for Y_k -----------
  constrRiP <- function(rho01, Q, vars){
    rho0q <- diag(rho01)
    SigmaP_Matrix <- diag(rho0q*vars)
    for(row in 1:Q){
      for(col in 1:Q){
        if(row != col){
          SigmaP_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*rho01[row,col]
        }
      }
    }
    # Check for matrix positive definite
    if(min(eigen(SigmaP_Matrix)$values) <= 1e-08){
      print("Warning: the resulting covariance matrix Sigma_phi is not positive definite. Check the input of rho01 and rho2.")
    }
  return(SigmaP_Matrix)
  }

  # Helper Function 3: Calculate covariance between betas ----------------------
  calCovbetas <- function(vars, rho01, rho2, cv, sigmaz.square, m, Q){
    sigmaE <- constrRiE(rho01, rho2, Q, vars)
    sigmaP <- constrRiP(rho01, Q, vars)
    tmp <- solve(diag(1,Q) - cv^2*(m*sigmaP %*% solve(sigmaE + m*sigmaP) %*% sigmaE %*% solve(sigmaE + m*sigmaP)))
    covMatrix <- 1/(m*sigmaz.square)*(sigmaE + m*sigmaP)%*%tmp
    covMatrix <- (covMatrix + t(covMatrix))/2  # symmerize the off-diagonal
    return(covMatrix)
  }

  # Helper Function 4: Calculate correlation between test statistics -----------
  calCorWks <- function(vars, rho01, rho2, sigmaz.square, cv, m, Q){
    top <- calCovbetas(vars, rho01, rho2, cv, sigmaz.square, m, Q)
    wCor <- diag(Q)
    for(row in 1:Q){
      for(col in 1:Q){
        if(row != col){
          wCor[row,col] <- top[row,col]/sqrt(top[row,row]*top[col,col])
        }
      }
    }
    return(wCor)
  }

  # Define necessary parameters
  sigmaks.sq <- diag(calCovbetas(vars = c(varY1, varY2),
                                 rho01 = matrix(c(rho01, rho1,
                                                  rho1, rho02),
                                                2, 2),
                                 rho2 = matrix(c(1, rho2,
                                                 rho2, 1),
                                               2, 2),
                                 cv = cv,
                                 sigmaz.square = sigmaz.square,
                                 m = m,
                                 Q = Q))
  meanVector <- sqrt(K_total)*(betas - deltas)/sqrt(sigmaks.sq)
  wCor <- calCorWks(vars = c(varY1, varY2),
                    rho01 = matrix(c(rho01, rho1,
                                     rho1, rho02),
                                   2, 2),
                    rho2 = matrix(c(1, rho2,
                                    rho2, 1),
                                  2, 2),
                    sigmaz.square = sigmaz.square,
                    cv = cv,
                    m = m,
                    Q = Q)

  if(dist == "T"){ # Using T-distribution
    # Calculate critical value and power
    criticalValue <- qt(p = (1 - alpha),
                        df = (K_total - 2*Q))
    power <- pmvt(lower = rep(criticalValue, Q),
                  upper = rep(Inf, Q),
                  df = (K_total - 2*Q),
                  sigma = wCor,
                  delta = meanVector)[1]
  } else if(dist == "MVN"){ # Using multivariate normal distribution
    # Calculate critical value and power
    criticalValue <- qnorm(p = 1 - alpha,
                           mean = 0,
                           sd = 1)
    power <- pmvnorm(lower = rep(criticalValue, Q),
                     upper = rep(Inf, Q),
                     corr = wCor,
                     mean = meanVector)[1]
  } else{
    stop("Please choose a valid input parameter for 'dist', either 'T' for T-distribution or 'MVN' for Multivariate Normal Distribution.")
  }
  return(round(power, 4))
} # End calc_pwr_conj_test()








#' Calculate required number of clusters per treatment group for a cluster randomized trial with co-primary endpoints using the conjunctive intersection-union test approach.
#'
#' @description
#' Allows user to calculate the required number of clusters per treatment group of a hybrid type 2 cluster randomized trial given a set of study design input values, including the statistical power, and cluster size. Uses the conjunctive intersection-union test approach.Code is adapted from "calSampleSize_ttestIU()" from https://github.com/siyunyang/coprimary_CRT written by Siyun Yang.
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
#' @param cv Cluster variation parameter, set to 0 if assuming all cluster sizes are equal; numeric.
#' @param deltas Vector of non-inferiority margins, set to delta_1 = delta_2 = 0; numeric vector.
#' @param dist Specification of which distribution to base calculation on, either 'T' for T-Distribution or 'MVN' for Multivariate Normal Distribution. Default is T-Distribution.
#' @returns A data frame of numerical values.
#' @examples
#' calc_K_conj_test(power = 0.8, m = 300, alpha = 0.05, beta1 = 0.1, beta2 = 0.1, varY1 = 0.23, varY2 = 0.25, rho01 = 0.025, rho02 = 0.025, rho1 = 0.01, rho2  = 0.05)
#' @export
calc_K_conj_test <- function(power,        # Desired statistical power
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
                             r = 1,        # Treatment allocation ratio
                             cv = 0,
                             deltas = c(0,0),
                             dist = "T"
                             ){

  # Check that input values are valid
  if(!is.numeric(c(power, m, alpha, beta1, beta2, varY1, varY2, rho01, rho02, rho1, rho2, r, cv))){
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

  # Some calculations require ratio be defined as K1/K rather than K2/K1,
  # so define new ratio variable based on the one that was inputted by user
  r_alt = 1/(r + 1)

  lowerBound <- 1
  upperBound <- 10000
  repeat{
    middle <- floor((lowerBound + upperBound)/2)
    power_temp <- suppressMessages(calc_pwr_conj_test(K = middle,
                                                      m = m,
                                                      alpha = alpha,
                                                      beta1 = beta1,
                                                      beta2 = beta2,
                                                      varY1 = varY1,
                                                      varY2 = varY2,
                                                      rho01 = rho01,
                                                      rho02 = rho02,
                                                      rho1 = rho1,
                                                      rho2 = rho2,
                                                      r = r,
                                                      cv = cv,
                                                      deltas = deltas,
                                                      dist = dist
                                                      ))
    if(power_temp < power){
      lowerBound <- middle
    }
    if(power_temp > power){
      upperBound <- middle
    }
    if(power_temp == power){
      return(middle)
      break
    }
    if((lowerBound-upperBound) == -1){
      K <- upperBound
      break
    }
  }

  if(r == 1){ # Case when treatment allocation ratio is equal
    K <- tibble(`Treatment (K)` = ceiling(K),
                `Control (K)` = ceiling(K))
  } else{ # Case for unequal treatment allocation, make sure we round up
    K1 <- ceiling(K) # ceiling(r_alt*K_total)
    K2 <- ceiling(r*K) # ceiling(r*K1)
    K <- tibble(`Treatment (K1)` = K1,
                `Control (K2)` = K2)
  }

  return(ceiling(K))
} # End calc_K_conj_test()









#' Calculate cluster size for a cluster randomized trial with co-primary endpoints using the conjunctive intersection-union test approach.
#'
#' @description
#' Allows user to calculate the cluster size of a hybrid type 2 cluster randomized trial given a set of study design input values, including the number of clusters in each trial arm, and statistical power. Uses the conjunctive intersection-union test approach.
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
#' @param cv Cluster variation parameter, set to 0 if assuming all cluster sizes are equal; numeric.
#' @param deltas Vector of non-inferiority margins, set to delta_1 = delta_2 = 0; numeric vector.
#' @param dist Specification of which distribution to base calculation on, either 'T' for T-Distribution or 'MVN' for Multivariate Normal Distribution. Default is T-Distribution.
#' @returns A numerical value.
#' @examples
#' calc_m_conj_test(power = 0.8, K = 15, alpha = 0.05, beta1 = 0.1, beta2 = 0.1, varY1 = 0.23, varY2 = 0.25, rho01 = 0.025, rho02 = 0.025, rho1 = 0.01, rho2  = 0.05)
#' @export
calc_m_conj_test <- function(power,        # Desired statistical power
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
                             r = 1,        # Treatment allocation ratio
                             cv = 0,
                             deltas = c(0, 0),
                             dist = "T"
                             ){

  # Check that input values are valid
  if(!is.numeric(c(power, K, alpha, beta1, beta2, varY1, varY2, rho01, rho02, rho1, rho2, r, cv))){
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

  # Function below requires ratio be defined as K1/K rather than K2/K1,
  # so define new ratio variable based on the one that was inputted by user
  r_alt = 1/(r + 1)
  K_total <- ceiling(K/r_alt)
  message("Using ", K_total, " as the total number of clusters for IU test.")

  # Initialize m and predictive power
  m <- 0
  pred.power <- 0
  while(pred.power < power){
    m <- m + 1
    pred.power <- suppressMessages(calc_pwr_conj_test(K = K,
                                                      m = m,
                                                      alpha = alpha,
                                                      beta1 = beta1,
                                                      beta2 = beta2,
                                                      varY1 = varY1,
                                                      varY2 = varY2,
                                                      rho01 = rho01,
                                                      rho02 = rho02,
                                                      rho1 = rho1,
                                                      rho2 = rho2,
                                                      r = r,
                                                      cv = cv,
                                                      deltas = deltas,
                                                      dist = dist
                                                      ))
    if(m > 100000){
      m <- Inf
      message("Cannot find large enough 'm' to reach study specifications for IU test. Please lower power or increase value for 'K'.")
      break
    }
  }
  return(ceiling(m))

} # End calc_m_conj_test()

