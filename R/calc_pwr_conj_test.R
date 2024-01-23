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

  # Function below requires ratio be defined as K1/K rather than K2/K1,
  # so define new ratio variable based on the one that was inputted by user
  r_alt <- 1/(r + 1)
  K_total <- ceiling(K/r_alt)
  betas = c(beta1, beta2)
  deltas = c(0, 0)
  Q = 2
  message("Using ", K_total, " as the total number of clusters for IU test.")

  # Variance of trt assignment
  sigmaz.square <- r_alt*(1 - r_alt)

  # Define small dependent functions -------------------------------------------
  # Dependent Function 1: Calculates covariance between betas
  calCovbetas <- function(vars, rho01, rho2, sigmaz.square, m, K){
    rho0k <- diag(rho01)
    sigmak.square <- (1+(m-1)*rho0k)*vars/(m*sigmaz.square)
    covMatrix <- diag(sigmak.square)
    for(row in 1:K ){
      for(col in 1:K){
        if(row != col){
          covMatrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*(rho2[row,col]+(m-1)*rho01[row,col])/(m*sigmaz.square)
        }
      }
    }
    return(covMatrix)
  }

  # Dependent Function 2: Calculates correlation between test statistics
  calCorWks <- function(vars, rho01, rho2, sigmaz.square, m, K){
    top <- calCovbetas(vars,rho01,rho2, sigmaz.square, m, K)
    wCor <- diag(K)
    for(row in 1:K){
      for(col in 1:K){
        if(row != col){
          wCor[row,col] <- top[row,col]/sqrt(top[row,row]*top[col,col])
        }
      }
    }
    return(wCor)
  }

  # Define sigma matrix
  sigmaks.sq <- diag(calCovbetas(vars = c(varY1, varY2),
                                 rho01 = matrix(c(rho01, rho1,
                                                  rho1, rho02),
                                                2, 2),
                                 rho2 = matrix(c(1, rho2,
                                                 rho2, 1),
                                               2, 2),
                                 sigmaz.square,
                                 m = m,
                                 K = Q # Number of outcomes
                                 ))

  # Define mean vector
  meanVector <- sqrt(K_total)*(betas - deltas)/sqrt(sigmaks.sq)
  wCor <- calCorWks(vars = c(varY1, varY2),
                    rho01 = matrix(c(rho01, rho1,
                                     rho1, rho02),
                                   2, 2),
                    rho2 = matrix(c(1, rho2,
                                    rho2, 1),
                                  2, 2),
                    sigmaz.square,
                    m = m,
                    K = Q # Number of outcomes
                    )

  # Calculate critical value
  criticalValue <- qt(p = (1 - alpha), df = (K_total - 2*Q))

  # Calculate power
  power <- pmvt(lower = rep(criticalValue, Q),
                upper = rep(Inf, Q), df = (K_total - 2*Q),
                sigma = wCor, delta = meanVector)[1]

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
                                                      r = r
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
                                                      r = r
                                                      ))
    if(m > 100000){
      m <- Inf
      message("Cannot find large enough 'm' to reach study specifications for IU test. Please lower power or increase value for 'K'.")
      break
    }
  }
  return(ceiling(m))

} # End calc_m_conj_test()

