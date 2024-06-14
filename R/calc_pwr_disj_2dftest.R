#' Calculate statistical power for a cluster-randomized trial with co-primary endpoints using a disjunctive 2-DF test approach.
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
#' @importFrom stats uniroot dchisq pchisq qchisq rchisq df pf qf rf dt pt qt rt
#'
#' @description
#' Allows user to calculate the statistical power of a cluster-randomized trial with two co-primary outcomes given a set of study design input values, including the number of clusters in each trial arm, and cluster size. Uses the disjunctive 2-DF test approach. Code is adapted from "calPower_omnibus()" from https://github.com/siyunyang/coprimary_CRT written by Siyun Yang.
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
#' calc_pwr_disj_2dftest(K = 15, m = 300, alpha = 0.05, beta1 = 0.1, beta2 = 0.1, varY1 = 0.23, varY2 = 0.25, rho01 = 0.025, rho02 = 0.025, rho1 = 0.01, rho2  = 0.05)
#' @export
calc_pwr_disj_2dftest <- function(dist = "Chi2",# Distribution to base calculation from
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
    stop("All input parameters must be numeric values (with the exception of 'dist').")
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

  # Define small dependent functions -------------------------------------------
  # Dependent Function 1: Calculates covariance between betas
  calCovbetas <- function(vars, rho01, rho2, cv, sigmaz.square, m, Q){
    sigmaE <- constrRiE(rho01, rho2, Q, vars)
    sigmaP <- constrRiP(rho01, Q, vars)
    tmp <- solve(diag(1,Q) - cv^2*(m*sigmaP %*% solve(sigmaE + m*sigmaP) %*% sigmaE %*% solve(sigmaE + m*sigmaP) ))
    covMatrix <- 1/(m*sigmaz.square)*(sigmaE + m*sigmaP)%*%tmp
    covMatrix <- (covMatrix + t(covMatrix))/2  # symmerize the off-diagonal
    return(covMatrix)
  }

  # Dependent Function 2: Constructs covariance matrix Sigma_E for Y_i
  constrRiE <- function(rho01, rho2, Q, vars){
    rho0k <- diag(rho01)
    SigmaE_Matrix <- diag((1 - rho0k)*vars)
    for(row in 1:Q) {
      for(col in 1:Q) {
        if(row != col){
          SigmaE_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*(rho2[row,col] - rho01[row,col])
        }
      }
    }
    return(SigmaE_Matrix)
  }

  # Dependent Function 3: Constructs covariance matrix Sigma_phi for Y_i
  constrRiP <- function(rho01, Q, vars) {
    rho0k <- diag(rho01)
    SigmaP_Matrix <- diag(rho0k*vars)
    for(row in 1:Q) {
      for(col in 1:Q) {
        if(row != col){
          SigmaP_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*rho01[row,col]
        }
      }
    }
    return(SigmaP_Matrix)
  }

  # Defining necessary parameters based on input values
  betas <- c(beta1, beta2) # Beta vector
  r_alt <- 1/(r + 1)
  Q <- 2 # Number of outcomes, could extend this to more than 2 in the future
  K_total <- ceiling(K/r_alt) # Total number of clusters
  vars <- c(varY1, varY2) # Vector of variances
  rho01_mat <- matrix(c(rho01, rho1, rho1, rho02), 2, 2)
  rho2_mat <- matrix(c(1, rho2, rho2, 1), 2, 2)
  clus_var <- 0 # cluster variation, placeholder for future extensions

  # Start calculating power
  sigmaz.square <- r_alt*(1 - r_alt) # Variance of treatment assignment
  omega <- calCovbetas(vars, rho01_mat, rho2_mat, clus_var, sigmaz.square, m, Q)
  tau <- K_total*t(betas) %*% solve(omega) %*% betas

  if(dist == "Chi2"){ # Using Chi2
    cv <- qchisq(1 - alpha, df = 2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
    power <- round(1 - pchisq(cv, df = 2, ncp = tau, lower.tail = TRUE), 4)
  } else if(dist == "F"){ # Using F
    Fscore <- qf(1 - alpha, df1 = Q, df2 = K_total - 2*Q, ncp = 0,
                 lower.tail = TRUE, log.p = FALSE)
    power <- round(1 - pf(Fscore, df1 = Q, df2 = K_total - 2*Q, tau,
                          lower.tail = TRUE, log.p = FALSE), 4)
  } else{
    stop("Please choose a valid input parameter for 'dist', either 'Chi2' for Chi-Square or 'F' for F-distribution.")
  }

  return(power)
} # End calc_pwr_disj_2dftest()








#' Calculate required number of clusters per treatment group for a cluster-randomized trial with co-primary endpoints using a disjunctive 2-DF test approach.
#'
#' @description
#' Allows user to calculate the number of clusters per treatment arm of a cluster-randomized trial with two co-primary outcomes given a set of study design input values, including the statistical power, and cluster size. Uses the disjunctive 2-DF test approach. Code is adapted from "calSampleSize_omnibus()" from https://github.com/siyunyang/coprimary_CRT.
#'
#' @param dist Specification of which distribution to base calculation on, either 'Chi2' for Chi-Squared or 'F' for F-Distribution.
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
#' calc_K_disj_2dftest(power = 0.8, m = 300, alpha = 0.05, beta1 = 0.1, beta2 = 0.1, varY1 = 0.23, varY2 = 0.25, rho01 = 0.025, rho02 = 0.025, rho1 = 0.01, rho2  = 0.05)
#' @export
calc_K_disj_2dftest <- function(dist = "Chi2",# Distribution to base calculation from
                                power,        # Desired statistical power
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
    stop("All input parameters must be numeric values (with the exception of 'dist').")
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

  # Define small dependent functions -------------------------------------------
  # Dependent Function 1: Calculates covariance between betas
  calCovbetas <- function(vars, rho01, rho2, cv, sigmaz.square, m, Q){
    sigmaE <- constrRiE(rho01, rho2, Q, vars)
    sigmaP <- constrRiP(rho01, Q, vars)
    tmp <- solve(diag(1,Q) - cv^2*(m*sigmaP %*% solve(sigmaE + m*sigmaP) %*% sigmaE %*% solve(sigmaE + m*sigmaP) ))
    covMatrix <- 1/(m*sigmaz.square)*(sigmaE + m*sigmaP)%*%tmp
    covMatrix <- (covMatrix + t(covMatrix))/2  # symmerize the off-diagonal
    return(covMatrix)
  }

  # Dependent Function 2: Constructs covariance matrix Sigma_E for Y_i
  constrRiE <- function(rho01, rho2, Q, vars){
    rho0k <- diag(rho01)
    SigmaE_Matrix <- diag((1 - rho0k)*vars)
    for(row in 1:Q) {
      for(col in 1:Q) {
        if(row != col){
          SigmaE_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*(rho2[row,col] - rho01[row,col])
        }
      }
    }
    return(SigmaE_Matrix)
  }

  # Dependent Function 3: Constructs covariance matrix Sigma_phi for Y_i
  constrRiP <- function(rho01, Q, vars) {
    rho0k <- diag(rho01)
    SigmaP_Matrix <- diag(rho0k*vars)
    for(row in 1:Q) {
      for(col in 1:Q) {
        if(row != col){
          SigmaP_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*rho01[row,col]
        }
      }
    }
    return(SigmaP_Matrix)
  }

  # Defining necessary parameters based on input values
  betas <- c(beta1, beta2) # Beta vector
  r_alt <- 1/(r + 1)
  Q <- 2 # Number of outcomes, could extend this to more than 2 in the future
  vars <- c(varY1, varY2) # Vector of variances
  rho01_mat <- matrix(c(rho01, rho1, rho1, rho02), 2, 2)
  rho2_mat <- matrix(c(1, rho2, rho2, 1), 2, 2)
  clus_var <- 0 # cluster variation, placeholder for future extensions
  pred.power <- 0 # Initializing predictive power
  sigmaz.square <- r_alt*(1-r_alt) # Variance of treatment assignment

  if(dist == "Chi2"){ # For Chi-Square
    # Initialize small number of clusters
    K_total <- 1
    while(pred.power < power){ # Iterate while the current predictive power is lower than desired power
      K_total <- K_total + 1
      omega <- calCovbetas(vars, rho01_mat, rho2_mat, clus_var, sigmaz.square, m, Q)
      tau <- K_total*t(betas) %*% solve(omega) %*% betas
      cv <- qchisq(1 - alpha, df = 2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
      pred.power <- round(1 - pchisq(cv, df = 2, ncp = tau, lower.tail = TRUE), 4)
    }
  } else if(dist == "F"){ # For F-Distribution
    # Initialize small number of clusters
    K_total <- 2*Q
    while(pred.power < power){ # Iterate while the current predictive power is lower than desired power
      K_total <- K_total + 1
      omega <- calCovbetas(vars, rho01_mat, rho2_mat, clus_var, sigmaz.square, m, Q)
      tau <- K_total*t(betas) %*% solve(omega) %*% betas
      Fscore <- qf(1 - alpha, df1 = Q, df2 = K_total - 2*Q, ncp = 0,
                   lower.tail = TRUE, log.p = FALSE)
      pred.power <- round(1 - pf(Fscore, df1 = Q, df2 = K_total - 2*Q, tau,
                                 lower.tail = TRUE, log.p = FALSE), 4)
    }
  } else{
    stop("Please choose a valid input parameter for 'dist', either 'Chi2' for Chi-Square or 'F' for F-distribution.")
  }

  if(r == 1){ # Case when treatment allocation ratio is equal
    K <- tibble(`Treatment (K)` = ceiling(K_total/2),
                `Control (K)` = ceiling(K_total/2))
  } else{ # Case for unequal treatment allocation, make sure we round up
    K1 <- ceiling(r_alt*K_total)
    K2 <- ceiling(r*K1)
    K <- tibble(`Treatment (K1)` = K1,
                `Control (K2)` = K2)
  }

  return(K)
} # End calc_K_disj_2dftest()









#' Calculate cluster size for a cluster-randomized trial with co-primary endpoints using a disjunctive 2-DF test approach.
#'
#' @description
#' Allows user to calculate the cluster size of a cluster-randomized trial with two co-primary outcomes given a set of study design input values, including the number of clusters in each trial arm, and statistical power. Uses the disjunctive 2-DF test approach.
#'
#' @param dist Specification of which distribution to base calculation on, either 'Chi2' for Chi-Squared or 'F' for F-Distribution.
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
#' calc_m_disj_2dftest(power = 0.8, K = 15, alpha = 0.05, beta1 = 0.1, beta2 = 0.1, varY1 = 0.23, varY2 = 0.25, rho01 = 0.025, rho02 = 0.025, rho1 = 0.01, rho2  = 0.05)
#' @export
calc_m_disj_2dftest <- function(dist = "Chi2",# Distribution to base calculation from
                                power,        # Desired statistical power
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
    stop("All input parameters must be numeric values  (with the exception of 'dist').")
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

  m <- 0 # Initializing a small cluster size
  pred.power <- 0 # Initializing predictive power
  # Iterate while the current predictive power is lower than the power
  # that we actually want (inputted power)

  while(pred.power < power){
    m <- m + 1
    pred.power <- calc_pwr_disj_2dftest(dist = dist,
                                        K = K,
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
                                        )
    if(m > 100000){
      m <- Inf
      message("Cannot find large enough 'm' to reach study specifications for disjunctive 2-DF test. Please lower power or increase value for 'K'.")
      break
    }
  }
  return(m)
} # End calc_m_disj_2dftest()


