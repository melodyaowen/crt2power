#' Calculate statistical power for a cluster-randomized trial with co-primary endpoints using the single 1-DF combined test approach.
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
#'
#' @description
#' Allows user to calculate the statistical power of a cluster-randomized trial with two co-primary endpoints given a set of study design input values, including the number of clusters in each trial arm, and cluster size. Uses the single 1-DF combined test approach for clustered data and two outcomes.
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
#' calc_pwr_single_1dftest(K = 15, m = 300, alpha = 0.05,
#' beta1 = 0.1, beta2 = 0.1, varY1 = 0.23, varY2 = 0.25,
#' rho01 = 0.025, rho02 = 0.025, rho1 = 0.01, rho2  = 0.05)
#' @export
calc_pwr_single_1dftest <- function(dist = "Chi2",# Distribution to base calculation from
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

  # Calculate test statistic for first and second outcome
  Z1.sq <- (beta1^2)/((((1 + 1/r)*varY1)/(K*m))*(1 + (m - 1)*rho01)) # Z1^2
  Z2.sq <- (beta2^2)/((((1 + 1/r)*varY2)/(K*m))*(1 + (m - 1)*rho02)) # Z2^2

  # Calculate correlation between test statistics
  CorrZ1Z2 <- (rho2 + (m - 1)*rho1)/sqrt((1 + (m - 1)*rho01)*(1 + (m - 1)*rho02))

  # Calculate lambda
  lambda <- ((sqrt(Z1.sq) + sqrt(Z2.sq))^2)/(2*(1 + CorrZ1Z2))

  if(dist == "Chi2"){
    # Calculate critical value
    cv <- qchisq(p = alpha, df = 1, lower.tail = FALSE)
    # Calculate power
    power <- round(1 - pchisq(cv, ncp = lambda, df = 1, lower.tail = TRUE), 4)
  } else if(dist == "F"){
    # Calculate critical value
    Fscore <- qf(alpha, df1 = 1, df2 = K_total - 2*Q, ncp = 0,
                 lower.tail = FALSE, log.p = FALSE)
    # Calculate power
    power <- round(1 - pf(Fscore, ncp = lambda,
                          df1 = 1, df2 = K_total - 2*Q,
                          lower.tail = TRUE, log.p = FALSE), 4)
  } else{
    stop("Please choose a valid input parameter for 'dist', either 'Chi2' for Chi-Square or 'F' for F-distribution.")
  }



  return(power)
} # End calc_pwr_single_1dftest()










#' Calculate required number of clusters per treatment group for a cluster-randomized trial with co-primary endpoints using the single 1-DF combined test approach.
#'
#' @description
#' Allows user to calculate the number of clusters per treatment arm of a cluster-randomized trial with two co-primary endpoints given a set of study design input values, including the statistical power, and cluster size. Uses the single 1-DF combined test approach for clustered data and two outcomes.
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
#' calc_K_single_1dftest(power = 0.8, m = 300, alpha = 0.05,
#' beta1 = 0.1, beta2 = 0.1, varY1 = 0.23, varY2 = 0.25,
#' rho01 = 0.025, rho02 = 0.025, rho1 = 0.01, rho2  = 0.05)
#' @export
calc_K_single_1dftest <- function(dist = "Chi2",# Distribution to base calculation from
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

  if(dist == "Chi2"){ # Using Chi2

    # Calculate correlation between test statistics
    CorrZ1Z2 <- (rho2 + (m - 1)*rho1)/sqrt((1 + (m - 1)*rho01)*(1 + (m - 1)*rho02))

    # Non-centrality parameter for given power and alpha level
    ncp <- calc_ncp_chi2(alpha, power, df = 1)

    # Calculate K
    K1 <- ceiling((2*ncp*(1 + CorrZ1Z2))/
                    (sqrt((beta1^2)/(((1 + 1/r)*varY1/m)*(1+(m-1)*rho01))) +
                       sqrt((beta2^2)/(((1 + 1/r)*varY2/m)*(1+(m-1)*rho02))))^2)

    if(r == 1){
      K <- tibble(`Treatment (K)` = K1,
                  `Control (K)` = K1)
    } else{
      K2 <- ceiling(r*K1)
      K <- tibble(`Treatment (K1)` = K1,
                  `Control (K2)` = K2)
    }

  } else if(dist == "F"){ # Using F

    # Initialize a wide range of values for m
    K_options <- 3:5000

    # Retrieve power for all choices of m
    res_list <- lapply(K_options, function(K) {
      out_vector <- calc_pwr_single_1dftest(
        m     = m,
        dist  = dist,
        K     = K,
        alpha = alpha,
        beta1 = beta1,
        beta2 = beta2,
        varY1 = varY1,
        varY2 = varY2,
        rho01 = rho01,
        rho02 = rho02,
        rho1  = rho1,
        rho2  = rho2,
        r     = r
      )

      # Add a column 'K'.
      out_K <- tibble(`Final Power` = out_vector, K)

      out_K
    })

    # Flatten to make one dataframe
    res_df <- do.call(rbind, res_list) %>%
      dplyr::rename(K1 = K) %>%
      mutate(K2 = ceiling(r*K1)) %>%
      mutate(K_total = K1 + K2)
    K_results <- dplyr::filter(res_df, `Final Power` >= power)

    # If there are no choices of "K" where the desired power is met, throw an error
    if(nrow(K_results) == 0){
      stop("Cannot find large enough 'K' to reach study specifications for single weighted 1-DF method. Please lower power or increase value for 'm'.")
    }

    K_final <- dplyr::filter(K_results, K_total == min(K_total))
    K1_final <- K_final$K1
    K2_final <- K_final$K2

    # K for Method 2
    if(r == 1){ # When treatment allocation is even
      K <- tibble(`Treatment (K)` = K1_final,
                  `Control (K)` = K2_final)
    } else{ # Unequal treatment allocation
      K <- tibble(`Treatment (K1)` = K1_final,
                  `Control (K2)` = K2_final)
    }

  } else{
    stop("Please choose a valid input parameter for 'dist', either 'Chi2' for Chi-Square or 'F' for F-distribution.")
  }

  return(K)
} # End calc_K_single_1dftest()










#' Calculate cluster size for a cluster-randomized trial with co-primary endpoints using the single 1-DF combined test approach.
#'
#' @description
#' Allows user to calculate the cluster size of a cluster-randomized trial with two co-primary endpoints given a set of study design input values, including the number of clusters in each trial arm, and statistical power. Uses the single 1-DF combined test approach for clustered data and two outcomes.
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
#' calc_m_single_1dftest(power = 0.8, K = 15, alpha = 0.05,
#' beta1 = 0.1, beta2 = 0.1, varY1 = 0.23, varY2 = 0.25,
#' rho01 = 0.025, rho02 = 0.025, rho1 = 0.01, rho2  = 0.05)
#' @export
calc_m_single_1dftest <- function(dist = "Chi2",# Distribution to base calculation from
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

  if(dist == "Chi2"){ # Using Chi2

    # Non-centrality parameter for given power and alpha level
    ncp <- calc_ncp_chi2(alpha, power, df = 1)

    # Function to solve for m
    power.eq.m <- function(m, para){ # Function for equation solve for m
      # Calculate test statistic for first and second outcome
      Z1.sq <- (para[1]^2)/(((1 + 1/para[11])*para[4]/(m*para[3]))*(1 + (m - 1)*para[6])) # Z1^2
      Z2.sq <- (para[2]^2)/(((1 + 1/para[11])*para[5]/(m*para[3]))*(1 + (m - 1)*para[7])) # Z1^2
      # Calculate Corr(Z1, Z2)
      CorrZ1Z2 <- (para[9] + (m - 1)*para[8])/
        sqrt((1+(m-1)*para[6])*(1+(m- 1)*para[7]))
      # Non-centrality parameter
      lambda <- ((sqrt(Z1.sq) + sqrt(Z2.sq))^2) / (2*(1 + CorrZ1Z2)) - para[10]
    }

    # Find m using multiroot function
    find.para.m <- multiroot(power.eq.m, start = 10,
                             para = c(beta1, beta2, K, varY1, varY2,
                                      rho01, rho02, rho1, rho2, ncp, r),
                             positive = TRUE)
    m <- ceiling(find.para.m$root)

    if(m < 1){
      stop("Cannot find large enough 'm' to reach study specifications for single weighted 1-DF method. Please lower power or increase value for 'K'.")
    }

  } else if(dist == "F"){ # Using F

    # Initialize a wide range of values for m
    m_options <- 1:10000

    # Retrieve power for all choices of m
    res_list <- lapply(m_options, function(m) {
      out_vector <- calc_pwr_single_1dftest(
        m     = m,
        dist  = dist,
        K     = K,
        alpha = alpha,
        beta1 = beta1,
        beta2 = beta2,
        varY1 = varY1,
        varY2 = varY2,
        rho01 = rho01,
        rho02 = rho02,
        rho1  = rho1,
        rho2  = rho2,
        r     = r
      )

      # Add a column 'm'.
      out_m <- tibble(`Final Power` = out_vector, m)

      out_m
    })

    # Flatten to make one dataframe
    res_df <- do.call(rbind, res_list)

    # Get vector of all "m" values that meet the desired power
    m_results <- dplyr::filter(res_df, `Final Power` >= power)

    # If there are no choices of "m" where the desired power is met, throw an error
    if(nrow(m_results) == 0){
      stop("Cannot find large enough 'm' to reach study specifications for single weighted 1-DF method. Please lower power or increase value for 'K'.")
    }

    m <- min(m_results$m)

  } else{
    stop("Please choose a valid input parameter for 'dist', either 'Chi2' for Chi-Square or 'F' for F-distribution.")
  }

  return(m)
} # End calc_m_single_1dftest()
