#' Calculate statistical power for a cluster randomized trial with co-primary endpoints using three common p-value adjustment methods
#'
#' @description
#' Allows user to calculate the statistical power of a hybrid type 2 cluster randomized trial given a set of study design input values, including the number of clusters in each trial arm, and cluster size. Uses three common p-value adjustment methods.
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
#' @param rho2 Correlation between the first and second outcomes for the same individual; numeric.
#' @returns A data frame of numerical values.
#' @examples
#' calc_pwr_pval_adj(K = 15, m = 300, alpha = 0.05, beta1 = 0.1, beta2 = 0.1, varY1 = 0.23, varY2 = 0.25, rho01 = 0.025, rho02 = 0.025, rho2  = 0.05)
calc_pwr_pval_adj <- function(K,            # Number of clusters in each arm
                              m,            # Individuals per cluster
                              alpha = 0.05, # Significance level
                              beta1,        # Effect for outcome 1
                              beta2,        # Effect for outcome 2
                              varY1,        # Variance for outcome 1
                              varY2,        # Variance for outcome 2
                              rho01,        # ICC for outcome 1
                              rho02,        # ICC for outcome 2
                              rho2          # Intra-subject between-endpoint ICC
                              ){

  # Adjusted p-values for three methods based on inputted alpha-level
  alpha_B <- alpha/2 # Bonferroni
  alpha_S <- 1 - (1 - alpha)^(1/2) # Sidak
  alpha_D <- 1 - (1 - alpha)^(1/(2^(1 - rho2))) # D/AP

  # Critical values for three p-value adjustment methods
  cv_B <- qchisq(1 - alpha_B, df = 1, ncp = 0) # Bonferroni
  cv_S <- qchisq(1 - alpha_S, df = 1, ncp = 0) # Sidak
  cv_D <- qchisq(1 - alpha_D, df = 1, ncp = 0) # D/AP

  # Power for each p-value adjustment method
  lambda1 <- (beta1^2)/(2*(varY1/(K*m))*(1 + (m-1)*rho01))
  lambda2 <- (beta2^2)/(2*(varY2/(K*m))*(1 + (m-1)*rho02))

  pwr_bonf  <- round(c(1 - pchisq(cv_B, 1, ncp = lambda1, lower.tail = TRUE),
                       1 - pchisq(cv_B, 1, ncp = lambda2, lower.tail = TRUE)), 4)
  pwr_sidak <- round(c(1 - pchisq(cv_S, 1, ncp = lambda1, lower.tail = TRUE),
                       1 - pchisq(cv_S, 1, ncp = lambda2, lower.tail = TRUE)), 4)
  pwr_dap   <- round(c(1 - pchisq(cv_D, 1, ncp = lambda1, lower.tail = TRUE),
                       1 - pchisq(cv_D, 1, ncp = lambda2, lower.tail = TRUE)), 4)

  # Final dataframe of power calculations by method and outcome
  pwr_dat <- tibble(`P-Value Adjustment Type` = c("Bonferroni", "Sidak", "D/AP"),
                    `Power (Y1)` = c(pwr_bonf[1], pwr_sidak[1], pwr_dap[1]),
                    `Power (Y2)` = c(pwr_bonf[2], pwr_sidak[2], pwr_dap[2]),
                    `Final Power` = c(min(pwr_bonf), min(pwr_sidak), min(pwr_dap)))

  return(pwr_dat) # Return dataframe of final power calculations
} # End calc_pwr_pval_adj()


#' Calculate required number of clusters per treatment group for a cluster randomized trial with co-primary endpoints using three common p-value adjustment methods
#'
#' @description
#' Allows user to calculate the number of clusters per treatment arm of a hybrid type 2 cluster randomized trial given a set of study design input values, including the statistical power, and cluster size. Uses three common p-value adjustment methods.
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
#' @param rho2 Correlation between the first and second outcomes for the same individual; numeric.
#' @returns A whole numer.
#' @examples
#' calc_K_pval_adj(power = 0.8, m = 300, alpha = 0.05, beta1 = 0.1, beta2 = 0.1, varY1 = 0.23, varY2 = 0.25, rho01 = 0.025, rho02 = 0.025, rho2  = 0.05)
calc_K_pval_adj <- function(power,        # Desired statistical power
                            m,            # Individuals per cluster
                            alpha = 0.05, # Significance level
                            beta1,        # Effect for outcome 1
                            beta2,        # Effect for outcome 2
                            varY1,        # Variance for outcome 1
                            varY2,        # Variance for outcome 2
                            rho01,        # ICC for outcome 1
                            rho02,        # ICC for outcome 2
                            rho2          # Intra-subject between-endpoint ICC
                            ){

  # Adjusted p-values for three methods based on inputted alpha-level
  alpha_B <- alpha/2 # Bonferroni
  alpha_S <- 1 - (1 - alpha)^(1/2) # Sidak
  alpha_D <- 1 - (1 - alpha)^(1/(2^(1 - rho2))) # D/AP

  # Lambda values for calculating K for three methods
  lambda_B <- calc_ncp_chi2(alpha_B, power, df = 1)
  lambda_S <- calc_ncp_chi2(alpha_S, power, df = 1)
  lambda_D <- calc_ncp_chi2(alpha_D, power, df = 1)

  # Power for each p-value adjustment method
  K_bonf  <- c((2*lambda_B*varY1*(1 + (m - 1)*rho01))/(m*(beta1^2)),
                       (2*lambda_B*varY2*(1 + (m - 1)*rho02))/(m*(beta2^2)))
  K_sidak <- c((2*lambda_S*varY1*(1 + (m - 1)*rho01))/(m*(beta1^2)),
                       (2*lambda_S*varY2*(1 + (m - 1)*rho02))/(m*(beta2^2)))
  K_dap   <- c((2*lambda_D*varY1*(1 + (m - 1)*rho01))/(m*(beta1^2)),
                       (2*lambda_D*varY2*(1 + (m - 1)*rho02))/(m*(beta2^2)))

  # Final dataframe of K calculations by method and outcome
  K_dat <- data.frame(type = c("Bonferroni", "Sidak", "D/AP"),
                      K_Y1 = ceiling(c(K_bonf[1], K_sidak[1], K_dap[1])),
                      K_Y2 = ceiling(c(K_bonf[2], K_sidak[2], K_dap[2])),
                      final = c(ceiling(max(K_bonf)),
                                ceiling(max(K_sidak)),
                                ceiling(max(K_dap))))

  return(K_dat) # Return dataframe of final K calculations
} # End calc_K_pval_adj()


#' Calculate cluster size for a cluster randomized trial with co-primary endpoints using three common p-value adjustment methods
#'
#'#' @description
#' Allows user to calculate the cluster size of a hybrid type 2 cluster randomized trial given a set of study design input values, including the number of clusters in each trial arm, and statistical power. Uses three common p-value adjustment methods.
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
#' @param rho2 Correlation between the first and second outcomes for the same individual; numeric.
#' @returns A whole number.
#' @examples
#' calc_m_pval_adj(power = 0.8, K = 15, alpha = 0.05, beta1 = 0.1, beta2 = 0.1, varY1 = 0.23, varY2 = 0.25, rho01 = 0.025, rho02 = 0.025, rho2  = 0.05)
calc_m_pval_adj <- function(power,        # Desired statistical power
                            K,            # Number of clusters in each arm
                            alpha = 0.05, # Significance level
                            beta1,        # Effect for outcome 1
                            beta2,        # Effect for outcome 2
                            varY1,        # Variance for outcome 1
                            varY2,        # Variance for outcome 2
                            rho01,        # ICC for outcome 1
                            rho02,        # ICC for outcome 2
                            rho2          # Intra-subject between-endpoint ICC
                            ){

  # Adjusted p-values for three methods based on inputted alpha-level
  alpha_B <- alpha/2 # Bonferroni
  alpha_S <- 1 - (1 - alpha)^(1/2) # Sidak
  alpha_D <- 1 - (1 - alpha)^(1/(2^(1 - rho2))) # D/AP

  # Lambda values for calculating m for three methods
  lambda_B <- calc_ncp_chi2(alpha_B, power, df = 1)
  lambda_S <- calc_ncp_chi2(alpha_S, power, df = 1)
  lambda_D <- calc_ncp_chi2(alpha_D, power, df = 1)

  # Power for each p-value adjustment method
  m_bonf  <- c((2*lambda_B*varY1*(1 - rho01))/((beta1^2)*K - 2*lambda_B*varY1*rho01),
               (2*lambda_B*varY2*(1 - rho02))/((beta2^2)*K - 2*lambda_B*varY2*rho02))
  m_sidak <- c((2*lambda_S*varY1*(1 - rho01))/((beta1^2)*K - 2*lambda_S*varY1*rho01),
               (2*lambda_S*varY2*(1 - rho02))/((beta2^2)*K - 2*lambda_S*varY2*rho02))
  m_dap   <- c((2*lambda_D*varY1*(1 - rho01))/((beta1^2)*K - 2*lambda_D*varY1*rho01),
               (2*lambda_D*varY2*(1 - rho02))/((beta2^2)*K - 2*lambda_D*varY2*rho02))

  # Final dataframe of m calculations by method and outcome
  m_dat <- data.frame(type = c("Bonferroni", "Sidak", "D/AP"),
                      m_Y1 = ceiling(c(m_bonf[1], m_sidak[1], m_dap[1])),
                      m_Y2 = ceiling(c(m_bonf[2], m_sidak[2], m_dap[2])),
                      final = c(ceiling(max(m_bonf)), ceiling(max(m_sidak)),
                                ceiling(max(m_dap))))

  return(m_dat) # Return dataframe of final m calculations
} # End calc_m_pval_adj()

