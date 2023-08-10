devtools::source_url("https://github.com/siyunyang/coprimary_CRT/blob/main/powerSampleCal_varCluster_ttest.R?raw=TRUE")

calc_pwr_conj_test <- function(K,            # Number of clusters in each arm
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

  power <-  calPower_ttestIU(betas = c(beta1, beta2),
                             K = 2, # In this function, K = # of outcomes
                             m = m,
                             deltas = c(0, 0),
                             vars = c(varY1, varY2),
                             rho01 = matrix(c(rho01, rho1,
                                              rho1, rho02),
                                            2, 2),
                             rho2 = matrix(c(1, rho2,
                                             rho2, 1),
                                           2, 2),
                             r = 0.5,
                             N = 2*K, # In this function, N = total clusters
                             alpha = alpha
                             )

  return(round(power, 4))

} # End calc_pwr_conj_test()



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
                             rho2          # Intra-subject between-endpoint ICC
                             ){

  K <- calSampleSize_ttestIU(betas = c(beta1, beta2),
                             m = m,
                             power = power,
                             deltas = c(0, 0),
                             vars = c(varY1, varY2),
                             rho01 = matrix(c(rho01, rho1,
                                              rho1, rho02),
                                            2, 2),
                             rho2 = matrix(c(1, rho2,
                                             rho2, 1),
                                           2, 2),
                             r = 0.5,
                             K = 2, # In this function, K = # of outcomes
                             alpha = alpha
                             )/2 # Divide by 2 because function returns
                                 # total number of clusters, which is 2K

  return(ceiling(K))
} # End calc_K_conj_test()



calc_m_conj_test <- function(power,        # Desired statistical power
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

  # Initialize m and predictive power
  m <- 0
  pred.power <- 0
  while(pred.power < power){
    m <- m + 1
    pred.power <- calPower_ttestIU(
                                   betas = c(beta1, beta2),
                                   K = 2, # In this function, K = # of outcomes
                                   m = m,
                                   deltas = c(0, 0),
                                   vars = c(varY1, varY2),
                                   rho01 = matrix(c(rho01, rho1,
                                                    rho1, rho02),
                                                  2, 2),
                                   rho2 = matrix(c(1, rho2,
                                                   rho2, 1),
                                                 2, 2),
                                   r = 0.5,
                                   N = 2*K, # In this function, N = total clusters
                                   alpha = alpha
                                   )
    if(m > 1000000){
      m <- Inf
      break
    }
  }
  return(ceiling(m))

} # End calc_m_conj_test

