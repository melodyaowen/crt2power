#' Find study design output specifications based on all five CRT co-primary design methods.
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
#' Allows user to calculate either statistical power, number of clusters per treatment group (K), or cluster size (m), given a set of input values for all five study design approaches.
#'
#' @param output Parameter to calculate, either "power", "K", or "m"; character.
#' @param power Desired statistical power; numeric.
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
#' @param r Treatment allocation ratio - K2 = rK1 where K1 is number of clusters in experimental group; numeric.
#' @returns A data frame of numerical values.
#' @examples
#' run_crt2_design(output = "power", K = 15, m = 300, alpha = 0.05,
#' beta1 = 0.1, beta2 = 0.1, varY1 = 0.23, varY2 = 0.25,
#' rho01 = 0.025, rho02 = 0.025, rho1 = 0.01, rho2  = 0.05)
#' @export
run_crt2_design <- function(output,       # Parameter to calculate
                            power = NA,   # Desired statistical power
                            K = NA,       # Number of clusters in each arm
                            m = NA,       # Individuals per cluster
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

  # Checks to make sure inputs are valid
  if(output == "power"){ # When output is power
    if(!is.na(power)){
      stop("'power' variable cannot be defined if desired study design output is 'power'.")
    }
    if(is.na(K)){
      stop("Must define 'K' in order to calculate power.")
    } else if(!is.numeric(K) | K < 1){
      stop("'K' must be a positive whole number.")
    } else if(K != round(K)) {
      stop("'K' must be a positive whole number.")
    }
    if(is.na(m)){
      stop("Must define 'm' in order to calculate power.")
    } else if(!is.numeric(m) | m < 1){
      stop("'m' must be a positive whole number.")
    } else if(m != round(m)){
      stop("'m' must be a positive whole number.")
    }
  } else if(output == "K"){ # When output is K
    if(!is.na(K)){
      stop("'K' variable cannot be defined if desired study design output is 'K'.")
    }
    if(is.na(power)){
      stop("Must define 'power' in order to calculate K.")
    } else if(!is.numeric(power) | power > 1 | power < 0){
      stop("'power' must be a number between 0 and 1.")
    }
    if(is.na(m)){
      stop("Must define 'm' in order to calculate K.")
    } else if(!is.numeric(m) | m < 1){
      stop("'m' must be a positive whole number.")
    } else if(m != round(m)){
      stop("'m' must be a positive whole number.")
    }
  } else if(output == "m"){ # When output is m
    if(!is.na(m)){
      stop("'m' variable cannot be defined if desired study design output is 'm'.")
    }
    if(is.na(K)){
      stop("Must define 'K' in order to calculate m.")
    } else if(!is.numeric(K) | K < 1){
      stop("'K' must be a positive whole number.")
    } else if(K != round(K)){
      stop("'K' must be a positive whole number.")
    }
    if(is.na(power)){
      stop("Must define 'power' in order to calculate m.")
    } else if(!is.numeric(power) | power > 1 | power < 0){
      stop("'power' must be a number between 0 and 1.")
    }
  } else {
    stop("'output' parameter must either be 'power', 'K', or 'm'.")
  }

  # Check to make sure variables are numeric
  if(!is.numeric(c(alpha, beta1, beta2, varY1, varY2, rho01, rho02, rho1, rho2))){
    stop("All input parameters must be numeric (with the exception of 'output').")
  }

  # When desired output is power
  if(output == "power"){
    # Method 1: P-Value Adjustments
    out1_Chi2 <- calc_pwr_pval_adj(K = K, m = m, alpha = alpha,
                                   beta1 = beta1, beta2 = beta2,
                                   varY1 = varY1, varY2 = varY2,
                                   rho01 = rho01, rho02 = rho02,
                                   rho2  = rho2, r = r, dist = "Chi2")
    out1_F <- calc_pwr_pval_adj(K = K, m = m, alpha = alpha,
                                beta1 = beta1, beta2 = beta2,
                                varY1 = varY1, varY2 = varY2,
                                rho01 = rho01, rho02 = rho02,
                                rho2  = rho2, r = r, dist = "F")

    # Method 2: Combined Outcome
    out2_Chi2 <- calc_pwr_comb_outcome(K = K, m = m, alpha = alpha,
                                       beta1 = beta1, beta2 = beta2,
                                       varY1 = varY1, varY2 = varY2,
                                       rho01 = rho01, rho02 = rho02,
                                       rho1  = rho1, rho2 = rho2,
                                       r = r, dist = "Chi2")
    out2_F <- calc_pwr_comb_outcome(K = K, m = m, alpha = alpha,
                                    beta1 = beta1, beta2 = beta2,
                                    varY1 = varY1, varY2 = varY2,
                                    rho01 = rho01, rho02 = rho02,
                                    rho1  = rho1, rho2 = rho2,
                                    r = r, dist = "F")

    # Method 3: Single 1-DF Test
    out3_Chi2 <- calc_pwr_single_1dftest(K = K, m = m, alpha = alpha,
                                         beta1 = beta1, beta2 = beta2,
                                         varY1 = varY1, varY2 = varY2,
                                         rho01 = rho01, rho02 = rho02,
                                         rho1  = rho1, rho2  = rho2,
                                         r = r, dist = "Chi2")
    out3_F <- calc_pwr_single_1dftest(K = K, m = m, alpha = alpha,
                                      beta1 = beta1, beta2 = beta2,
                                      varY1 = varY1, varY2 = varY2,
                                      rho01 = rho01, rho02 = rho02,
                                      rho1  = rho1, rho2  = rho2,
                                      r = r, dist = "F")

    # Method 4: Disjunctive 2-DF Test
    out4_Chi2 <- calc_pwr_disj_2dftest(K = K, m = m, alpha = alpha,
                                       beta1 = beta1, beta2 = beta2,
                                       varY1 = varY1, varY2 = varY2,
                                       rho01 = rho01, rho02 = rho02,
                                       rho1  = rho1, rho2  = rho2,
                                       r = r, dist = "Chi2")
    out4_F <- calc_pwr_disj_2dftest(K = K, m = m, alpha = alpha,
                                    beta1 = beta1, beta2 = beta2,
                                    varY1 = varY1, varY2 = varY2,
                                    rho01 = rho01, rho02 = rho02,
                                    rho1  = rho1, rho2  = rho2,
                                    r = r, dist = "F")

    # Method 5: Conjunctive IU Test
    out5_T_1sided <- calc_pwr_conj_test(K = K, m = m, alpha = alpha,
                                        beta1 = beta1, beta2 = beta2,
                                        varY1 = varY1, varY2 = varY2,
                                        rho01 = rho01, rho02 = rho02,
                                        rho1  = rho1, rho2  = rho2,
                                        r = r, dist = "T", two_sided = FALSE)
    out5_T_2sided <- calc_pwr_conj_test(K = K, m = m, alpha = alpha,
                                        beta1 = beta1, beta2 = beta2,
                                        varY1 = varY1, varY2 = varY2,
                                        rho01 = rho01, rho02 = rho02,
                                        rho1  = rho1, rho2  = rho2,
                                        r = r, dist = "T", two_sided = TRUE)
    out5_MVN_1sided <- calc_pwr_conj_test(K = K, m = m, alpha = alpha,
                                          beta1 = beta1, beta2 = beta2,
                                          varY1 = varY1, varY2 = varY2,
                                          rho01 = rho01, rho02 = rho02,
                                          rho1  = rho1, rho2  = rho2,
                                          r = r, dist = "MVN", two_sided = FALSE)
    out5_MVN_2sided <- calc_pwr_conj_test(K = K, m = m, alpha = alpha,
                                          beta1 = beta1, beta2 = beta2,
                                          varY1 = varY1, varY2 = varY2,
                                          rho01 = rho01, rho02 = rho02,
                                          rho1  = rho1, rho2  = rho2,
                                          r = r, dist = "MVN", two_sided = TRUE)
  }

  # When desired output is K
  if(output == "K"){
    # Method 1: P-Value Adjustments
    out1_Chi2 <- calc_K_pval_adj(power = power, m = m, alpha = alpha,
                                 beta1 = beta1, beta2 = beta2,
                                 varY1 = varY1, varY2 = varY2,
                                 rho01 = rho01, rho02 = rho02,
                                 rho2  = rho2, r = r, dist = "Chi2")
    out1_F <- calc_K_pval_adj(power = power, m = m, alpha = alpha,
                              beta1 = beta1, beta2 = beta2,
                              varY1 = varY1, varY2 = varY2,
                              rho01 = rho01, rho02 = rho02,
                              rho2  = rho2, r = r, dist = "F")

    # Method 2: Combined Outcome
    out2_Chi2 <- calc_K_comb_outcome(power = power, m = m, alpha = alpha,
                                     beta1 = beta1, beta2 = beta2,
                                     varY1 = varY1, varY2 = varY2,
                                     rho01 = rho01, rho02 = rho02,
                                     rho1  = rho1, rho2 = rho2,
                                     r = r, dist = "Chi2")
    out2_F <- calc_K_comb_outcome(power = power, m = m, alpha = alpha,
                                  beta1 = beta1, beta2 = beta2,
                                  varY1 = varY1, varY2 = varY2,
                                  rho01 = rho01, rho02 = rho02,
                                  rho1  = rho1, rho2 = rho2,
                                  r = r, dist = "F")

    # Method 3: Single 1-DF Test
    out3_Chi2 <- calc_K_single_1dftest(power = power, m = m, alpha = alpha,
                                       beta1 = beta1, beta2 = beta2,
                                       varY1 = varY1, varY2 = varY2,
                                       rho01 = rho01, rho02 = rho02,
                                       rho1  = rho1, rho2  = rho2,
                                       r = r, dist = "Chi2")
    out3_F <- calc_K_single_1dftest(power = power, m = m, alpha = alpha,
                                    beta1 = beta1, beta2 = beta2,
                                    varY1 = varY1, varY2 = varY2,
                                    rho01 = rho01, rho02 = rho02,
                                    rho1  = rho1, rho2  = rho2,
                                    r = r, dist = "F")

    # Method 4: Disjunctive 2-DF Test
    out4_Chi2 <- calc_K_disj_2dftest(power = power, m = m, alpha = alpha,
                                     beta1 = beta1, beta2 = beta2,
                                     varY1 = varY1, varY2 = varY2,
                                     rho01 = rho01, rho02 = rho02,
                                     rho1  = rho1, rho2  = rho2,
                                     r = r, dist = "Chi2")
    out4_F <- calc_K_disj_2dftest(power = power, m = m, alpha = alpha,
                                  beta1 = beta1, beta2 = beta2,
                                  varY1 = varY1, varY2 = varY2,
                                  rho01 = rho01, rho02 = rho02,
                                  rho1  = rho1, rho2  = rho2,
                                  r = r, dist = "F")

    # Method 5: Conjunctive IU Test
    out5_T_1sided <- calc_K_conj_test(power = power, m = m, alpha = alpha,
                                      beta1 = beta1, beta2 = beta2,
                                      varY1 = varY1, varY2 = varY2,
                                      rho01 = rho01, rho02 = rho02,
                                      rho1  = rho1, rho2  = rho2,
                                      r = r, dist = "T", two_sided = FALSE)
    out5_T_2sided <- calc_K_conj_test(power = power, m = m, alpha = alpha,
                                      beta1 = beta1, beta2 = beta2,
                                      varY1 = varY1, varY2 = varY2,
                                      rho01 = rho01, rho02 = rho02,
                                      rho1  = rho1, rho2  = rho2,
                                      r = r, dist = "T", two_sided = TRUE)
    out5_MVN_1sided <- calc_K_conj_test(power = power, m = m, alpha = alpha,
                                        beta1 = beta1, beta2 = beta2,
                                        varY1 = varY1, varY2 = varY2,
                                        rho01 = rho01, rho02 = rho02,
                                        rho1  = rho1, rho2  = rho2,
                                        r = r, dist = "MVN", two_sided = FALSE)
    out5_MVN_2sided <- calc_K_conj_test(power = power, m = m, alpha = alpha,
                                        beta1 = beta1, beta2 = beta2,
                                        varY1 = varY1, varY2 = varY2,
                                        rho01 = rho01, rho02 = rho02,
                                        rho1  = rho1, rho2  = rho2,
                                        r = r, dist = "MVN", two_sided = TRUE)
  }

  # When desired output is m
  if(output == "m"){
    # Method 1: P-Value Adjustments
    out1_Chi2 <- calc_m_pval_adj(K = K, power = power, alpha = alpha,
                                 beta1 = beta1, beta2 = beta2,
                                 varY1 = varY1, varY2 = varY2,
                                 rho01 = rho01, rho02 = rho02,
                                 rho2  = rho2, r = r, dist = "Chi2")
    out1_F <- calc_m_pval_adj(K = K, power = power, alpha = alpha,
                              beta1 = beta1, beta2 = beta2,
                              varY1 = varY1, varY2 = varY2,
                              rho01 = rho01, rho02 = rho02,
                              rho2  = rho2, r = r, dist = "F")

    # Method 2: Combined Outcome
    out2_Chi2 <- calc_m_comb_outcome(K = K, power = power, alpha = alpha,
                                     beta1 = beta1, beta2 = beta2,
                                     varY1 = varY1, varY2 = varY2,
                                     rho01 = rho01, rho02 = rho02,
                                     rho1  = rho1, rho2 = rho2,
                                     r = r, dist = "Chi2")
    out2_F <- calc_m_comb_outcome(K = K, power = power, alpha = alpha,
                                  beta1 = beta1, beta2 = beta2,
                                  varY1 = varY1, varY2 = varY2,
                                  rho01 = rho01, rho02 = rho02,
                                  rho1  = rho1, rho2 = rho2,
                                  r = r, dist = "F")

    # Method 3: Single 1-DF Test
    out3_Chi2 <- calc_m_single_1dftest(K = K, power = power, alpha = alpha,
                                       beta1 = beta1, beta2 = beta2,
                                       varY1 = varY1, varY2 = varY2,
                                       rho01 = rho01, rho02 = rho02,
                                       rho1  = rho1, rho2  = rho2,
                                       r = r, dist = "Chi2")
    out3_F <- calc_m_single_1dftest(K = K, power = power, alpha = alpha,
                                    beta1 = beta1, beta2 = beta2,
                                    varY1 = varY1, varY2 = varY2,
                                    rho01 = rho01, rho02 = rho02,
                                    rho1  = rho1, rho2  = rho2,
                                    r = r, dist = "F")

    # Method 4: Disjunctive 2-DF Test
    out4_Chi2 <- calc_m_disj_2dftest(K = K, power = power, alpha = alpha,
                                     beta1 = beta1, beta2 = beta2,
                                     varY1 = varY1, varY2 = varY2,
                                     rho01 = rho01, rho02 = rho02,
                                     rho1  = rho1, rho2  = rho2,
                                     r = r, dist = "Chi2")
    out4_F <- calc_m_disj_2dftest(K = K, power = power, alpha = alpha,
                                  beta1 = beta1, beta2 = beta2,
                                  varY1 = varY1, varY2 = varY2,
                                  rho01 = rho01, rho02 = rho02,
                                  rho1  = rho1, rho2  = rho2,
                                  r = r, dist = "F")

    # Method 5: Conjunctive IU Test
    out5_T_1sided <- calc_m_conj_test(K = K, power = power, alpha = alpha,
                                      beta1 = beta1, beta2 = beta2,
                                      varY1 = varY1, varY2 = varY2,
                                      rho01 = rho01, rho02 = rho02,
                                      rho1  = rho1, rho2  = rho2,
                                      r = r, dist = "T", two_sided = FALSE)
    out5_T_2sided <- calc_m_conj_test(K = K, power = power, alpha = alpha,
                                      beta1 = beta1, beta2 = beta2,
                                      varY1 = varY1, varY2 = varY2,
                                      rho01 = rho01, rho02 = rho02,
                                      rho1  = rho1, rho2  = rho2,
                                      r = r, dist = "T", two_sided = TRUE)
    out5_MVN_1sided <- calc_m_conj_test(K = K, power = power, alpha = alpha,
                                        beta1 = beta1, beta2 = beta2,
                                        varY1 = varY1, varY2 = varY2,
                                        rho01 = rho01, rho02 = rho02,
                                        rho1  = rho1, rho2  = rho2,
                                        r = r, dist = "MVN", two_sided = FALSE)
    out5_MVN_2sided <- calc_m_conj_test(K = K, power = power, alpha = alpha,
                                        beta1 = beta1, beta2 = beta2,
                                        varY1 = varY1, varY2 = varY2,
                                        rho01 = rho01, rho02 = rho02,
                                        rho1  = rho1, rho2  = rho2,
                                        r = r, dist = "MVN", two_sided = TRUE)
  }



  # Rename "Output" column to be named what the output parameter is
  if(output == "power"){
    outputTable <- tibble(`Design Method` = c("1. P-Value Adjustments",
                                              "a. Bonferroni",
                                              "b. Sidak",
                                              "c. D/AP",
                                              "2. Combined Outcomes",
                                              "3. Single Weighted 1-DF Combined Test",
                                              "4. Disjunctive 2-DF Test",
                                              "5. Conjunctive IU Test (MVN and t-dist.)",
                                              "a. 1-Sided Hypothesis",
                                              "b. 2-Sided Hypothesis"),
                          `Power (Chi2-distribution)` = c(NA,
                                                          pull(out1_Chi2[, ncol(out1_Chi2)]),
                                                          out2_Chi2,
                                                          out3_Chi2,
                                                          out4_Chi2,
                                                          NA,
                                                          out5_MVN_1sided,
                                                          out5_MVN_2sided),
                          `Power (F-distribution)` = c(NA,
                                                       pull(out1_F[, ncol(out1_F)]),
                                                       out2_F,
                                                       out3_F,
                                                       out4_F,
                                                       NA,
                                                       out5_T_1sided,
                                                       out5_T_2sided))
  } else if(output == "K"){
    outputTable <- tibble(`Design Method` = c("1. P-Value Adjustments",
                                              "a. Bonferroni",
                                              "b. Sidak",
                                              "c. D/AP",
                                              "2. Combined Outcomes",
                                              "3. Single Weighted 1-DF Combined Test",
                                              "4. Disjunctive 2-DF Test",
                                              "5. Conjunctive IU Test (MVN and t-dist.)",
                                              "a. 1-Sided Hypothesis",
                                              "b. 2-Sided Hypothesis"),
                          `Treatment K1 (Chi2-Distribution)` = c(NA,
                                   pull(dplyr::select(out1_Chi2, contains("Final Treatment"))),
                                   pull(dplyr::select(out2_Chi2, contains("Treatment"))),
                                   pull(dplyr::select(out3_Chi2, contains("Treatment"))),
                                   pull(dplyr::select(out4_Chi2, contains("Treatment"))),
                                   NA,
                                   pull(dplyr::select(out5_MVN_1sided, contains("Treatment"))),
                                   pull(dplyr::select(out5_MVN_2sided, contains("Treatment")))),
                          `Control K2 (Chi2-Distribution)` = c(NA,
                                   pull(dplyr::select(out1_Chi2, contains("Final Control"))),
                                   pull(dplyr::select(out2_Chi2, contains("Control"))),
                                   pull(dplyr::select(out3_Chi2, contains("Control"))),
                                   pull(dplyr::select(out4_Chi2, contains("Control"))),
                                   NA,
                                   pull(dplyr::select(out5_MVN_1sided, contains("Control"))),
                                   pull(dplyr::select(out5_MVN_2sided, contains("Control")))),
                          `Treatment K1 (F-Distribution)` = c(NA,
                                   pull(dplyr::select(out1_F, contains("Final Treatment"))),
                                   pull(dplyr::select(out2_F, contains("Treatment"))),
                                   pull(dplyr::select(out3_F, contains("Treatment"))),
                                   pull(dplyr::select(out4_F, contains("Treatment"))),
                                   NA,
                                   pull(dplyr::select(out5_T_1sided, contains("Treatment"))),
                                   pull(dplyr::select(out5_T_2sided, contains("Treatment")))),
                          `Control K2 (F-Distribution)` = c(NA,
                                   pull(dplyr::select(out1_F, contains("Final Control"))),
                                   pull(dplyr::select(out2_F, contains("Control"))),
                                   pull(dplyr::select(out3_F, contains("Control"))),
                                   pull(dplyr::select(out4_F, contains("Control"))),
                                   NA,
                                   pull(dplyr::select(out5_T_1sided, contains("Control"))),
                                   pull(dplyr::select(out5_T_2sided, contains("Control"))))
                          )

  } else if(output == "m"){
    outputTable <- tibble(`Design Method` = c("1. P-Value Adjustments",
                                              "a. Bonferroni",
                                              "b. Sidak",
                                              "c. D/AP",
                                              "2. Combined Outcomes",
                                              "3. Single Weighted 1-DF Combined Test",
                                              "4. Disjunctive 2-DF Test",
                                              "5. Conjunctive IU Test (MVN and t-dist.)",
                                              "a. 1-Sided Hypothesis",
                                              "b. 2-Sided Hypothesis"),
                          `m (Chi2-distribution)` = c(NA,
                                                      pull(out1_Chi2[, ncol(out1_Chi2)]),
                                                      out2_Chi2,
                                                      out3_Chi2,
                                                      out4_Chi2,
                                                      NA,
                                                      out5_MVN_1sided,
                                                      out5_MVN_2sided),
                          `m (F-distribution)` = c(NA,
                                                   pull(out1_F[, ncol(out1_F)]),
                                                   out2_F,
                                                   out3_F,
                                                   out4_F,
                                                   NA,
                                                   out5_T_1sided,
                                                   out5_T_2sided))
  }

  return(outputTable)
} # End run_crt2_design()
