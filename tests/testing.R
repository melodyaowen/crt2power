# Script for testing
rm(list = ls(all.names = TRUE))
setwd("/Users/melodyowen/Desktop/GitHub/hybrid2power/tests")
source("../R/calc_ncp_chi2.R")
source("../R/calc_pwr_comb_outcome.R")
source("../R/calc_pwr_conj_test.R")
source("../R/calc_pwr_disj_2dftest.R")
source("../R/calc_pwr_pval_adj.R")
source("../R/calc_pwr_single_1dftest.R")

# Method 1: P-Value Adjustments
calc_pwr_pval_adj(K     = 15,
                  m     = 300,
                  alpha = 0.05,
                  beta1 = 0.1,
                  beta2 = 0.1,
                  varY1 = 0.23,
                  varY2 = 0.25,
                  rho01 = 0.025,
                  rho02 = 0.025,
                  rho2  = 0.05)
calc_K_pval_adj(power = 0.8,
                m     = 300,
                alpha = 0.05,
                beta1 = 0.1,
                beta2 = 0.1,
                varY1 = 0.23,
                varY2 = 0.25,
                rho01 = 0.025,
                rho02 = 0.025,
                rho2  = 0.05)
calc_m_pval_adj(power = 0.8,
                K     = 15,
                alpha = 0.05,
                beta1 = 0.1,
                beta2 = 0.1,
                varY1 = 0.23,
                varY2 = 0.25,
                rho01 = 0.025,
                rho02 = 0.025,
                rho2  = 0.05)


# Method 2: Combined Outcome
calc_pwr_comb_outcome(K     = 15,
                      m     = 300,
                      alpha = 0.05,
                      beta1 = 0.1,
                      beta2 = 0.1,
                      varY1 = 0.23,
                      varY2 = 0.25,
                      rho01 = 0.025,
                      rho02 = 0.025,
                      rho1  = 0.01,
                      rho2  = 0.05)
calc_K_comb_outcome(power = 0.8,
                    m     = 300,
                    alpha = 0.05,
                    beta1 = 0.1,
                    beta2 = 0.1,
                    varY1 = 0.23,
                    varY2 = 0.25,
                    rho01 = 0.025,
                    rho02 = 0.025,
                    rho1  = 0.01,
                    rho2  = 0.05)
calc_m_comb_outcome(power = 0.8,
                    K     = 15,
                    alpha = 0.05,
                    beta1 = 0.1,
                    beta2 = 0.1,
                    varY1 = 0.23,
                    varY2 = 0.25,
                    rho01 = 0.025,
                    rho02 = 0.025,
                    rho1  = 0.01,
                    rho2  = 0.05)

# Method 3: Single 1-DF Test
calc_pwr_single_1dftest(K     = 15,
                        m     = 300,
                        alpha = 0.05,
                        beta1 = 0.1,
                        beta2 = 0.1,
                        varY1 = 0.23,
                        varY2 = 0.25,
                        rho01 = 0.025,
                        rho02 = 0.025,
                        rho1  = 0.01,
                        rho2  = 0.05)
calc_K_single_1dftest(power = 0.8,
                      m     = 300,
                      alpha = 0.05,
                      beta1 = 0.1,
                      beta2 = 0.1,
                      varY1 = 0.23,
                      varY2 = 0.25,
                      rho01 = 0.025,
                      rho02 = 0.025,
                      rho1  = 0.01,
                      rho2  = 0.05)
calc_m_single_1dftest(power = 0.8,
                      K     = 15,
                      alpha = 0.05,
                      beta1 = 0.1,
                      beta2 = 0.1,
                      varY1 = 0.23,
                      varY2 = 0.25,
                      rho01 = 0.025,
                      rho02 = 0.025,
                      rho1  = 0.01,
                      rho2  = 0.05)

# Method 4: Disjunctive 2-DF Test
calc_pwr_disj_2dftest(K     = 15,
                      m     = 300,
                      alpha = 0.05,
                      beta1 = 0.1,
                      beta2 = 0.1,
                      varY1 = 0.23,
                      varY2 = 0.25,
                      rho01 = 0.025,
                      rho02 = 0.025,
                      rho1  = 0.01,
                      rho2  = 0.05)
calc_K_disj_2dftest(power = 0.8,
                    m     = 300,
                    alpha = 0.05,
                    beta1 = 0.1,
                    beta2 = 0.1,
                    varY1 = 0.23,
                    varY2 = 0.25,
                    rho01 = 0.025,
                    rho02 = 0.025,
                    rho1  = 0.01,
                    rho2  = 0.05)
calc_m_disj_2dftest(power = 0.8,
                    K     = 15,
                    alpha = 0.05,
                    beta1 = 0.1,
                    beta2 = 0.1,
                    varY1 = 0.23,
                    varY2 = 0.25,
                    rho01 = 0.025,
                    rho02 = 0.025,
                    rho1  = 0.01,
                    rho2  = 0.05)

# Method 5: Conjunctive IU Test
calc_pwr_conj_test(K     = 15,
                   m     = 300,
                   alpha = 0.05,
                   beta1 = 0.1,
                   beta2 = 0.1,
                   varY1 = 0.23,
                   varY2 = 0.25,
                   rho01 = 0.025,
                   rho02 = 0.025,
                   rho1  = 0.01,
                   rho2  = 0.05)
calc_K_conj_test(power = 0.8,
                 m     = 300,
                 alpha = 0.05,
                 beta1 = 0.1,
                 beta2 = 0.1,
                 varY1 = 0.23,
                 varY2 = 0.25,
                 rho01 = 0.025,
                 rho02 = 0.025,
                 rho1  = 0.01,
                 rho2  = 0.05)
calc_m_conj_test(power = 0.8,
                 K     = 15,
                 alpha = 0.05,
                 beta1 = 0.1,
                 beta2 = 0.1,
                 varY1 = 0.23,
                 varY2 = 0.25,
                 rho01 = 0.025,
                 rho02 = 0.025,
                 rho1  = 0.01,
                 rho2  = 0.05)



rm(list = ls(all.names = TRUE))
setwd("/Users/melodyowen/Desktop/GitHub/hybrid2power/tests")
source("../R/calc_ncp_chi2.R")
source("../R/calc_pwr_comb_outcome.R")
source("../R/calc_pwr_conj_test.R")
source("../R/calc_pwr_disj_2dftest.R")
source("../R/calc_pwr_pval_adj.R")
source("../R/calc_pwr_single_1dftest.R")
run_hybrid2_design(output = "m",
                   power = 0.8,
                   K     = 15,
                   #m     = 300,
                   alpha = 0.05,
                   beta1 = 0.1,
                   beta2 = 0.1,
                   varY1 = 0.23,
                   varY2 = 0.25,
                   rho01 = 0.025,
                   rho02 = 0.025,
                   rho1  = 0.01,
                   rho2  = 0.05)


power = 0.8; K = 15; m = 300; alpha = 0.05
beta1 = 0.1; beta2 = 0.1; varY1 = 0.23; varY2 = 0.25
rho01 = 0.025; rho02 = 0.025; rho1  = 0.01; rho2  = 0.05












power = 0.8
K     = 15
alpha = 0.05
beta1 = 0.1
beta2 = 0.1
varY1 = 0.23
varY2 = 0.25
rho01 = 0.025
rho02 = 0.025
rho1  = 0.01
rho2  = 0.05



