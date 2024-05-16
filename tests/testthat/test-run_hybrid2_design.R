# Testing script the P-Value Adjustment functions
# `run_hybrid2_design` - Calculates all study design metrics for all methods

# Cases for when we want to calculate K ----------------------------------------
test_that("Inconsistent inputs results in error", {
  expect_error(run_hybrid2_design(output = "K",       # Parameter to calculate
                                  power = 0.9,   # Desired statistical power
                                  K = 10,       # Number of clusters in each arm
                                  m = 10,       # Individuals per cluster
                                  alpha = 0.05, # Significance level
                                  beta1 = 0.2,        # Effect for outcome 1
                                  beta2 = 0.1,        # Effect for outcome 2
                                  varY1 = 0.2,        # Variance for outcome 1
                                  varY2 = 0.1,        # Variance for outcome 2
                                  rho01 = 0.03,        # ICC for outcome 1
                                  rho02 = 0.02,        # ICC for outcome 2
                                  rho1 = 0.01,         # Inter-subject between-endpoint ICC
                                  rho2 = 0.04,         # Intra-subject between-endpoint ICC
                                  r = 1         # Treatment allocation ratio
  ), "'K' variable cannot be defined if desired study design output is 'K'.")
})

test_that("Inconsistent inputs results in error", {
  expect_error(run_hybrid2_design(output = "K",       # Parameter to calculate
                                  power = NA,   # Desired statistical power
                                  K = NA,       # Number of clusters in each arm
                                  m = 10,       # Individuals per cluster
                                  alpha = 0.05, # Significance level
                                  beta1 = 0.2,        # Effect for outcome 1
                                  beta2 = 0.1,        # Effect for outcome 2
                                  varY1 = 0.2,        # Variance for outcome 1
                                  varY2 = 0.1,        # Variance for outcome 2
                                  rho01 = 0.03,        # ICC for outcome 1
                                  rho02 = 0.02,        # ICC for outcome 2
                                  rho1 = 0.01,         # Inter-subject between-endpoint ICC
                                  rho2 = 0.04,         # Intra-subject between-endpoint ICC
                                  r = 1         # Treatment allocation ratio
  ), "Must define 'power' in order to calculate K.")
})

test_that("Inconsistent inputs results in error", {
  expect_error(run_hybrid2_design(output = "K",       # Parameter to calculate
                                  power = 0.9,   # Desired statistical power
                                  K = NA,       # Number of clusters in each arm
                                  m = NA,       # Individuals per cluster
                                  alpha = 0.05, # Significance level
                                  beta1 = 0.2,        # Effect for outcome 1
                                  beta2 = 0.1,        # Effect for outcome 2
                                  varY1 = 0.2,        # Variance for outcome 1
                                  varY2 = 0.1,        # Variance for outcome 2
                                  rho01 = 0.03,        # ICC for outcome 1
                                  rho02 = 0.02,        # ICC for outcome 2
                                  rho1 = 0.01,         # Inter-subject between-endpoint ICC
                                  rho2 = 0.04,         # Intra-subject between-endpoint ICC
                                  r = 1         # Treatment allocation ratio
  ), "Must define 'm' in order to calculate K.")
})

# Cases for when we want to calculate power ------------------------------------
test_that("Inconsistent inputs results in error", {
  expect_error(run_hybrid2_design(output = "power",       # Parameter to calculate
                                  power = 0.9,   # Desired statistical power
                                  K = 10,       # Number of clusters in each arm
                                  m = 10,       # Individuals per cluster
                                  alpha = 0.05, # Significance level
                                  beta1 = 0.2,        # Effect for outcome 1
                                  beta2 = 0.1,        # Effect for outcome 2
                                  varY1 = 0.2,        # Variance for outcome 1
                                  varY2 = 0.1,        # Variance for outcome 2
                                  rho01 = 0.03,        # ICC for outcome 1
                                  rho02 = 0.02,        # ICC for outcome 2
                                  rho1 = 0.01,         # Inter-subject between-endpoint ICC
                                  rho2 = 0.04,         # Intra-subject between-endpoint ICC
                                  r = 1         # Treatment allocation ratio
  ), "'power' variable cannot be defined if desired study design output is 'power'.")
})

test_that("Inconsistent inputs results in error", {
  expect_error(run_hybrid2_design(output = "power",       # Parameter to calculate
                                  power = NA,   # Desired statistical power
                                  K = NA,       # Number of clusters in each arm
                                  m = 10,       # Individuals per cluster
                                  alpha = 0.05, # Significance level
                                  beta1 = 0.2,        # Effect for outcome 1
                                  beta2 = 0.1,        # Effect for outcome 2
                                  varY1 = 0.2,        # Variance for outcome 1
                                  varY2 = 0.1,        # Variance for outcome 2
                                  rho01 = 0.03,        # ICC for outcome 1
                                  rho02 = 0.02,        # ICC for outcome 2
                                  rho1 = 0.01,         # Inter-subject between-endpoint ICC
                                  rho2 = 0.04,         # Intra-subject between-endpoint ICC
                                  r = 1         # Treatment allocation ratio
  ), "Must define 'K' in order to calculate power.")
})

test_that("Inconsistent inputs results in error", {
  expect_error(run_hybrid2_design(output = "power",       # Parameter to calculate
                                  power = NA,   # Desired statistical power
                                  K = 10,       # Number of clusters in each arm
                                  m = NA,       # Individuals per cluster
                                  alpha = 0.05, # Significance level
                                  beta1 = 0.2,        # Effect for outcome 1
                                  beta2 = 0.1,        # Effect for outcome 2
                                  varY1 = 0.2,        # Variance for outcome 1
                                  varY2 = 0.1,        # Variance for outcome 2
                                  rho01 = 0.03,        # ICC for outcome 1
                                  rho02 = 0.02,        # ICC for outcome 2
                                  rho1 = 0.01,         # Inter-subject between-endpoint ICC
                                  rho2 = 0.04,         # Intra-subject between-endpoint ICC
                                  r = 1         # Treatment allocation ratio
  ), "Must define 'm' in order to calculate power.")
})

# Cases for when we want to calcualte m ----------------------------------------

test_that("Inconsistent inputs results in error", {
  expect_error(run_hybrid2_design(output = "m",       # Parameter to calculate
                                  power = 0.9,   # Desired statistical power
                                  K = 10,       # Number of clusters in each arm
                                  m = 10,       # Individuals per cluster
                                  alpha = 0.05, # Significance level
                                  beta1 = 0.2,        # Effect for outcome 1
                                  beta2 = 0.1,        # Effect for outcome 2
                                  varY1 = 0.2,        # Variance for outcome 1
                                  varY2 = 0.1,        # Variance for outcome 2
                                  rho01 = 0.03,        # ICC for outcome 1
                                  rho02 = 0.02,        # ICC for outcome 2
                                  rho1 = 0.01,         # Inter-subject between-endpoint ICC
                                  rho2 = 0.04,         # Intra-subject between-endpoint ICC
                                  r = 1         # Treatment allocation ratio
  ), "'m' variable cannot be defined if desired study design output is 'm'.")
})

test_that("Inconsistent inputs results in error", {
  expect_error(run_hybrid2_design(output = "m",       # Parameter to calculate
                                  power = NA,   # Desired statistical power
                                  K = 10,       # Number of clusters in each arm
                                  m = NA,       # Individuals per cluster
                                  alpha = 0.05, # Significance level
                                  beta1 = 0.2,        # Effect for outcome 1
                                  beta2 = 0.1,        # Effect for outcome 2
                                  varY1 = 0.2,        # Variance for outcome 1
                                  varY2 = 0.1,        # Variance for outcome 2
                                  rho01 = 0.03,        # ICC for outcome 1
                                  rho02 = 0.02,        # ICC for outcome 2
                                  rho1 = 0.01,         # Inter-subject between-endpoint ICC
                                  rho2 = 0.04,         # Intra-subject between-endpoint ICC
                                  r = 1         # Treatment allocation ratio
  ), "Must define 'power' in order to calculate m.")
})

test_that("Inconsistent inputs results in error", {
  expect_error(run_hybrid2_design(output = "m",       # Parameter to calculate
                                  power = 0.9,   # Desired statistical power
                                  K = NA,       # Number of clusters in each arm
                                  m = NA,       # Individuals per cluster
                                  alpha = 0.05, # Significance level
                                  beta1 = 0.2,        # Effect for outcome 1
                                  beta2 = 0.1,        # Effect for outcome 2
                                  varY1 = 0.2,        # Variance for outcome 1
                                  varY2 = 0.1,        # Variance for outcome 2
                                  rho01 = 0.03,        # ICC for outcome 1
                                  rho02 = 0.02,        # ICC for outcome 2
                                  rho1 = 0.01,         # Inter-subject between-endpoint ICC
                                  rho2 = 0.04,         # Intra-subject between-endpoint ICC
                                  r = 1         # Treatment allocation ratio
  ), "Must define 'K' in order to calculate m.")
})

