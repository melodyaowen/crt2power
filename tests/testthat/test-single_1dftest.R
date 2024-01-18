# Testing script the Single 1-DF Combined Test functions
# `calc_pwr_single_1dftest()` - calculates power
# `calc_K_single_1dftest()`   - calculates # of clusters per treatment group
# `calc_m_single_1dftest()`   - calculates cluster size

# `calc_pwr_single_1dftest()` --------------------------------------------------

# Test for single 1 df test power
test_that("Single 1-DF test power works", {
  expect_equal(calc_pwr_single_1dftest(K = 15,
                                       m = 300,
                                       alpha = 0.05,
                                       beta1 = 0.1,
                                       beta2 = 0.1,
                                       varY1 = 0.23,
                                       varY2 = 0.25,
                                       rho01 = 0.025,
                                       rho02 = 0.025,
                                       rho1 = 0.01,
                                       rho2 = 0.05,
                                       r = 1),
               0.9811)
})

# `calc_K_single_1dftest()` ----------------------------------------------------

# Test for single 1 df test treatment group K
test_that("Single 1-DF test treatment group 'K' works", {
  expect_equal(calc_K_single_1dftest(power = 0.8,
                                     m = 300,
                                     alpha = 0.05,
                                     beta1 = 0.1,
                                     beta2 = 0.1,
                                     varY1 = 0.23,
                                     varY2 = 0.25,
                                     rho01 = 0.025,
                                     rho02 = 0.025,
                                     rho1 = 0.01,
                                     rho2 = 0.05,
                                     r = 1)$`Treatment (K)`,
               8)
})

# Test for single 1 df test control group K
test_that("Single 1-DF test control group 'K' works", {
  expect_equal(calc_K_single_1dftest(power = 0.8,
                                     m = 300,
                                     alpha = 0.05,
                                     beta1 = 0.1,
                                     beta2 = 0.1,
                                     varY1 = 0.23,
                                     varY2 = 0.25,
                                     rho01 = 0.025,
                                     rho02 = 0.025,
                                     rho1 = 0.01,
                                     rho2 = 0.05,
                                     r = 1)$`Control (K)`,
               8)
})

# `calc_m_single_1dftest()` ----------------------------------------------------

# Test for single 1 df test m
test_that("Single 1-DF test 'm' works", {
  expect_equal(calc_m_single_1dftest(power = 0.8,
                                     K = 15,
                                     alpha = 0.05,
                                     beta1 = 0.1,
                                     beta2 = 0.1,
                                     varY1 = 0.23,
                                     varY2 = 0.25,
                                     rho01 = 0.025,
                                     rho02 = 0.025,
                                     rho1 = 0.01,
                                     rho2 = 0.05,
                                     r = 1),
               23)
})


# Testing that all 3 functions align -------------------------------------------

# K aligns with power
test_that("Single 1-DF test K calculation aligns with power", {
  result <- calc_K_single_1dftest(power = calc_pwr_single_1dftest(K = 15,
                                                                  m = 300,
                                                                  alpha = 0.05,
                                                                  beta1 = 0.1,
                                                                  beta2 = 0.1,
                                                                  varY1 = 0.23,
                                                                  varY2 = 0.25,
                                                                  rho01 = 0.025,
                                                                  rho02 = 0.025,
                                                                  rho1 = 0.01,
                                                                  rho2 = 0.05,
                                                                  r = 1),
                                  m = 300,
                                  alpha = 0.05,
                                  beta1 = 0.1,
                                  beta2 = 0.1,
                                  varY1 = 0.23,
                                  varY2 = 0.25,
                                  rho01 = 0.025,
                                  rho02 = 0.025,
                                  rho1 = 0.01,
                                  rho2 = 0.05,
                                  r = 1)$`Treatment (K)`
  expect_true(result >= 15 & result <= 16,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# K aligns with m
test_that("Single 1-DF test K calculation aligns with power", {
  result <- calc_K_single_1dftest(power = 0.8,
                                  m = calc_m_single_1dftest(power = 0.8,
                                                            K = 15,
                                                            alpha = 0.05,
                                                            beta1 = 0.1,
                                                            beta2 = 0.1,
                                                            varY1 = 0.23,
                                                            varY2 = 0.25,
                                                            rho01 = 0.025,
                                                            rho02 = 0.025,
                                                            rho1 = 0.01,
                                                            rho2 = 0.05,
                                                            r = 1),
                                  alpha = 0.05,
                                  beta1 = 0.1,
                                  beta2 = 0.1,
                                  varY1 = 0.23,
                                  varY2 = 0.25,
                                  rho01 = 0.025,
                                  rho02 = 0.025,
                                  rho1 = 0.01,
                                  rho2 = 0.05,
                                  r = 1)$`Treatment (K)`
  expect_true(result >= 15 & result <= 16,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# Power aligns with K
test_that("Single 1-DF test K calculation aligns with power", {
  result <- calc_m_single_1dftest(power = calc_pwr_single_1dftest(m = 300,
                                                                  K = 15,
                                                                  alpha = 0.05,
                                                                  beta1 = 0.1,
                                                                  beta2 = 0.1,
                                                                  varY1 = 0.23,
                                                                  varY2 = 0.25,
                                                                  rho01 = 0.025,
                                                                  rho02 = 0.025,
                                                                  rho1 = 0.01,
                                                                  rho2 = 0.05,
                                                                  r = 1),
                                  K = 15,
                                  alpha = 0.05,
                                  beta1 = 0.1,
                                  beta2 = 0.1,
                                  varY1 = 0.23,
                                  varY2 = 0.25,
                                  rho01 = 0.025,
                                  rho02 = 0.025,
                                  rho1 = 0.01,
                                  rho2 = 0.05,
                                  r = 1)
  expect_true(result >= 299 & result <= 301,
              info = "The result should be 300, but with rounding between 299 and 301")
})
