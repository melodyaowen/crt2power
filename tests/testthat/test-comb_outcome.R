# Testing script the Combined Outcomes Approach functions
# `calc_pwr_comb_outcome()`  - calculates power
# `calc_K_comb_outcome()`    - calculates # of clusters per treatment group
# `calc_m_comb_outcome()`    - calculates cluster size


# `calc_pwr_comb_outcome()` ---------------------------------------------------

# Test for combined outcomes power
test_that("Combined outcomes power works", {
  expect_equal(calc_pwr_comb_outcome(K = 15,
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
               0.9810)
})

# `calc_K_comb_outcome()` ------------------------------------------------------

# Test for combined outcomes treatment group K
test_that("Combined outcomes treatment group 'K' works", {
  expect_equal(calc_K_comb_outcome(power = 0.8,
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

# Test for combined outcomes control group K
test_that("Combined outcomes control group 'K' works", {
  expect_equal(calc_K_comb_outcome(power = 0.8,
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

# `calc_m_comb_outcome()` ------------------------------------------------------

# Test for combined outcomes m
test_that("Combined outcomes 'm' works", {
  expect_equal(calc_m_comb_outcome(power = 0.8,
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
test_that("Combined outcomes K calculation aligns with power", {
  result <- calc_K_comb_outcome(power = calc_pwr_comb_outcome(K = 15,
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
test_that("Combined outcomes K calculation aligns with power", {
  result <- calc_K_comb_outcome(power = 0.8,
                                m = calc_m_comb_outcome(power = 0.8,
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
test_that("Combined outcomes K calculation aligns with power", {
  result <- calc_m_comb_outcome(power = calc_pwr_comb_outcome(m = 300,
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



