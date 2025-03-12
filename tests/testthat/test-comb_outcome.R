# Testing script the Combined Outcomes Approach functions
# `calc_pwr_comb_outcome()`  - calculates power
# `calc_K_comb_outcome()`    - calculates # of clusters per treatment group
# `calc_m_comb_outcome()`    - calculates cluster size


# `calc_pwr_comb_outcome()` ---------------------------------------------------

# Test for combined outcomes power
test_that("Combined outcomes power works (Chi2 dist)", {
  expect_equal(calc_pwr_comb_outcome(dist = "Chi2",
                                     K = 15,
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

test_that("Combined outcomes power works (F dist)", {
  expect_equal(calc_pwr_comb_outcome(dist = "F",
                                     K = 15,
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
               0.9727)
})

# `calc_K_comb_outcome()` ------------------------------------------------------

# Test for combined outcomes treatment group K
test_that("Combined outcomes treatment group 'K' works (Chi2 dist)", {
  expect_equal(calc_K_comb_outcome(dist = "Chi2",
                                   power = 0.8,
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

test_that("Combined outcomes treatment group 'K' works (F dist)", {
  expect_equal(calc_K_comb_outcome(dist = "F",
                                   power = 0.8,
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
               9)
})

# Test for combined outcomes control group K
test_that("Combined outcomes control group 'K' works (Chi2 dist)", {
  expect_equal(calc_K_comb_outcome(dist = "Chi2",
                                   power = 0.8,
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

test_that("Combined outcomes control group 'K' works (F dist)", {
  expect_equal(calc_K_comb_outcome(dist = "F",
                                   power = 0.8,
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
               9)
})

# `calc_m_comb_outcome()` ------------------------------------------------------

# Test for combined outcomes m
test_that("Combined outcomes 'm' works (Chi2 dist)", {
  expect_equal(calc_m_comb_outcome(dist = "Chi2",
                                   power = 0.8,
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

test_that("Combined outcomes 'm' works (F dist)", {
  expect_equal(calc_m_comb_outcome(dist = "F",
                                   power = 0.8,
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
               27)
})

# Throws an error when 'm' cannot be large enough to attain desired power
test_that("Combined outcome final 'm' throws error when desired power cannot be attained (F dist)", {
  expect_error(calc_m_comb_outcome(dist = "F",
                                   power = 0.9,
                                   K = 4,
                                   alpha = 0.05,
                                   beta1 = 0.1,
                                   beta2 = 0.1,
                                   varY1 = 0.23,
                                   varY2 = 0.25,
                                   rho01 = 0.025,
                                   rho02 = 0.025,
                                   rho1 = 0.01,
                                   rho2 = 0.05,
                                   r = 1))
})

test_that("Combined outcome final 'm' throws error when desired power cannot be attained (Chi2 dist)", {
  expect_error(calc_m_comb_outcome(dist = "Chi2",
                                   power = 0.9,
                                   K = 4,
                                   alpha = 0.05,
                                   beta1 = 0.1,
                                   beta2 = 0.1,
                                   varY1 = 0.23,
                                   varY2 = 0.25,
                                   rho01 = 0.025,
                                   rho02 = 0.025,
                                   rho1 = 0.01,
                                   rho2 = 0.05,
                                   r = 1))
})

# Testing that all 3 functions align -------------------------------------------

# K aligns with power for Chi2 distribution
test_that("Combined outcomes K calculation aligns with power (Chi2 dist)", {
  result <- calc_K_comb_outcome(dist = "Chi2",
                                power = calc_pwr_comb_outcome(dist = "Chi2",
                                                              K = 15,
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

# K aligns with m for Chi2 distribution
test_that("Combined outcomes K calculation aligns with power (Chi2 dist)", {
  result <- calc_K_comb_outcome(dist = "Chi2",
                                power = 0.8,
                                m = calc_m_comb_outcome(dist = "Chi2",
                                                        power = 0.8,
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

# Power aligns with K for Chi2 distribution
test_that("Combined outcomes K calculation aligns with power (Chi2 dist)", {
  result <- calc_m_comb_outcome(dist = "Chi2",
                                power = calc_pwr_comb_outcome(dist = "Chi2",
                                                              m = 300,
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




# K aligns with power for F distribution
test_that("Combined outcomes K calculation aligns with power (F dist)", {
  result <- calc_K_comb_outcome(dist = "F",
                                power = calc_pwr_comb_outcome(dist = "F",
                                                              K = 15,
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

# K aligns with m for F distribution
test_that("Combined outcomes K calculation aligns with power (F dist)", {
  result <- calc_K_comb_outcome(dist = "F",
                                power = 0.8,
                                m = calc_m_comb_outcome(dist = "F",
                                                        power = 0.8,
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

# Power aligns with K for F distribution
test_that("Combined outcomes K calculation aligns with power (F dist)", {
  result <- calc_m_comb_outcome(dist = "F",
                                power = calc_pwr_comb_outcome(dist = "F",
                                                              m = 300,
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
  expect_true(result >= 298 & result <= 301,
              info = "The result should be 300, but with rounding between 298 and 301")
})



