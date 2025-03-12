# Testing script the Disjunctive 2-DF Test functions
# `calc_pwr_disj_2dftest()`  - calculates power
# `calc_K_disj_2dftest()`    - calculates # of clusters per treatment group
# `calc_m_disj_2dftest()`    - calculates cluster size

# `calc_pwr_disj_2dftest()` ----------------------------------------------------

# Test for Disjunctive 2-DF test power for Chi2 distribution
test_that("Disjunctive 2-DF test power works (Chi2 dist)", {
  expect_equal(calc_pwr_disj_2dftest(dist = "Chi2",
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
               0.9601)
})

# Test for Disjunctive 2-DF test power for F distribution
test_that("Disjunctive 2-DF test power works (F dist)", {
  expect_equal(calc_pwr_disj_2dftest(dist = "F",
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
               0.9363)
})

# `calc_K_disj_2dftest()` ------------------------------------------------------

# Test for Disjunctive 2-DF test treatment group K for Chi2 distribution
test_that("Disjunctive 2-DF test treatment group 'K' works (Chi2 dist)", {
  expect_equal(calc_K_disj_2dftest(dist = "Chi2",
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

# Test for Disjunctive 2-DF test control group K for Chi2 distribution
test_that("Disjunctive 2-DF test control group 'K' works (Chi2 dist)", {
  expect_equal(calc_K_disj_2dftest(dist = "Chi2",
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

# Test for Disjunctive 2-DF test treatment group K for F distribution
test_that("Disjunctive 2-DF test treatment group 'K' works (F dist)", {
  expect_equal(calc_K_disj_2dftest(dist = "F",
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
               11)
})

# Test for Disjunctive 2-DF test control group K for F distribution
test_that("Disjunctive 2-DF test control group 'K' works (F dist)", {
  expect_equal(calc_K_disj_2dftest(dist = "F",
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
               11)
})

# `calc_m_disj_2dftest()` ------------------------------------------------------

# Test for Disjunctive 2-DF test m for Chi2 distribution
test_that("Disjunctive 2-DF test 'm' works (Chi2 dist)", {
  expect_equal(calc_m_disj_2dftest(dist = "Chi2",
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
               34)
})

# Test for Disjunctive 2-DF test m for F distribution
test_that("Disjunctive 2-DF test 'm' works (F dist)", {
  expect_equal(calc_m_disj_2dftest(dist = "F",
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
               45)
})


# Testing that all 3 functions align -------------------------------------------

# K aligns with power for Chi2 distribution
test_that("Disjunctive 2-DF test K calculation aligns with power (Chi2 dist)", {
  result <- calc_K_disj_2dftest(dist = "Chi2",
                                power = calc_pwr_disj_2dftest(dist = "Chi2",
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
test_that("Disjunctive 2-DF test K calculation aligns with power (Chi2 dist)", {
  result <- calc_K_disj_2dftest(dist = "Chi2",
                                power = 0.8,
                                m = calc_m_disj_2dftest(dist = "Chi2",
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
test_that("Disjunctive 2-DF test K calculation aligns with power (Chi2 dist)", {
  result <- calc_m_disj_2dftest(dist = "Chi2",
                                power = calc_pwr_disj_2dftest(dist = "Chi2",
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
test_that("Disjunctive 2-DF test K calculation aligns with power (F dist)", {
  result <- calc_K_disj_2dftest(dist = "F",
                                power = calc_pwr_disj_2dftest(dist = "F",
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
test_that("Disjunctive 2-DF test K calculation aligns with power (F dist)", {
  result <- calc_K_disj_2dftest(dist = "F",
                                power = 0.8,
                                m = calc_m_disj_2dftest(dist = "F",
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
test_that("Disjunctive 2-DF test K calculation aligns with power (F dist)", {
  result <- calc_m_disj_2dftest(dist = "F",
                                power = calc_pwr_disj_2dftest(dist = "F",
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



