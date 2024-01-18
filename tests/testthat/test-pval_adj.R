# Testing script the P-Value Adjustment functions
# `calc_pwr_pval_adj()`  - calculates power
# `calc_K_pval_adj()`    - calculates # of clusters per treatment group
# `calc_m_pval_adj()`    - calculates cluster size


# Tests for `calc_pwer_pval_adj()` ---------------------------------------------

# First outcome power
test_that("P-value power for Y1 works", {
  expect_equal(calc_pwr_pval_adj(K = 15,
                                 m = 300,
                                 alpha = 0.05,
                                 beta1 = 0.1,
                                 beta2 = 0.1,
                                 varY1 = 0.23,
                                 varY2 = 0.25,
                                 rho01 = 0.025,
                                 rho02 = 0.025,
                                 rho2 = 0.05,
                                 r = 1)$`Power (Y1)`,
               c(0.8762, 0.8772, 0.8799))
})

# Second outcome power
test_that("P-value power for Y2 works", {
  expect_equal(calc_pwr_pval_adj(K = 15,
                                 m = 300,
                                 alpha = 0.05,
                                 beta1 = 0.1,
                                 beta2 = 0.1,
                                 varY1 = 0.23,
                                 varY2 = 0.25,
                                 rho01 = 0.025,
                                 rho02 = 0.025,
                                 rho2 = 0.05,
                                 r = 1)$`Power (Y2)`,
               c(0.8455, 0.8467, 0.8498))
})

# Final power
test_that("P-value final power works", {
  expect_equal(calc_pwr_pval_adj(K = 15,
                                 m = 300,
                                 alpha = 0.05,
                                 beta1 = 0.1,
                                 beta2 = 0.1,
                                 varY1 = 0.23,
                                 varY2 = 0.25,
                                 rho01 = 0.025,
                                 rho02 = 0.025,
                                 rho2 = 0.05,
                                 r = 1)$`Final Power`,
               c(0.8455, 0.8467, 0.8498))
})

# Tests for `calc_K_pval_adj()` ------------------------------------------------

# Treatment group K for Y1
test_that("P-value treatment group 'K' for Y1 works", {
  expect_equal(calc_K_pval_adj(power = 0.8,
                               m = 300,
                               alpha = 0.05,
                               beta1 = 0.1,
                               beta2 = 0.1,
                               varY1 = 0.23,
                               varY2 = 0.25,
                               rho01 = 0.025,
                               rho02 = 0.025,
                               rho2 = 0.05,
                               r = 1)$`Treatment (K) for Y1`,
               c(13, 13, 13))
})

# Treatment group K for Y2
test_that("P-value treatment group 'K' for Y2 works", {
  expect_equal(calc_K_pval_adj(power = 0.8,
                               m = 300,
                               alpha = 0.05,
                               beta1 = 0.1,
                               beta2 = 0.1,
                               varY1 = 0.23,
                               varY2 = 0.25,
                               rho01 = 0.025,
                               rho02 = 0.025,
                               rho2 = 0.05,
                               r = 1)$`Treatment (K) for Y2`,
               c(14, 14, 14))
})

# Control group K for Y1
test_that("P-value control group 'K' for Y1 works", {
  expect_equal(calc_K_pval_adj(power = 0.8,
                               m = 300,
                               alpha = 0.05,
                               beta1 = 0.1,
                               beta2 = 0.1,
                               varY1 = 0.23,
                               varY2 = 0.25,
                               rho01 = 0.025,
                               rho02 = 0.025,
                               rho2 = 0.05,
                               r = 1)$`Control (K) for Y1`,
               c(13, 13, 13))
})

# Control group K for Y2
test_that("P-value control group 'K' for Y2 works", {
  expect_equal(calc_K_pval_adj(power = 0.8,
                               m = 300,
                               alpha = 0.05,
                               beta1 = 0.1,
                               beta2 = 0.1,
                               varY1 = 0.23,
                               varY2 = 0.25,
                               rho01 = 0.025,
                               rho02 = 0.025,
                               rho2 = 0.05,
                               r = 1)$`Control (K) for Y2`,
               c(14, 14, 14))
})

# Final treatment group K
test_that("P-value final treatment group 'K' works", {
  expect_equal(calc_K_pval_adj(power = 0.8,
                               m = 300,
                               alpha = 0.05,
                               beta1 = 0.1,
                               beta2 = 0.1,
                               varY1 = 0.23,
                               varY2 = 0.25,
                               rho01 = 0.025,
                               rho02 = 0.025,
                               rho2 = 0.05,
                               r = 1)$`Final Treatment (K)`,
               c(14, 14, 14))
})

# Final control group K
test_that("P-value final control group 'K' works", {
  expect_equal(calc_K_pval_adj(power = 0.8,
                               m = 300,
                               alpha = 0.05,
                               beta1 = 0.1,
                               beta2 = 0.1,
                               varY1 = 0.23,
                               varY2 = 0.25,
                               rho01 = 0.025,
                               rho02 = 0.025,
                               rho2 = 0.05,
                               r = 1)$`Final Control (K)`,
               c(14, 14, 14))
})

# Tests for `calc_m_pval_adj()` ------------------------------------------------

# First outcome m
test_that("P-value 'm' for Y1 works", {
  expect_equal(calc_m_pval_adj(power = 0.8,
                               K = 15,
                               alpha = 0.05,
                               beta1 = 0.1,
                               beta2 = 0.1,
                               varY1 = 0.23,
                               varY2 = 0.25,
                               rho01 = 0.025,
                               rho02 = 0.025,
                               rho2 = 0.05,
                               r = 1)$`m (Y1)`,
               c(105, 104, 101))
})

# Second outcome m
test_that("P-value 'm' for Y2 works", {
  expect_equal(calc_m_pval_adj(power = 0.8,
                               K = 15,
                               alpha = 0.05,
                               beta1 = 0.1,
                               beta2 = 0.1,
                               varY1 = 0.23,
                               varY2 = 0.25,
                               rho01 = 0.025,
                               rho02 = 0.025,
                               rho2 = 0.05,
                               r = 1)$`m (Y2)`,
               c(149, 147, 141))
})

# Final m
test_that("P-value final 'm' works", {
  expect_equal(calc_m_pval_adj(power = 0.8,
                               K = 15,
                               alpha = 0.05,
                               beta1 = 0.1,
                               beta2 = 0.1,
                               varY1 = 0.23,
                               varY2 = 0.25,
                               rho01 = 0.025,
                               rho02 = 0.025,
                               rho2 = 0.05,
                               r = 1)$`Final m`,
               c(149, 147, 141))
})

# Testing that all 3 functions align -------------------------------------------

# Power and `m` align for Bonferroni
test_that("P-value power aligns with m for Bonferroni", {
  result <- calc_m_pval_adj(power = calc_pwr_pval_adj(K = 15,
                                                      m = 300,
                                                      alpha = 0.05,
                                                      beta1 = 0.1,
                                                      beta2 = 0.1,
                                                      varY1 = 0.23,
                                                      varY2 = 0.25,
                                                      rho01 = 0.025,
                                                      rho02 = 0.025,
                                                      rho2 = 0.05,
                                                      r = 1)$`Final Power`[1],
                            K = 15,
                            alpha = 0.05,
                            beta1 = 0.1,
                            beta2 = 0.1,
                            varY1 = 0.23,
                            varY2 = 0.25,
                            rho01 = 0.025,
                            rho02 = 0.025,
                            rho2 = 0.05,
                            r = 1)$`Final m`[1]
  expect_true(result >= 300 & result <= 301,
               info = "The result should be 300, but with rounding between 300 and 301")
})

# Power and `m` align for Sidak
test_that("P-value power aligns with m for Sidak", {
  result <- calc_m_pval_adj(power = calc_pwr_pval_adj(K = 15,
                                                      m = 300,
                                                      alpha = 0.05,
                                                      beta1 = 0.1,
                                                      beta2 = 0.1,
                                                      varY1 = 0.23,
                                                      varY2 = 0.25,
                                                      rho01 = 0.025,
                                                      rho02 = 0.025,
                                                      rho2 = 0.05,
                                                      r = 1)$`Final Power`[2],
                            K = 15,
                            alpha = 0.05,
                            beta1 = 0.1,
                            beta2 = 0.1,
                            varY1 = 0.23,
                            varY2 = 0.25,
                            rho01 = 0.025,
                            rho02 = 0.025,
                            rho2 = 0.05,
                            r = 1)$`Final m`[2]
  expect_true(result >= 300 & result <= 301,
              info = "The result should be 300, but with rounding between 300 and 301")
})

# Power and `m` align for D/AP
test_that("P-value power aligns with m for D/AP", {
  result <- calc_m_pval_adj(power = calc_pwr_pval_adj(K = 15,
                                                      m = 300,
                                                      alpha = 0.05,
                                                      beta1 = 0.1,
                                                      beta2 = 0.1,
                                                      varY1 = 0.23,
                                                      varY2 = 0.25,
                                                      rho01 = 0.025,
                                                      rho02 = 0.025,
                                                      rho2 = 0.05,
                                                      r = 1)$`Final Power`[3],
                            K = 15,
                            alpha = 0.05,
                            beta1 = 0.1,
                            beta2 = 0.1,
                            varY1 = 0.23,
                            varY2 = 0.25,
                            rho01 = 0.025,
                            rho02 = 0.025,
                            rho2 = 0.05,
                            r = 1)$`Final m`[3]
  expect_true(result >= 300 & result <= 301,
              info = "The result should be 300, but with rounding between 300 and 301")
})

# Power and `K` align for Bonferroni
test_that("P-value power aligns with K for Bonferroni", {
  result <- calc_K_pval_adj(power = calc_pwr_pval_adj(K = 15,
                                                      m = 300,
                                                      alpha = 0.05,
                                                      beta1 = 0.1,
                                                      beta2 = 0.1,
                                                      varY1 = 0.23,
                                                      varY2 = 0.25,
                                                      rho01 = 0.025,
                                                      rho02 = 0.025,
                                                      rho2 = 0.05,
                                                      r = 1)$`Final Power`[1],
                            m = 300,
                            alpha = 0.05,
                            beta1 = 0.1,
                            beta2 = 0.1,
                            varY1 = 0.23,
                            varY2 = 0.25,
                            rho01 = 0.025,
                            rho02 = 0.025,
                            rho2 = 0.05,
                            r = 1)$`Final Treatment (K)`[1]
  expect_true(result >= 15 & result <= 16,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# Power and `K` align for Sidak
test_that("P-value power aligns with K for Sidak", {
  result <- calc_K_pval_adj(power = calc_pwr_pval_adj(K = 15,
                                                      m = 300,
                                                      alpha = 0.05,
                                                      beta1 = 0.1,
                                                      beta2 = 0.1,
                                                      varY1 = 0.23,
                                                      varY2 = 0.25,
                                                      rho01 = 0.025,
                                                      rho02 = 0.025,
                                                      rho2 = 0.05,
                                                      r = 1)$`Final Power`[2],
                            m = 300,
                            alpha = 0.05,
                            beta1 = 0.1,
                            beta2 = 0.1,
                            varY1 = 0.23,
                            varY2 = 0.25,
                            rho01 = 0.025,
                            rho02 = 0.025,
                            rho2 = 0.05,
                            r = 1)$`Final Treatment (K)`[2]
  expect_true(result >= 15 & result <= 16,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# Power and `K` align for D/AP
test_that("P-value power aligns with K for D/AP", {
  result <- calc_K_pval_adj(power = calc_pwr_pval_adj(K = 15,
                                                      m = 300,
                                                      alpha = 0.05,
                                                      beta1 = 0.1,
                                                      beta2 = 0.1,
                                                      varY1 = 0.23,
                                                      varY2 = 0.25,
                                                      rho01 = 0.025,
                                                      rho02 = 0.025,
                                                      rho2 = 0.05,
                                                      r = 1)$`Final Power`[3],
                            m = 300,
                            alpha = 0.05,
                            beta1 = 0.1,
                            beta2 = 0.1,
                            varY1 = 0.23,
                            varY2 = 0.25,
                            rho01 = 0.025,
                            rho02 = 0.025,
                            rho2 = 0.05,
                            r = 1)$`Final Treatment (K)`[3]
  expect_true(result >= 15 & result <= 16,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# `m` and `K` align for Bonferroni
test_that("P-value m aligns with K for Bonferroni", {
  result <- calc_K_pval_adj(power = 0.8,
                            m = calc_m_pval_adj(K = 15,
                                                power = 0.8,
                                                alpha = 0.05,
                                                beta1 = 0.1,
                                                beta2 = 0.1,
                                                varY1 = 0.23,
                                                varY2 = 0.25,
                                                rho01 = 0.025,
                                                rho02 = 0.025,
                                                rho2 = 0.05,
                                                r = 1)$`Final m`[1],
                            alpha = 0.05,
                            beta1 = 0.1,
                            beta2 = 0.1,
                            varY1 = 0.23,
                            varY2 = 0.25,
                            rho01 = 0.025,
                            rho02 = 0.025,
                            rho2 = 0.05,
                            r = 1)$`Final Treatment (K)`[1]
  expect_true(result >= 15 & result <= 16,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# `m` and `K` align for Sidak
test_that("P-value m aligns with K for Sidak", {
  result <- calc_K_pval_adj(power = 0.8,
                            m = calc_m_pval_adj(K = 15,
                                                power = 0.8,
                                                alpha = 0.05,
                                                beta1 = 0.1,
                                                beta2 = 0.1,
                                                varY1 = 0.23,
                                                varY2 = 0.25,
                                                rho01 = 0.025,
                                                rho02 = 0.025,
                                                rho2 = 0.05,
                                                r = 1)$`Final m`[2],
                            alpha = 0.05,
                            beta1 = 0.1,
                            beta2 = 0.1,
                            varY1 = 0.23,
                            varY2 = 0.25,
                            rho01 = 0.025,
                            rho02 = 0.025,
                            rho2 = 0.05,
                            r = 1)$`Final Treatment (K)`[2]
  expect_true(result >= 15 & result <= 16,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# `m` and `K` align for D/AP
test_that("P-value m aligns with K for D/AP", {
  result <- calc_K_pval_adj(power = 0.8,
                            m = calc_m_pval_adj(K = 15,
                                                power = 0.8,
                                                alpha = 0.05,
                                                beta1 = 0.1,
                                                beta2 = 0.1,
                                                varY1 = 0.23,
                                                varY2 = 0.25,
                                                rho01 = 0.025,
                                                rho02 = 0.025,
                                                rho2 = 0.05,
                                                r = 1)$`Final m`[3],
                            alpha = 0.05,
                            beta1 = 0.1,
                            beta2 = 0.1,
                            varY1 = 0.23,
                            varY2 = 0.25,
                            rho01 = 0.025,
                            rho02 = 0.025,
                            rho2 = 0.05,
                            r = 1)$`Final Treatment (K)`[3]
  expect_true(result >= 15 & result <= 16,
              info = "The result should be 15, but with rounding between 15 and 16")
})
