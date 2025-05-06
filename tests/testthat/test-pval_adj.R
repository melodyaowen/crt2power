# Testing script the P-Value Adjustment functions
# `calc_pwr_pval_adj()`  - calculates power
# `calc_K_pval_adj()`    - calculates # of clusters per treatment group
# `calc_m_pval_adj()`    - calculates cluster size


# Tests for `calc_pwer_pval_adj()` ---------------------------------------------

# Power based on Chi2 distribution
test_that("P-value final power works with Chi2 distribution", {
  expect_equal(calc_pwr_pval_adj(dist = "Chi2",
                                 K = 15,
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

# Power based on F-distribution
test_that("P-value final power works with F distribution", {
  expect_equal(calc_pwr_pval_adj(dist = "F",
                                 K = 15,
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
               c(0.8045, 0.8061, 0.8102))
})

# Throws error when distribution is invalid
test_that("P-value final power throws error correctly when distribution is invalid", {
  expect_error(calc_pwr_pval_adj(dist = 4,
                                 K = 15,
                                 m = 300,
                                 alpha = 0.05,
                                 beta1 = 0.1,
                                 beta2 = 0.1,
                                 varY1 = 0.23,
                                 varY2 = 0.25,
                                 rho01 = 0.025,
                                 rho02 = 0.025,
                                 rho2 = 0.05,
                                 r = 1))
})

# Tests for `calc_K_pval_adj()` ------------------------------------------------

# Final treatment group K for Chi2 distribution
test_that("P-value final treatment group 'K' works for Chi2 distribution", {
  expect_equal(calc_K_pval_adj(dist = "Chi2",
                               power = 0.8,
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

# Final control group K for Chi2 distribution
test_that("P-value final control group 'K' works for Chi2 distribution", {
  expect_equal(calc_K_pval_adj(dist = "Chi2",
                               power = 0.8,
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

# Final treatment group K for F distribution
test_that("P-value final treatment group 'K' works for F distribution", {
  expect_equal(calc_K_pval_adj(dist = "F",
                               power = 0.8,
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
               c(15, 15, 15))
})

# Final control group K for F distribution
test_that("P-value final control group 'K' works for F distribution", {
  expect_equal(calc_K_pval_adj(dist = "F",
                               power = 0.8,
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
               c(15, 15, 15))
})

# Throws error when power can't be reached for F distribution
test_that("P-value final treatment group 'K' throws error when power can't be reached", {
  expect_error(calc_K_pval_adj(dist = "F",
                               power = 0.9,
                               m = 30,
                               alpha = 0.001,
                               beta1 = 0.01,
                               beta2 = 0.01,
                               varY1 = 0.23,
                               varY2 = 0.25,
                               rho01 = 0.025,
                               rho02 = 0.025,
                               rho2 = 0.05,
                               r = 1))
})

# Tests for `calc_m_pval_adj()` ------------------------------------------------

# Final m for Chi2 distribution
test_that("P-value final 'm' works", {
  expect_equal(calc_m_pval_adj(dist = "Chi2",
                               power = 0.8,
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

# Final m for F distribution
test_that("P-value final 'm' works", {
  expect_equal(calc_m_pval_adj(dist = "F",
                               power = 0.8,
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
               c(275, 267, 248))
})

# Throws an error when 'm' cannot be large enough to attain desired power
test_that("P-value final 'm' throws error when desired power cannot be attained (F dist)", {
  expect_error(calc_m_pval_adj(dist = "F",
                               power = 0.9,
                               K = 4,
                               alpha = 0.05,
                               beta1 = 0.1,
                               beta2 = 0.1,
                               varY1 = 0.23,
                               varY2 = 0.25,
                               rho01 = 0.025,
                               rho02 = 0.025,
                               rho2 = 0.05,
                               r = 1))
})

test_that("P-value final 'm' throws error when desired power cannot be attained (Chi2 dist)", {
  expect_error(calc_m_pval_adj(dist = "Chi2",
                               power = 0.9,
                               K = 4,
                               alpha = 0.05,
                               beta1 = 0.1,
                               beta2 = 0.1,
                               varY1 = 0.23,
                               varY2 = 0.25,
                               rho01 = 0.025,
                               rho02 = 0.025,
                               rho2 = 0.05,
                               r = 1))
})

# Testing that all 3 functions align -------------------------------------------

# Power and `m` align for Bonferroni (Chi2 distribution)
test_that("P-value power aligns with m for Bonferroni (Chi2 dist)", {
  result <- calc_m_pval_adj(dist = "Chi2",
                            power = calc_pwr_pval_adj(dist = "Chi2",
                                                      K = 15,
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
  expect_equal(result, 300, tolerance = 5,
               info = "The result should be 300, but with rounding between 300 and 301")
})

# Power and `m` align for Sidak (Chi2 distribution)
test_that("P-value power aligns with m for Sidak (Chi2 dist)", {
  result <- calc_m_pval_adj(dist = "Chi2",
                            power = calc_pwr_pval_adj(dist = "Chi2",
                                                      K = 15,
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
  expect_equal(result, 300, tolerance = 5,
              info = "The result should be 300, but with rounding between 300 and 301")
})

# Power and `m` align for D/AP (Chi2 distribution)
test_that("P-value power aligns with m for D/AP (Chi2 dist)", {
  result <- calc_m_pval_adj(dist = "Chi2",
                            power = calc_pwr_pval_adj(dist = "Chi2",
                                                      K = 15,
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
  expect_equal(result, 300, tolerance = 5,
              info = "The result should be 300, but with rounding between 300 and 301")
})

# Power and `K` align for Bonferroni (Chi2 distribution)
test_that("P-value power aligns with K for Bonferroni (Chi2 Dist)", {
  result <- calc_K_pval_adj(dist = "Chi2",
                            power = calc_pwr_pval_adj(dist = "Chi2",
                                                      K = 15,
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
  expect_equal(result, 15, tolerance = 2,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# Power and `K` align for Sidak (Chi2 distribution)
test_that("P-value power aligns with K for Sidak (Chi2 dist)", {
  result <- calc_K_pval_adj(dist = "Chi2",
                            power = calc_pwr_pval_adj(dist = "Chi2",
                                                      K = 15,
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
  expect_equal(result, 15, tolerance = 2,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# Power and `K` align for D/AP (Chi2 distribution)
test_that("P-value power aligns with K for D/AP (Chi2 dist)", {
  result <- calc_K_pval_adj(dist = "Chi2",
                            power = calc_pwr_pval_adj(dist = "Chi2",
                                                      K = 15,
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
  expect_equal(result, 15, tolerance = 2,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# `m` and `K` align for Bonferroni (Chi2 distribution)
test_that("P-value m aligns with K for Bonferroni (Chi2 dist)", {
  result <- calc_K_pval_adj(dist = "Chi2",
                            power = 0.8,
                            m = calc_m_pval_adj(dist = "Chi2",
                                                K = 15,
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
  expect_equal(result, 15, tolerance = 2,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# `m` and `K` align for Sidak (Chi2 distribution)
test_that("P-value m aligns with K for Sidak (Chi2 dist)", {
  result <- calc_K_pval_adj(dist = "Chi2",
                            power = 0.8,
                            m = calc_m_pval_adj(dist = "Chi2",
                                                K = 15,
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
  expect_equal(result, 15, tolerance = 2,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# `m` and `K` align for D/AP (Chi2 distribution)
test_that("P-value m aligns with K for D/AP (Chi2 dist)", {
  result <- calc_K_pval_adj(dist = "Chi2",
                            power = 0.8,
                            m = calc_m_pval_adj(dist = "Chi2",
                                                K = 15,
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
  expect_equal(result, 15, tolerance = 2,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# Power and `m` align for Bonferroni (F distribution)
test_that("P-value power aligns with m for Bonferroni (F dist)", {
  result <- calc_m_pval_adj(dist = "F",
                            power = calc_pwr_pval_adj(dist = "F",
                                                      K = 15,
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
  expect_equal(result, 300, tolerance = 5,
              info = "The result should be 300, but with rounding between 300 and 301")
})

# Power and `m` align for Sidak (F distribution)
test_that("P-value power aligns with m for Sidak (F dist)", {
  result <- calc_m_pval_adj(dist = "F",
                            power = calc_pwr_pval_adj(dist = "F",
                                                      K = 15,
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
  expect_equal(result, 300, tolerance = 5,
              info = "The result should be 300, but with rounding between 300 and 301")
})

# Power and `m` align for D/AP (F distribution)
test_that("P-value power aligns with m for D/AP (F dist)", {
  result <- calc_m_pval_adj(dist = "F",
                            power = calc_pwr_pval_adj(dist = "F",
                                                      K = 15,
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
  expect_equal(result, 300, tolerance = 5,
              info = "The result should be 300, but with rounding between 300 and 301")
})

# Power and `K` align for Bonferroni (F distribution)
test_that("P-value power aligns with K for Bonferroni (F Dist)", {
  result <- calc_K_pval_adj(dist = "F",
                            power = calc_pwr_pval_adj(dist = "F",
                                                      K = 15,
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
  expect_equal(result, 15, tolerance = 2,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# Power and `K` align for Sidak (F distribution)
test_that("P-value power aligns with K for Sidak (F dist)", {
  result <- calc_K_pval_adj(dist = "F",
                            power = calc_pwr_pval_adj(dist = "F",
                                                      K = 15,
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
  expect_equal(result, 15, tolerance = 2,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# Power and `K` align for D/AP (F distribution)
test_that("P-value power aligns with K for D/AP (F dist)", {
  result <- calc_K_pval_adj(dist = "F",
                            power = calc_pwr_pval_adj(dist = "F",
                                                      K = 15,
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
  expect_equal(result, 15, tolerance = 2,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# `m` and `K` align for Bonferroni (F distribution)
test_that("P-value m aligns with K for Bonferroni (F dist)", {
  result <- calc_K_pval_adj(dist = "F",
                            power = 0.8,
                            m = calc_m_pval_adj(dist = "F",
                                                K = 15,
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
  expect_equal(result, 15, tolerance = 2,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# `m` and `K` align for Sidak (F distribution)
test_that("P-value m aligns with K for Sidak (F dist)", {
  result <- calc_K_pval_adj(dist = "F",
                            power = 0.8,
                            m = calc_m_pval_adj(dist = "F",
                                                K = 15,
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
  expect_equal(result, 15, tolerance = 2,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# `m` and `K` align for D/AP (F distribution)
test_that("P-value m aligns with K for D/AP (F dist)", {
  result <- calc_K_pval_adj(dist = "F",
                            power = 0.8,
                            m = calc_m_pval_adj(dist = "F",
                                                K = 15,
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
  expect_equal(result, 15, tolerance = 2,
              info = "The result should be 15, but with rounding between 15 and 16")
})

