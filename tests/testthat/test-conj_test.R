# Testing script the Conjunctive IU Test functions
# `calc_pwr_conj_test()`  - calculates power
# `calc_K_conj_test()`    - calculates # of clusters per treatment group
# `calc_m_conj_test()`    - calculates cluster size

# `calc_pwr_conj_test()` ----------------------------------------------------

# Test for Conjunctive IU test power
test_that("Conjunctive IU test power works", {
  result <- calc_pwr_conj_test(K = 15,
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
                               r = 1)
  expect_true(round(result, 3) == 0.899, info = paste(result))
})

# `calc_K_conj_test()` ------------------------------------------------------

# Test for Conjunctive IU test treatment group K
test_that("Conjunctive IU test treatment group 'K' works", {
  expect_equal(calc_K_conj_test(power = 0.8,
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
               12)
})

# Test for Conjunctive IU test control group K
test_that("Conjunctive IU test control group 'K' works", {
  expect_equal(calc_K_conj_test(power = 0.8,
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
               12)
})

# `calc_m_conj_test()` ------------------------------------------------------

# Test for Conjunctive IU test m
test_that("Conjunctive IU test 'm' works", {
  expect_equal(calc_m_conj_test(power = 0.8,
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
               86)
})


# Testing that all 3 functions align -------------------------------------------

# K aligns with power
test_that("Conjunctive IU test K calculation aligns with power", {
  result <- calc_K_conj_test(power = calc_pwr_conj_test(K = 15,
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
  expect_true(result >= 14.9 & result <= 16.1,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# K aligns with m
test_that("Conjunctive IU test K calculation aligns with power", {
  result <- calc_K_conj_test(power = 0.8,
                             m = calc_m_conj_test(power = 0.8,
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
  expect_true(result >= 14.9 & result <= 16.1,
              info = "The result should be 15, but with rounding between 15 and 16")
})

# Power aligns with K
test_that("Conjunctive IU test K calculation aligns with power", {
  result <- calc_m_conj_test(power = calc_pwr_conj_test(m = 300,
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
  expect_true(result >= 298 & result <= 303,
              info = "The result should be 300, but with rounding between 298 and 303")
})

# Check that it aligns with the source code ------------------------------------

# test_that("Conjunctive IU power calculation aligns with GitHub source code", {
#
#   devtools::source_url("https://github.com/siyunyang/coprimary_CRT/blob/main/powerSampleCal_varCluster_ttest.R?raw=TRUE")
#
#   result1 <- round(calc_pwr_conj_test(m = 300,
#                                 K = 15,
#                                 alpha = 0.05,
#                                 beta1 = 0.1,
#                                 beta2 = 0.1,
#                                 varY1 = 0.23,
#                                 varY2 = 0.25,
#                                 rho01 = 0.025,
#                                 rho02 = 0.025,
#                                 rho1 = 0.01,
#                                 rho2 = 0.05,
#                                 r = 1), 3)
#
#   result2 <- round(calPower_ttestIU(betas = c(0.1, 0.1),
#                               deltas = c(0, 0),
#                               vars = c(0.23, 0.25),
#                               rho01 = matrix(c(0.025, 0.01,
#                                                0.01, 0.025),
#                                              2, 2),
#                               rho2 = matrix(c(1, 0.05,
#                                               0.05, 1),
#                                             2, 2),
#                               N = 30,
#                               r = .5,
#                               m = 300,
#                               K = 2,
#                               alpha = 0.05), 3)
#
#   expect_true(result1 == result2,
#               info = paste0("Does not align with the original source code (", result1, " vs. ", result2, ")"))
# })
