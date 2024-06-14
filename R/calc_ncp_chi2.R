#' Find the non-centrality parameter corresponding to Type I error rate and statistical power
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
#' Allows user to find the corresponding non-centrality parameter for power analysis based on the Type I error rate, statistical power, and degrees of freedom.
#'
#' @param alpha Type I error rate; numeric.
#' @param power Desired statistical power in decimal form; numeric.
#' @param df Degrees of freedom; numeric.
#' @returns A number.
#' @examples
#' calc_ncp_chi2(alpha = 0.05, power = 0.8, df = 1)
#' @export
calc_ncp_chi2 <- function(alpha, power, df = 1){

  # First compute critical value based on alpha level and DF
  crit_val <- qchisq(alpha, df = df, lower.tail = FALSE)

  # Start with a small value and increase until it results in sufficient power
  ncp_start <- 1
  while(pchisq(crit_val, df = df, lower.tail = FALSE, ncp = ncp_start) < power){
    ncp_start <- ncp_start + 0.01
  }

  # Find exact NCP value to return
  ncp_out <- uniroot(function(ncp)
    return(pchisq(crit_val, df = df, lower.tail = FALSE, ncp = ncp) - power),
    c(0, ncp_start))$root

  # Return exact NCP value
  return(ncp_out)
} # End calc_ncp_chi2()
