# @description
# `calc_ncp_chi2()`
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
}
