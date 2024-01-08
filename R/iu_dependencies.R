#' Calculate power given design configurations based on intersection-union test, written by Siyun Yang
#'
#' @description
#' Internal function to calculate power given design configurations based on the intersection-union test, written by Siyun Yang with more details available at https://github.com/siyunyang/coprimary_CRT/blob/main/powerSampleCal_varCluster_ttest.R
#'
#' @returns A number.
#' @keywords internal
calPower_ttestIU <- function(betas, deltas, vars, rho01, rho2, N, r, m, K, alpha){
  # Variance of trt assignment
  sigmaz.square <- r*(1-r)

  # Define function to calculate correlation between betas
  calCovbetas <- function(vars,rho01,rho2, sigmaz.square, m, K){
    rho0k <- diag(rho01)
    sigmak.square <-(1+(m-1)*rho0k)*vars/(m*sigmaz.square)
    covMatrix <- diag(sigmak.square)
    for(row in 1:K )
    {
      for(col in 1:K)
      {
        if(row != col){
          covMatrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*(rho2[row,col]+(m-1)*rho01[row,col])/(m*sigmaz.square)
        }
      }
    }
    return(covMatrix)
  }

  # Define function to calculate correlation between test statistics
  calCorWks <-  function(vars,rho01,rho2, sigmaz.square, m, K)
  {
    top <- calCovbetas(vars,rho01,rho2, sigmaz.square, m, K)
    wCor <- diag(K)
    for(row in 1:K )
    {
      for(col in 1:K)
      {
        if(row != col){
          wCor[row,col] <- top[row,col]/sqrt(top[row,row]*top[col,col])
        }
      }
    }
    return(wCor)
  }
  sigmaks.sq <- diag(calCovbetas(vars,rho01,rho2, sigmaz.square, m, K))
  meanVector <- sqrt(N)*(betas-deltas)/sqrt(sigmaks.sq)
  wCor <- calCorWks(vars,rho01,rho2, sigmaz.square, m, K)
  criticalValue <- qt(p=(1-alpha), df=(N-2*K))
  pred.power <- pmvt(lower = rep(criticalValue,K),upper=rep(Inf,K),df = (N-2*K) , sigma = wCor,delta=meanVector)[1]
  return(pred.power)
} # End calPower_ttestIU()






#' Calculate sample size given design configurations based on intersection-union test, written by Siyun Yang
#'
#' @description
#' Internal function to calculate sample size given design configurations based on the intersection-union test, written by Siyun Yang with more details available at https://github.com/siyunyang/coprimary_CRT/blob/main/powerSampleCal_varCluster_ttest.R
#'
#' @returns A number.
#' @keywords internal
calSampleSize_ttestIU <- function(betas, deltas, vars, rho01, rho2, m, r, K, alpha, power)
{
  lowerBound <- 1
  upperBound <- 1000
  repeat{
    middle <- floor((lowerBound+upperBound)/2)
    power_temp <- calPower_ttestIU(betas,deltas,vars,rho01,rho2,N=middle,r,m,K,alpha)
    if(power_temp < power)
    {
      lowerBound <- middle
    }
    if(power_temp > power)
    {
      upperBound <- middle
    }
    if(power_temp == power)
    {
      return(middle)
      break
    }
    if((lowerBound-upperBound) == -1)
    {
      return(upperBound)
      break
    }
  }
}
