#' Calculate Wasserstein distance between uniform distribution and target distribution
#' @description Calculate Wasserstein distance between uniform distribution on (0,1) and
#' the target distribution from 0 to 1 (for example a sequence of p-values)
#'
#' @param p a sequence of values from 0 to 1
#'
#' @return 2 times the Wasserstein distance
#' @export
#' @importFrom zoo rollapply na.locf
#' @examples
#' set.seed(123)
#' library(csgc)
#' p = runif(100)
#' wasserstein_uniform(p)
wasserstein_uniform <- function(p){
  # required package: zoo
  # INPUT
  # p: a sequence of p-value
  #
  # OUTPUT
  # area: 2*\int_0^1 |ecdf(p) - x| dx
  # 2 times the absolute difference between empirical cdf of p-values and discrete uniform
  # Note: area \in [0,1]
  #
  u = seq(0,1,1/length(p))
  s = sort(c(p,u))
  au = rollapply(s,2,diff) * rollapply(s,2,sum) / 2
  aph = rep(NA,length(au))
  aph[which(s %in% sort(c(0,p)))] = u
  aph = na.locf(aph,fromLast = F)
  ap = rollapply(s,2,diff)*aph
  area = 2*sum(abs(au-ap))
  return(area)
}
