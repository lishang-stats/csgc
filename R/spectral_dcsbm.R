#' Spectral clustering for Degree Corrected Stochastic Block Model
#'
#' @param A adjacency matrix
#' @param k number of communities
#'
#' @return estimated labels
#' @export
#' @importFrom irlba irlba
#' @importFrom stats kmeans
#' @description Tianxi Li, Elizaveta Levina, Ji Zhu and Can M. Le (2022). randnet:
#' Random Network Model Estimation, Selection and Parameter Tuning. R
#' package version 0.5. https://CRAN.R-project.org/package=randnet
#' @examples
#' library(csgc)
#' k = 4
#' n = 200
#' prob = c(0.8, 0.1)
#' K = matrix(prob[2],k,k)
#' diag(K) = prob[1]
#' Theta = runif(n,0.8,1)
#' z = rep(1:k,each=n/k)
#' A = gen_adj_dcsbm(K,Theta,z)$A
#' spectral_dcsbm(A,k)
spectral_dcsbm <- function(A,k){
  # required package: irlba, stats
  # INPUT
  # A = adjacency matrix
  # k = number of communities
  #
  # OUTPUT
  # labels: estimated labels
  V = irlba(A+mean(colSums(A))/nrow(A),k)$v
  vr = apply(V,1,function(x) sqrt(sum(x^2)))
  reg.V = diag(1/vr) %*% V
  labels = kmeans(reg.V,k,nstart=30,iter.max=100)$cluster
  return(labels)
}
