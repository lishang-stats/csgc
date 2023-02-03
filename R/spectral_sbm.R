#' Spectral clustering for Stochastic Block Model
#'
#' @param A adjacency matrix
#' @param k number of communities
#'
#' @return estimated labels
#' @export
#' @importFrom irlba irlba
#' @importFrom stats kmeans
#'
#' @examples
#' library(csgc)
#' k = 4
#' n = 200
#' prob = c(0.8, 0.1)
#' K = matrix(prob[2],k,k)
#' diag(K) = prob[1]
#' z = rep(1:k,each=n/k)
#' A = gen_adj_sbm(K,z)$A
#' spectral_sbm(A,k)
spectral_sbm <- function(A,k){
  # required package: irlba, stats
  # INPUT
  # A = adjacency matrix
  # k = number of communities
  #
  # OUTPUT
  # labels: estimated labels
  V = irlba(A+mean(colSums(A))/nrow(A),k)$v
  labels = kmeans(V,k)$cluster
  return(labels)
}
