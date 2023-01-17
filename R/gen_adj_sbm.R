#' Generate adjacency matrix from Stochastic Block Model
#'
#' @param K connection probability matrix
#' @param z community labels
#'
#' @return adjacency matrix (A) and probability matrix (P)
#' @export
#'
#' @examples
#' library(csgc)
#' k = 4
#' n = 200
#' prob = c(0.8, 0.1)
#' K = matrix(prob[2],k,k)
#' diag(K) = prob[1]
#' z = rep(1:k,each=n/k)
#' gen_adj_sbm(K,z)

gen_adj_sbm <- function(K,z){
  # INPUT
  # K = connection probability matrix
  # z = labels
  #
  # OUTPUT
  # A:  adjacency matrix
  # P:  probability matrix
  k = dim(K)[1]
  n = length(z)
  P = K[z,z]
  diag(P) = 0
  A = matrix(0,n,n)
  A[upper.tri(A)] = rbinom(n*(n-1)/2,1,P[upper.tri(P)])
  A = A+t(A)
  return(list(A=A,P=P))
}
