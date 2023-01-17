#' Generate adjacency matrix from Degree Corrected Stochastic Block Model
#'
#' @param K connection probability matrix
#' @param Theta degree parameter vector
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
#' Theta = runif(n,0.8,1)
#' z = rep(1:k,each=n/k)
#' gen_adj_dcsbm(K,Theta,z)

gen_adj_dcsbm <- function(K,Theta,z){
  # INPUT
  # K = connection probability matrix
  # Theta = degree parameter vector
  # z = community labels
  #
  # OUTPUT
  # A:  adjacency matrix
  # P:  probability matrix
  k = dim(K)[1]
  n = length(z)
  P = diag(Theta)%*%K[z,z]%*%diag(Theta)
  diag(P) = 0
  A = matrix(0,n,n)
  A[upper.tri(A)] = rpois(n*(n-1)/2,P[upper.tri(P)])
  A = A+t(A)
  return(list(A=A,P=P))
}
