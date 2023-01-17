#' Calculate Maximum Likelihood Estimator of Stochastic Block Model
#' @description Calculate Maximum Likelihood Estimator of Stochastic Block Model,
#' using binomial model assuming simple graph without self-loop.
#' @param A adjacency matrix
#' @param z community labels
#'
#' @return estimated connection probability matrix (K) and estimated probability matrix (P)
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
#' A = gen_adj_sbm(K,z)$A
#' sbm_mle(A,z)

sbm_mle <- function(A,z){
  # INPUT
  # A = adjacency matrix
  # z = community labels
  #
  # OUTPUT
  # K:   estimated connection probability matrix
  # P:   estimated probability matrix
  k = length(unique(z))
  n = dim(A)[1]
  K = matrix(0,k,k)
  if (k==1){ K[1,1] = sum(A)*2 / ((n-1)*n) }
  else {
    for (i in 1:(k-1)){
      for (j in (i+1):k){
        K[i,j] = K[j,i] = sum(A[which(z==i),which(z==j)]) / (length(
          which(z==i))*length(which(z==j)))
      }
    }
    for (i in 1:k){
      K[i,i] = ifelse(length(which(z==i))>1, sum(A[which(z==i),which(z==i)]) / (length(
        which(z==i))*(length(which(z==i))-1)), 0)
    }
  }
  P = K[z,z]
  return(list(K=K,P=P))
}
