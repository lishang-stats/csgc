#' Calculate Maximum Likelihood Estimator of Degree Corrected Stochastic Block Model
#' @description Calculate Maximum Likelihood Estimator of Degree Corrected Stochastic
#'  Block Model, using poisson model assuming simple graph without self-loop.
#' @param A adjacency matrix
#' @param z community labels
#'
#' @return estimated connection probability matrix (K) and estimated probability matrix (P)
#' @export
#'
#' @importFrom rootSolve multiroot
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
#' A = gen_adj_dcsbm(K,Theta,z)$A
#' dcsbm_mle(A,z)
dcsbm_mle <- function(A,z){
  # required package: rootSolve
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
  for (i in 1:k){
    for (j in i:k){
      K[i,j] = K[j,i] = sum(A[which(z==i),which(z==j)]) / (length(
        which(z==i))*length(which(z==j)))
    }
  }
  deg = rowSums(A)
  Theta0 = Theta = numeric(n)
  for (j in 1:k){
    Theta0[z==j] = ifelse(sum(deg[z==j])==0,0,deg[z==j]/sum(deg[z==j]))
  }
  f = function(x) {
    Thetasq_r = sum(x^2)
    ret = (k_r-kappa_r*x)*(1-Thetasq_r) - x*m_rr*Thetasq_r+m_rr*x^2
  }
  for (r in 1:k) {
    init = Theta0[z==r]
    k_r = deg[z==r]
    kappa_r = sum(k_r)
    m_rr = diag(K)[r]
    res = multiroot(f,init)
    Theta[z==r] = res$root
    diag(K)[r] = diag(K)[r] / (1-sum(res$root^2))
  }
  P = diag(Theta)%*%K[z,z]%*%diag(Theta)
  diag(P) = 0
  return(list(K=K,Theta=Theta,P=P))
}
