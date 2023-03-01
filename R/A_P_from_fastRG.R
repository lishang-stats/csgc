#' Extract adjacency matrix and probability matrix from fastRG package
#' @description Only applicable for stochastic block model and degree corrected
#' stochastic block model. Probability matrix (P) is obtained by maximising likelihood.
#'
#' @param model can be objects from fastRG package including: "undirected_dcsbm",
#' "directed_dcsbm", "undirected_sbm"
#'
#' @return adjacency matrix (A), probability matrix (P) and model type ("poisson" or "binomial")
#' @importFrom fastRG sample_sparse
#' @export
#'
#' @examples
#' set.seed(27)
#' k <- 5
#' n <- 100
#' B <- matrix(stats::runif(k * k), nrow = k, ncol = k)
#' theta <- round(stats::rlnorm(n, 2))
#' pi <- c(1, 2, 4, 1, 1)
#' m1 <- fastRG::dcsbm(theta = theta, B = B, pi = pi, expected_degree = 50)
#' out = A_P_from_fastRG(model = m1)
#' out$A
#' out$P
#' out$modeltype
A_P_from_fastRG <- function(model){
  if (length(grep("_dcsbm", attr(model,"class")))>0){
    A = as.matrix(sample_sparse(model))
    diag(A) = 0
    A[A>1] = 1
    z = as.numeric(model$z)
    P = if (length(grep("_sbm", attr(model,"class")))>0) sbm_mle(A,z)$P else dcsbm_mle(A,z)$P
    modeltype = if (m1$poisson_edges==TRUE) "poisson" else "binomial"
  } else{
    stop("The model does not belong to SBM/DCSBM, cannot calculate csgc statistics!")
  }
  return(list(A=A,P=P,modeltype=modeltype))
}
