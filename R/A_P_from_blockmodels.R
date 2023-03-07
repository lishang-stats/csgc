#' Extract adjacency and probability matrix from blockmodels object
#'
#' @param model a blockmodels object
#'
#' @return adjacency matrix (A), probability matrix (P) and model type
#' ("bernoulli" or "poisson")
#' @export
#'
#' @examples
#' library(blockmodels)
#' npc <- 30 # nodes per class
#' Q <- 3 # classes
#' n <- npc * Q # nodes
#' Z<-diag(Q)%x%matrix(1,npc,1)
#' P<-matrix(runif(Q*Q),Q,Q)
#' M<-1*(matrix(runif(n*n),n,n)<Z%*%P%*%t(Z)) ## adjacency matrix
#' fit <- BM_bernoulli("SBM",M, plotting='')
#' out = A_P_from_blockmodels(model=fit)
#' out$A
#' out$P
#' out$modeltype
A_P_from_blockmodels <- function(model){
  if (model$membership_name %in% c("SBM", "SBM_sym") &
      model$model_name %in% c("bernoulli", "poisson")){
    A = model$adj
    if (model$membership_name == "SBM"){
      A[upper.tri(A)] = pmax(A[upper.tri(A)], t(A)[upper.tri(A)])
      A[lower.tri(A)] = t(A)[lower.tri(A)]
      warning("Convert directed graph to undirected graph.")
    }
    sink('NUL'); model$estimate(); sink()
    k = which.max(model$ICL)
    z = apply(model$memberships[[k]]$Z,1,which.max)
    K = model$model_parameters[[k]]$pi
    P = K[z,z]
    modeltype = model$model_name
    return(list(A=A,P=P,modeltype=modeltype))
  } else{
    stop("Input model must be stochastic block model with distribution bernoulli
         or poisson!")
  }
}
