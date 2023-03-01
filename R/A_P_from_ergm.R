#' Extract adjacency and probability matrix from ergm object
#'
#' @param model an ergm object
#'
#' @return adjacency matrix (A) and probability matrix (P)
#' @importFrom ergm ergm
#' @export
#'
#' @examples
#' library(ergm)
#' data("sampson")
#' fit <- ergm(samplike ~ edges + cycle(4,semi=TRUE))
#' out = A_P_from_ergm(model=fit)
#' out$A
#' out$P
A_P_from_ergm <- function(model){
  if (attr(model,"class") == "ergm"){
    A = as.matrix(get(model$formula[[2]]))
    if (get.network.attribute(samplike,"hyper")==TRUE){
      stop("Hypergraph cannot be handled by csgc package!")
    }
    if (get.network.attribute(samplike,"loops")==TRUE){
      diag(A) = 0
      message("Remove self-loops.")
    }
    if (get.network.attribute(samplike,"directed")==TRUE){
      A[upper.tri(A)] = pmax(A[upper.tri(A)], t(A)[upper.tri(A)])
      A[lower.tri(A)] = t(A)[lower.tri(A)]
      message("Convert directed graph to undirected graph.")
    }
    if (get.network.attribute(samplike,"multiple")==TRUE){
      A[A>1] = 1
      message("Convert multigraph to simple graph.")
    }
    P = predict(model,output="matrix")
  } else{
    stop("Input should be ergm object!")
  }
  return(list(A=A,P=P))
}
