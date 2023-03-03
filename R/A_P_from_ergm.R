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
#' data(faux.dixon.high)
#' fit <- ergm(faux.dixon.high ~ edges + mutual)
#' out = A_P_from_ergm(model=fit)
#' out$A
#' out$P
A_P_from_ergm <- function(model){
  message("Please check that the model satisfies the conditional independence of the edges.")
  if (attr(model,"class") == "ergm"){
    data = get(model$formula[[2]])
    A = as.matrix(data)
    if (dim(A)[1]<100){
      warning("Too few vertices, convergence is not guaranteed for csgc function!")
    }
    if (get.network.attribute(data,"hyper")==TRUE){
      stop("Hypergraph cannot be handled by csgc package!")
    }
    if (get.network.attribute(data,"loops")==TRUE){
      diag(A) = 0
      warning("Remove self-loops.")
    }
    if (get.network.attribute(data,"directed")==TRUE){
      A[upper.tri(A)] = pmax(A[upper.tri(A)], t(A)[upper.tri(A)])
      A[lower.tri(A)] = t(A)[lower.tri(A)]
      warning("Convert directed graph to undirected graph.")
    }
    if (get.network.attribute(data,"multiple")==TRUE){
      A[A>1] = 1
      warning("Convert multigraph to simple graph by collapsing edges.")
    }
    P = predict(model,output="matrix")
  } else{
    stop("Input should be ergm object!")
  }
  return(list(A=A,P=P))
}
