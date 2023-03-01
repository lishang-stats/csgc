#' Extract adjacency and probability matrix from igraph object
#'
#' @description Only applicable for stochastic block model. Number of communities
#' is obtained from greedy optimization of modularity. Community labels are obtained
#' from spectral clustering algorithm.
#' @param graph an igraph object with name "Stochastic block model"
#'
#' @return adjacency matrix (A) and probability matrix (P)
#' @importFrom igraph as_adjacency_matrix cluster_fast_greedy
#' @export
#'
#' @examples
#' num = 1000
#' pm = matrix(c(.5, .1, .1, .5), 2, 2)
#' bs = c(300, 700)
#' g = igraph::sample_sbm(num, pm, bs)
#' out = A_P_from_igraph(graph = g)
#' out$A
#' out$P
A_P_from_igraph <- function(graph){
  if (g$name == "Stochastic block model"){
    A = as_adjacency_matrix(graph,sparse=FALSE)
    diag(A) = 0
    A[upper.tri(A)] = pmax(A[upper.tri(A)], t(A)[upper.tri(A)])
    A[lower.tri(A)] = t(A)[lower.tri(A)]
    kh = max(cluster_fast_greedy(graph)$membership)
    z = spectral_sbm(A,kh)
    P = sbm_mle(A,z)$P
  } else{
    stop("Input graph does not come from stochastic block model!")
  }
  return(list(A=A,P=P))
}
