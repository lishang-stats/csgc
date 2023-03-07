#' Extract adjacency and probability matrix from SbmFit and DcSbmFit object
#'
#' @description Extract adjacency and probability matrix from SbmFit and DcSbmFit object
#' from "greed" package
#' @param data network connection data with type "matrix", "dgCmatrix" or "igraph"
#' @param blockmodel an SbmFit or DcSbmFit object
#'
#' @return adjacency matrix (A), probability matrix (P) and model type ("bernoulli"
#' or "poisson")
#' @export
#'
#' @importFrom igraph as_adj
#' @examples
#' library(greed)
#' data(Books)
#' sbm = greed(Books$X, model = Sbm())
#' out1 = A_P_from_greed(data=Books$X, blockmodel=sbm)

#' library(igraphdata)
#' data(karate)
#' dcsbm = greed(karate, model=DcSbm())
#' out2 = A_P_from_greed(data=karate, blockmodel=dcsbm)
A_P_from_greed <- function(data, blockmodel){
  if (class(data)[1] == "igraph"){
    A = as_adj(data,sparse=FALSE)
  } else{
    A = as.matrix(data)
  }
  if (class(blockmodel@model)[1] == "Sbm"){
    z = blockmodel@cl
    par = coef(blockmodel)
    P = par$thetakl[z,z]
    modeltype = "bernoulli"
  } else if (class(blockmodel@model)[1] == "DcSbm"){
    z = blockmodel@cl
    zc = blockmodel@obs_stats$counts
    par = coef(blockmodel)
    k = blockmodel@K
    K = matrix(0,k,k)
    for (i in 1:k){
      for (j in 1:k){
        K[i,j] = par$thetakl[i,j] * zc[i] * zc[j]
      }
    }
    if (prior(blockmodel)@type == "directed"){
      Thetain = as.vector(par$gammain)
      Thetaout = as.vector(par$gammaout)
      P = diag(Thetaout) %*% K[z,z] %*% diag(Thetain)
    } else{
      Theta = as.vector(par$gamma)
      P = diag(Theta) %*% K[z,z] %*% diag(Theta)
      modeltype = "poisson"
    }
  } else{
    stop("Input model must be SBM or DCSBM!")
  }
  return(list(A=A,P=P,modeltype=modeltype))
}
