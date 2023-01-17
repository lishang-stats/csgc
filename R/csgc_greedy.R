#' Csgc greedy algorithm
#' @description Use csgc greedy algorithm to obtain optimised community labels
#' and the corresponding chi-square statistics, under stochastic block model
#' @param A adjacency matrix
#' @param z0 initial community labels
#' @param var.structure allow choose of binomial/poisson model, default is binomial
#' @param parallel allow use parallel computing, default is FALSE
#'
#' @return initial community labels/chi-square statistic/csgc statistics, optimised
#' community labels/chi-square statistic/csgc statistics
#'
#' @importFrom purrr map
#' @import foreach
#' @import doParallel
#'
#' @export
#'
#' @examples
#' library(csgc)
#' library(randnet)
#' k = 2
#' n = 400
#' K = matrix(c(0.8, 0.5, 0.5, 0.1), k, k)
#' z = rep(1:k,each=n/k)
#' A = gen_adj_sbm(K,z)$A
#' z0 = reg.SP(A,k)$cluster
#' greedy(A,z0)
csgc_greedy <- function(A,z0,var.structure="binomial",parallel=F){
  # required package: purrr, foreach, doParallel
  # INPUT
  # A  = adjacency matrix
  # z0 = estimated labels using other methods (e.g. spectral clustering)
  # var.structure = allow "binomial" or "poisson" when using csgc function
  # parallel = whether to use multiple CPUs
  #
  # OUTPUT
  # zout     = adjusted labels after greedy algorithm
  # chisqout = corresponding chisq stats for csgc using zout
  k = length(unique(z0))
  n = dim(A)[1]
  zadj = z0
  stats0 = csgc(A,sbm_mle(A,z0)$P,var.structure)$t
  chisq0 = sum(stats0^2)
  chisqadj = chisq0
  chisqdiff = .01
  zadj_list = chisqadj_list = list() #save labels and chisq value at each step
  if (!parallel){
    while ( chisqdiff > 0 ){
      chisq_list = numeric()
      zhat_list = list()
      k_list = 1:k
      for (iter in 1:(n*(k-1))){
        r = ceiling(iter/(k-1))
        s = ifelse(iter%%(k-1)==0, iter%%(k-1)+(k-1), iter%%(k-1))
        zhat = zadj
        zhat[r] = k_list[k_list!=zhat[r]][s]
        zhat_list[[iter]] = zhat
        Phat = sbm_mle(A,zhat)$P
        chisq_list[iter] = sum((csgc(A,Phat,var.structure)$t)^2)
      }
      chisqdiff = max(chisqadj-chisq_list)
      update_ind = which.min(chisq_list)
      chisqadj = min(chisq_list)
      message("chisq stats: ", chisqadj, ", change location: ", ceiling(update_ind/(k-1)))
      zadj = zhat_list[[update_ind]]
      message("adjusted labels: ", zadj)
      zadj_list = append(zadj_list, list(zadj))
      chisqadj_list = append(chisqadj_list, chisqadj)
    }
  } else {
    while ( chisqdiff > 0 ){
      chisq_list = numeric()
      zhat_list = list()
      k_list = 1:k
      cl <- makeCluster(detectCores())
      registerDoParallel(cl)
      out <- foreach(iter=1:(n*(k-1)),
                     .noexport = c("A","k","chisq_list","zhat_list"),
                     .export = ls(environment(sbm_mle))) %dopar% {
                       res = list()
                       r = ceiling(iter/(k-1))
                       s = ifelse(iter%%(k-1)==0, iter%%(k-1)+(k-1), iter%%(k-1))
                       zhat = zadj
                       zhat[r] = k_list[k_list!=zhat[r]][s]
                       res$z = zhat_list[[iter]] = zhat
                       Phat = sbm_mle(A,zhat)$P
                       res$chi = chisq_list[iter] = sum((csgc(A,Phat,var.structure)$t)^2)
                       return(res)
                     }
      stopCluster(cl)
      zhat_list = map(out,1)
      chisq_list = unlist(map(out,2))
      chisqdiff = max(chisqadj-chisq_list)
      update_ind = which.min(chisq_list)
      chisqadj = min(chisq_list)
      message("chisq stats: ", chisqadj, ", change location: ", ceiling(update_ind/(k-1)))
      zadj = zhat_list[[update_ind]]
      message("adjusted labels: ", zadj)
      zadj_list = append(zadj_list, list(zadj))
      chisqadj_list = append(chisqadj_list, chisqadj)
    }
  }
  # retrieve updated labels from second last step
  zout = zadj_list[[length(zadj_list)-1]]
  chisqout = chisqadj_list[[length(chisqadj_list)-1]]
  Pout = sbm_mle(A,zout)$P
  stats = csgc(A,Pout,var.structure)$t
  ans = list(zin=z0, chisqin=chisq0,statsin=stats0,
             zout=zout, chisqout=chisqout, statsout=stats)
  names(ans) = c("zin","chisqin","statsin",
                 "zout", "chisqout", "statsout")
  return(ans)
}