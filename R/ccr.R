#' Calculate Correct Classification Rate
#' @param label_true true community labels
#' @param label_est estimated community labels
#' @param parallel allow use parallel computing, default is FALSE
#'
#' @return Correct Classification Rate
#' @export
#'
#' @importFrom iterpc iterpc getall
#' @import foreach
#' @import doParallel
#' @import parallel
#'
#' @examples
#' library(csgc)
#' label_true = sample(1:5,100,T)
#' label_est = sample(1:5,100,T)
#' ccr(label_true,label_est)
ccr <- function(label_true, label_est, parallel=F){
  # required package: iterpc, foreach, doParallel
  # INPUT
  # label_true = true labels
  # label_est = estimated labels
  # parallel = whether to use multiple CPUs
  #
  # OUTPUT
  # ccr:   correct classification rate
  k = length(unique(label_true))
  n = length(label_true)
  I = iterpc(k,k,ordered = T)
  all_perm = getall(I)
  if (!parallel){
    ccr_all = numeric()
    label_est = factor(label_est)
    for (i in 1:factorial(k)){
      levels(label_est) = all_perm[i,]
      ccr_all[i] = sum(label_est==label_true)
    }
  } else {
    num_perm <- nrow(all_perm)
    member_mat <- matrix(0, n, k)
    for (i in 1:n) {
      member_mat[i, label_est[i]] <- 1
    }
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    ccr_all <- foreach(iter=1:num_perm,.combine ='rbind') %dopar% {
      perm <- all_perm[iter, ]
      member_mat_perm <- member_mat[, perm]
      label_perm <- apply(member_mat_perm, 1, function(x) which(x == 1))
      return(sum(label_perm==label_true))
    }
    stopCluster(cl)
  }
  ccr = max(ccr_all)/n
  return(ccr)
}
