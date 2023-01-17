#' Return the trace of a matrix
#'
#' @param A matrix
#'
#' @return trace of matrix A
#' @export
#'
#' @examples
#' library(csgc)
#' A = matrix(c(1,2,3,4),2,2)
#' trace(A)
trace = function(A) {
  # Return trace of a matrix
  return(sum(diag(A)))
}
