#' Convert dist and matrix type data into dataframe
#'
#' @param inDist a dist or matrix type data
#'
#' @return a dataframe
#' @export
#'
#' @examples
dist2df <- function(inDist) {
  if (class(inDist) == "matrix") inDist <- as.dist(inDist) # convert Matrix to to dist class
  if (class(inDist) != "dist") stop("wrong input type, should be a matrix or dist type data")
  A <- attr(inDist, "Size")
  B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
  if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
  if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
  data.frame(
    dist.x = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
    dist.y = rep(B[-length(B)], (length(B)-1):1),
    value = as.vector(inDist))
}
