#' functions from MCMCpack to collapse and
#' re-expand symmetric matrices

vech = function (x) 
{
  x <- as.matrix(x)
  if (dim(x)[1] != dim(x)[2]) {
    stop("Non-square matrix passed to vech().\n")
  }
  output <- x[lower.tri(x, diag = TRUE)]
  dim(output) <- NULL
  return(output)
}

xpnd = function (x, nrow = NULL) 
{
  dim(x) <- NULL
  if (is.null(nrow)) 
    nrow <- (-1 + sqrt(1 + 8 * length(x)))/2
  output <- matrix(0, nrow, nrow)
  output[lower.tri(output, diag = TRUE)] <- x
  hold <- output
  hold[upper.tri(hold, diag = TRUE)] <- 0
  output <- output + t(hold)
  return(output)
}
