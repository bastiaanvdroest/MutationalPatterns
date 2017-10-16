#' Cosine similarity function
#' 
#' Calculate the cosine similarity between two vectors
#' 
#' @param x Vector 1
#' @param y Vector 2
#' @return Cosine similarity value
#' 
#' @example 
#' 
#' 
#' @export


cos_sim = function(x, y)
{
  res = x %*% y / (sqrt(x %*% x) * sqrt(y %*% y))
  # coerce matrix to numeric
  res = as.numeric(res)
  return(res)
}