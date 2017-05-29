#' Cosine similarity matrix function
#' 
#' Calculate the cosine similarity between all columns in a matrix
#' 
#' @param matrix Matrix
#' @noRd
#' @return Cosine similarity matrix

cos_sim_matrix = function(matrix)
{
  mat = t(matrix)
  n = nrow(mat)
  res = matrix(nrow = n, ncol = n)
  
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      res[i,j] = cos_sim(mat[i,], mat[j,])
    }
  }
  rownames(res) = rownames(mat)
  colnames(res) = rownames(mat)
  return(res)
}