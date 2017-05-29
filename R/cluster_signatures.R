#' Signature clustering
#' 
#' Hierarchical clustering of signatures based on cosine similarity
#' 
#' @param signatures Matrix with usually 96 mutational features (rows) and any number of signatures (columns)
#' @noRd
#' @return hclust object


cluster_signatures = function(signatures)
{
  # construct cosine similarity matrix
  sim = cos_sim_matrix(signatures)
  # transform to distance
  dist = as.dist(1 - sim)
  # perform hierarchical clustering
  hc_sig_cos = hclust(dist)
  return(hc_sig_cos)
}