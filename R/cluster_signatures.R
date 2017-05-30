#' Signature clustering function
#' 
#' Hierarchical clustering of signatures based on cosine similarity
#' 
#' @param signatures Matrix with usually 96 trinucleotides (rows) and any number of signatures (columns)
#' @return hclust object
#' 
#' @examples
#' ## You can download the signatures from the pan-cancer study by
#' ## Alexandrov et al:
#' # http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt
#' 
#' ## We copied the file into our package for your convenience.
#' filename <- system.file("extdata/signatures_probabilities.txt",
#'                         package="MutationalPatterns")
#'
#' cancer_signatures <- read.table(filename, sep = "\t", header = TRUE)
#'
#' ## Reorder the columns to make the order of the trinucleotide changes compatible.
#' cancer_signatures <- cancer_signatures[order(cancer_signatures[,1]),]
#'
#' ## Include only the signatures in the matrix.
#' cancer_signatures <- as.matrix(cancer_signatures[,4:33])
#'
#' ## Hierarchically luster the cancer signatures based on cosine similarity
#' hclust_cancer_signatures = cluster_signatures(cancer_signatures)
#' 
#' ## Plot dendrogram
#' plot(hclust_cancer_signatures)
#' 
#' ## Save the signature names in the order of the clustering
#' sig_order = colnames(cancer_signatures)[hclust_cancer_signatures$order]
#' 
#' @seealso
#' \code{\link{plot_contribution_heatmap}}
#' 
#' @export



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