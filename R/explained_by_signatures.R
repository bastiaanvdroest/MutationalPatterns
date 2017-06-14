#' Determine how much of a mutational profile can be explained by a mutational signature.
#' 
#' Calculates the cosine similarity between each mutation profile and mutational signature in the input matrices. 
#' The cosine similarity is a value between 0 (distinct) and 1 (identical) and indicates how much of the 96 
#' mutation profile can be explained by an individual signature.   
#' 
#' @param mut_matrix 96 mutation count matrix (dimensions: 96 mutations X n samples)
#' @param signatures Signature matrix (dimensions: 96 mutations X n signatures)
#' @return Matrix with pairwise cosine similarities (dimensions: n samples X n signatures)
#' 
#' @usage 
#' explained_by_signatures(mut_matrix, cancer_signatures)
#'
#' @examples
#' ## You can download the signatures from the pan-cancer study by
#' ## Alexandrov et al:
#' # http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt
#' 
#' ## We copied the file into our package for your convenience.
#' 
#' filename <- system.file("extdata/signatures_probabilities.txt",
#'                         package="MutationalPatterns")
#' cancer_signatures <- read.table(filename, sep = "\t", header = TRUE)
#' ## Reorder the columns to make the order of the trinucleotide changes compatible.
#' cancer_signatures <- cancer_signatures[order(cancer_signatures[,1]),]
#' ## Include only the signatures in the matrix.
#' cancer_signatures <- as.matrix(cancer_signatures[,4:33])
#'
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                     package="MutationalPatterns"))
#'
#' ## Calculate the cosine similarity between each signature and each 96 mutational profile
#' 
#' explained_by_signatures(mut_mat, cancer_signatures)
#' 
#' @seealso
#' \code{\link{mut_matrix}},
#' \code{\link{fit_to_signatures}},
#' \code{\link{plot_cosine_heatmap}}
#' 
#' @export


explained_by_signatures = function(mut_matrix, signatures)
{
  n_samples = ncol(mut_matrix)
  n_sigs = ncol(signatures)
  explained_matrix = matrix(nrow = n_samples, ncol = n_sigs)
  
  for(s in 1:n_samples)
  {
    signal = mut_matrix[,s]
    cos_sim_vector = c()
    for(i in 1:n_sigs)
    {
      signature = signatures[,i]
      cos_sim_vector[i] = cos_sim(signature, signal)
    }
    explained_matrix[s,] = cos_sim_vector
  }
  
  colnames(explained_matrix) = colnames(signatures)
  rownames(explained_matrix) = colnames(mut_matrix)
  return(explained_matrix)
}