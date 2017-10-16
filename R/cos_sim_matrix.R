#' Creat a mutational matrix with cosine similarities
#' 
#' Calculates the pairwise cosine similarity between the mutational profiles provided in the two mutation count matrices. 
#' The cosine similarity is a value between 0 (distinct) and 1 (identical) and indicates how much two vectors are alike.
#' 
#' @param mut_matrix1 96 mutation count matrix (dimensions: 96 mutations X n samples)
#' @param mut_matrix2 96 mutation count matrix (dimensions: 96 mutations X m samples)
#' @return Matrix with pairwise cosine similarities (dimensions: n mutational profiles X m mutational profiles)
#' 
#' @usage 
#' cos_sim_matrix(mut_matrix, cancer_signatures)
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
#' cos_sim_matrix(mut_mat, cancer_signatures)
#' 
#' @seealso
#' \code{\link{mut_matrix}},
#' \code{\link{fit_to_signatures}},
#' \code{\link{plot_cosine_heatmap}}
#' 
#' @export


cos_sim_matrix = function(mut_matrix1, mut_matrix2)
{
  n_samples1 = ncol(mut_matrix1)
  n_samples2 = ncol(mut_matrix2)
  res_matrix = matrix(nrow = n_samples1, ncol = n_samples2)
  
  for(s in 1:n_samples1)
  {
    signal1 = mut_matrix1[,s]
    cos_sim_vector = c()
    for(i in 1:n_samples2)
    {
      signal2 = mut_matrix2[,i]
      cos_sim_vector[i] = cos_sim(signal1, signal2)
    }
    res_matrix[s,] = cos_sim_vector
  }
  rownames(res_matrix) = colnames(mut_matrix1)
  colnames(res_matrix) = colnames(mut_matrix2)
  
  return(res_matrix)
}

explained_by_signatures = function(mut_matrix, signatures)
{
  .Defunct("explained_by_signatures", package="MutationalPatterns",
           msg=paste("This function has been renamed to",
                     "'cos_sim_matrix'."))
}