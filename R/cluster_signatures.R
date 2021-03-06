#' Signature clustering function
#' 
#' Hierarchical clustering of signatures based on cosine similarity
#' 
#' @param signatures Matrix of all variants from a mutation type (rows) and any number of
#' signatures (columns)
#' @param method     (Optional) The agglomeration method to be used for hierarchical
#' clustering. This should be one of "ward.D", "ward.D2", "single", "complete",
#' "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or
#' "centroid" (= UPGMC). Default = "complete".
#' @return hclust object
#' 
#' @examples
#' ## You can download mutational signatures from the COSMIC database:
#' # sp_url = http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt
#' # cancer_signatures = read.table(sp_url, sep = "\t", header = T)
#' 
#' ## We copied the file into our package for your convenience.
#' filename <- system.file("extdata/signatures_probabilities.txt",
#'                         package="MutationalPatterns")
#' cancer_signatures <- read.table(filename, sep = "\t", header = TRUE)
#'
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                     package="MutationalPatterns"))
#' 
#' ## Match the order to MutationalPatterns standard of mutation matrix
#' order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
#' ## Reorder cancer signatures dataframe
#' cancer_signatures = cancer_signatures[order,]
#' ## Use trinucletiode changes names as row.names
#' ## row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
#' ## Keep only 96 contributions of the signatures in matrix
#' cancer_signatures = as.matrix(cancer_signatures[,4:33])
#' ## Rename signatures to number only
#' colnames(cancer_signatures) = as.character(1:30)
#' 
#' ## Hierarchically cluster the cancer signatures based on cosine similarity
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

cluster_signatures = function(signatures, method = "complete")
{
    # construct cosine similarity matrix
    sim = cos_sim_matrix(signatures, signatures)
    # transform to distance
    dist = as.dist(1 - sim)
    # perform hierarchical clustering
    hc_sig_cos = hclust(dist, method = method)
    return(hc_sig_cos)
}
