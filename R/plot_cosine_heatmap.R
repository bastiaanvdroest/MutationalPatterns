#' Plot cosine similarity heatmap
#' 
#' Plot cosine similarity between mutational profiles and signatures in a heatmap. 
#' 
#' 
#' @param explained Matrix with pairwise cosine similarities (dimensions: n samples X n signatures). Result from \code{\link{explained_by_signatures}}
#' @param sig_order Character vector with the desired order of the signature names for plotting. Optional.
#' @param cluster_samples Hierarchically cluster samples based on eucledian distance. Default = T.
#' @param method The agglomeration method to be used for hierarchical clustering. This should be one of 
#' "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) 
#' or "centroid" (= UPGMC). Default = "complete".
#'
#' @return Heatmap with cosine similarity of each signature with each 96 mutational profile
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom ggdendro dendro_data segment theme_dendro
#' @importFrom cowplot plot_grid
#'
#' @usage
#' plot_cosine_heatmap(explained_matrix, sig_order, cluster_samples = T)
#' 
#' @details 
#' The cosine similarity is a value between 0 (distinct) and 1 (identical) and indicates how much of the 96 
#' mutation profile can be explained by an individual signature. Similar mutational signatures will be equally 
#' good at explaining a mutational profile. For this reason it is recommended to cluster the signatures based 
#' on cosine similarity with \code{\link{cluster_signatures}} and plot the signatures in this order, see example.
#'
#' @examples
#' 
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
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                     package="MutationalPatterns"))
#'
#' ## Calculate the cosine similarity between each signature and each 96 mutational profile
#' explained = explained_by_signatures(mut_mat, cancer_signatures)
#' 
#' ## Cluster signatures based on cosine similarity 
#' sig_hclust = cluster_signatures(cancer_signatures)
#' sig_order = colnames(cancer_signatures)[sig_hclust$order]
#' 
#' ## Plot the cosine similarity between each signature and each sample with hierarchical 
#' ## sample clustering and signatures order based on similarity
#' 
#' plot_cosine_heatmap(explained, sig_order, cluster_samples = T)
#' 
#' 
#' @seealso
#' \code{\link{mut_matrix}},
#' \code{\link{explained_by_signatures}},
#' \code{\link{cluster_signatures}} 
#' 
#' @export

plot_cosine_heatmap = function(explained_matrix, sig_order, cluster_samples = T, method = "complete")
{
  # check explained argument
  if(class(explained_matrix) != "matrix")
  {stop("explained_matrix must be a matrix")}
  # if no signature order is provided, use the order as in the input matrix
  if(missing(sig_order))
  {
    sig_order = colnames(explained_matrix)
  }
  # check sig_order argument
  if(class(sig_order) != "character")
  {stop("sig_order must be a character vector")}
  if(length(sig_order) != ncol(explained_matrix))
  {stop("sig_order must have the same length as the number of signatures in the explained matrix")}
  
  # if cluster samples is TRUE, perform clustering
  if(cluster_samples == T)
  {
  # cluster samples based on eucledian distance between relative contribution
  hc.sample = hclust(dist(explained_matrix), method = method)
  # order samples according to clustering
  sample_order = rownames(explained_matrix)[hc.sample$order]
  }
  else
  {
    sample_order = rownames(explained_matrix)
  }
  
  # melt
  explained_matrix.m = melt(explained_matrix)
  # assign variable names
  colnames(explained_matrix.m) = c("Sample", "Signature", "Cosine.sim")
  
  # change factor levels to the correct order for plotting
  explained_matrix.m$Signature = factor(explained_matrix.m$Signature, levels = sig_order)
  explained_matrix.m$Sample = factor(explained_matrix.m$Sample, levels = sample_order)
  # plot heatmap
  heatmap = ggplot(explained_matrix.m, aes(x=Signature, y=Sample, fill=Cosine.sim)) + 
    geom_tile(color = "white") +
    scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Cosine \nsimilarity", limits = c(0,1)) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x=NULL, y=NULL)
  
  # if cluster samples is TRUE, make dendrogram
  if(cluster_samples == T)
  {
    # get dendrogram
    dhc = as.dendrogram(hc.sample)
    # rectangular lines
    ddata = dendro_data(dhc, type = "rectangle")
    # plot dendrogram of hierachical clustering
    dendrogram = ggplot(segment(ddata)) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
      coord_flip() + 
      scale_y_reverse(expand = c(0.2, 0)) + 
      theme_dendro()
    # combine plots
    plot_final = plot_grid(dendrogram, heatmap, align='h', rel_widths=c(0.3,1))
  }
  else
  {
    plot_final = heatmap
  }
  
  return(plot_final)
}