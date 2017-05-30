#' Plot signature contribution heatmap
#' 
#' Plot relative contribution of signatures in a heatmap
#' 
#' @param contribution Signature contribution matrix
#' @param sig_order Character vector with the desired order of the signature names for plotting. Optional.
#' @param cluster_samples Hierarchically cluster samples based on eucledian distance. Default = T.
#'
#' @return Heatmap with relative contribution of each signature for each sample
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom ggdendro dendro_data segment theme_dendro
#' @importFrom cowplot plot_grid
#'
#' @usage
#' plot_contribution_heatmap(contribution, sig_order, cluster_samples = T)
#'
#' @examples
#'#' ## See the 'mut_matrix()' example for how we obtained the following
#' ## mutation matrix.
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                                 package="MutationalPatterns"))
#'
#' ## Extracting signatures can be computationally intensive, so
#' ## we use pre-computed data generated with the following command:
#' 
#' # nmf_res <- extract_signatures(mut_mat, rank = 2)
#'
#' nmf_res <- readRDS(system.file("states/nmf_res_data.rds",
#'                     package="MutationalPatterns"))
#'
#' ## Set signature names as row names in the contribution matrix
#' rownames(nmf_res$contribution) = c("Signature A", "Signature B")
#' 
#' ## Define signature order for plotting
#' sig_order = c("Signature B", "Signature A")
#'
#'
#' ## Contribution heatmap with signatures in defined order
#' plot_contribution_heatmap(nmf_res$contribution, 
#'                           sig_order = c("Signature B", "Signature A"))
#' 
#' ## Contribution heatmap without sample clustering
#' plot_contribution_heatmap(nmf_res$contribution, 
#'                           sig_order = c("Signature B", "Signature A"), 
#'                           cluster_samples = F)
#'
#' @seealso
#' \code{\link{extract_signatures}},
#' \code{\link{mut_matrix}},
#' \code{\link{plot_contribution}}
#'
#' @export

# plotting function for relative contribution of signatures in heatmap
plot_contribution_heatmap = function(contribution, sig_order, cluster_samples = T)
{
  # check contribution argument
  if(class(contribution) != "matrix")
    {stop("contribution must be a matrix")}
  # if no signature order is provided, use the order as in the input matrix
  if(missing(sig_order))
  {
    sig_order = rownames(contribution)
  }
  # check sig_order argument
  if(class(sig_order) != "character")
    {stop("sig_order must be a character vector")}
  if(length(sig_order) != nrow(contribution))
    {stop("sig_order must have the same length as the number of signatures in the contribution matrix")}
    
  # transpose
  contribution = t(contribution)
  # relative contribution
  contribution_norm = contribution / rowSums(contribution)
  
  # if cluster samples is TRUE, perform clustering
  if(cluster_samples == T)
  {
    # hiearchically cluster samples based on eucledian distance between relative contribution
    hc.sample = hclust(dist(contribution_norm))
    # order of samples according to hierarchical clustering
    sample_order = rownames(contribution)[hc.sample$order]
  } 
  else
  {
    sample_order = rownames(contribution)
  }

  # melt data frame
  contribution_norm.m = melt(contribution_norm)
  # assign variable names
  colnames(contribution_norm.m) = c("Sample", "Signature", "Contribution")
  # change factor levels to the order for plotting
  contribution_norm.m$Sample = factor(contribution_norm.m$Sample, levels = sample_order)
  contribution_norm.m$Signature = factor(contribution_norm.m$Signature, levels = sig_order)
  
  # plot heatmap
  heatmap = ggplot(contribution_norm.m, aes(x=Signature, y=Sample, fill=Contribution)) + 
    geom_tile(color = "white") +
    scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Relative \ncontribution", limits = c(0,1) ) +
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
    plot_final = cowplot::plot_grid(dendrogram, heatmap, align = 'h', rel_widths = c(0.3, 1))
  } 
  else
  {
    plot_final = heatmap
  }

  return(plot_final)
}