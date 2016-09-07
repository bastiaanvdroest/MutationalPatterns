#' Plot signature contribution
#' 
#' Plot contribution of signatures
#' 
#' @param contribution Signature contribution matrix
#' @param signatures Signature matrix
#' @param index optional sample subset parameter
#' @param coord_flip Flip X and Y coordinates, default = FALSE
#' @param mode "relative" or "absolute"; to plot the relative contribution or absolute number of mutations, default = "relative"
#' @return Stacked barplot with contribution of each signatures for each sample
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_fill_discrete
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @export
#' 

plot_contribution = function(contribution, signatures, index=c(), coord_flip = FALSE, mode = "relative")
{
  # check mode parameter
  if(!(mode == "relative" | mode == "absolute")){stop("mode parameter should be either 'relative' or 'absolute' ")}
  # optional subsetting if index parameter is provided
  if(length(index > 0)){contribution = contribution[,index]}

  # These variables will be available at run-time, but not at compile-time.
  # To avoid compiling trouble, we initialize them to NULL.
  Sample = NULL
  Contribution = NULL
  Signature = NULL

  # if mode is relative
  if(mode == "relative")
  {
    # Plot contribution
    m_contribution = melt(contribution)
    colnames(m_contribution) = c("Signature", "Sample", "Contribution")
    
    plot = ggplot(m_contribution, aes(x = factor(Sample), y = Contribution, fill = factor(Signature), order = Sample)) + 
      geom_bar(position = "fill", stat="identity", colour="black")  +  
      # make sure sample ordering is correct
      xlim(rev(levels(factor(m_contribution$Sample)))) +
      # ylabel
      labs(x = "", y = "Relative contribution") +  
      scale_fill_discrete(name="Signature") +
      # white background
      theme_bw() +
      # no gridlines
      theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) +
      theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank())
  }

  # if mode is absolute
  if(mode == "absolute")
  {
    if(missing(signatures)){stop("For contribution plotting in mode 'absolute': also provide signatures matrix")}
    # total number of mutations per siganture
    total_signatures = colSums(signatures) 
    # calculate signature contribution in absolute number of signatures
    abs_contribution = contribution * total_signatures
    
    # Plot contribution
    m_contribution = melt(abs_contribution)
    colnames(m_contribution) = c("Signature", "Sample", "Contribution")

    plot = ggplot(m_contribution, aes(x = factor(Sample), y = Contribution, fill = factor(Signature), order = Sample)) + 
      geom_bar(stat="identity", colour = "black")  +  
      # make sure sample ordering is correct
      xlim(rev(levels(factor(m_contribution$Sample)))) +
      # ylabel
      labs(x = "", y = "Absolute contribution \n (no. mutations)") +  
      scale_fill_discrete(name="Signature") +
      # white background
      theme_bw() +
      # no gridlines
      theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) +
      theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank())
  }
  
  # optional coordinate flipping
  if(coord_flip == TRUE){plot = plot + coord_flip()}
  return(plot)
}
