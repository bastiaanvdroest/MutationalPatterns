#' Plot strand bias
#'
#' A function to plot the absolute or relative strand bias
#' 
#' @param strand_bias_df the data frame containing strand bias data
#' @param mode Either "absolute" or "relative"
#' @param colors the colors to use in the plot
#' @return a plot of the strand bias data
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_alpha_discrete
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 theme_bw
#' @export

plot_strand_bias = function(strand_bias_df, mode = "relative", colors)
{
  if(missing(colors)){colors=COLORS6}
  if(mode == "relative")
  {
    plot = ggplot(strand_bias_df, aes(x=strand, y=relative_contribution, group=strand, fill=type, alpha=strand)) +
      geom_bar(stat="identity", position = "dodge", colour="black", cex=0.5) + 
      scale_fill_manual(values= colors) +
      scale_alpha_discrete(range = c(1, 0.4)) +
      ylab("Relative contribution within group") +
      facet_grid(group ~ type) +
      theme_bw() +
      scale_x_discrete(breaks=NULL) +
      xlab("")
  }
  if(mode == "absolute")
  {
    plot = ggplot(strand_bias_df, aes(x=strand, y=no_mutations, group=strand, fill=type, alpha=strand)) +
      geom_bar(stat="identity", position = "dodge", colour="black", cex=0.5) + 
      scale_fill_manual(values= colors) +
      scale_alpha_discrete(range = c(1, 0.4)) +
      ylab("Total number of mutations") +
      facet_grid(group ~ type) +
      theme_bw() +
      scale_x_discrete(breaks=NULL) +
      xlab("")
  }
  return(plot)
}
