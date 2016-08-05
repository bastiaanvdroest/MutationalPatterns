#' Plot strand bias per base substitution type
#' 
#' @description For each base substitution type and transcriptional strand the total number of mutations
#' and the relative contribution within a group is returned
#' @param strand_bias_df Data.frame, result from strand_bias function
#' @param mode Either "absolute" for absolute number of mutations, or "relative" for relative contribution, default = "relative"
#' @param colors Optional color vector for plotting with 6 values
#' @return Barplot
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
