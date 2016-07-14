#' 
#' 
#' @param strand_bias_df
#' @param mode 
#' @param colors 
#' @return
#' @importFrom ggplot2 ggplot
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
