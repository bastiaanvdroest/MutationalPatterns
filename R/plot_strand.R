#' Plot strand per base substitution type
#' 
#' @description For each base substitution type and transcriptional strand the total number of mutations
#' and the relative contribution within a group is returned
#' @param strand_bias_df Data.frame, result from strand_bias function
#' @param mode Either "absolute" for absolute number of mutations, or "relative" for relative contribution, default = "relative"
#' @param colors Optional color vector for plotting with 6 values
#' @return Barplot
#' @export


plot_strand = function(strand_bias_df, mode = "relative", colors)
{
  # if colors parameter not provided, set to default colors
  if(missing(colors)){colors=COLORS6}
  # plot relative contribution within each group
  if(mode == "relative")
  {
    plot = ggplot(strand_bias_df, aes(x=type, y=relative_contribution, fill=type, alpha=strand)) +
      geom_bar(stat="identity", position = "dodge", colour="black", cex=0.5) + 
      scale_fill_manual(values= colors) +
      scale_alpha_discrete(range = c(1, 0.4)) +
      ylab("Relative contribution") +
      facet_grid(. ~ group) +
      theme_bw() +
      scale_x_discrete(breaks=NULL) +
      xlab("")
  }
  # plot absolute contribution within each group
  if(mode == "absolute")
  {
    plot = ggplot(strand_bias_df, aes(x=type, y=no_mutations, fill=type, alpha=strand)) +
      geom_bar(stat="identity", position = "dodge", colour="black", cex=0.5) + 
      scale_fill_manual(values= colors) +
      scale_alpha_discrete(range = c(1, 0.4)) +
      ylab("Total number of mutations") +
      facet_grid(. ~ group) +
      theme_bw() +
      scale_x_discrete(breaks=NULL) +
      xlab("")
  }
  return(plot)
}