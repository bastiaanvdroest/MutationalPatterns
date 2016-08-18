#' Plot strand bias per base substitution type
#' 
#' @description For each base substitution type and transcriptional strand the total number of mutations
#' and the relative contribution within a group is returned
#' @param strand_bias Data.frame, result from strand_bias function
#' @param colors Optional color vector for plotting with 6 values
#' @return Barplot
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_alpha_discrete
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 position_dodge
#' @export


plot_strand_bias = function(strand_bias, colors)
{
  # if colors parameter not provided, set to default colors
  if(missing(colors)){colors=COLORS6}
  # plot strand bias with poisson test results
  plot = ggplot(strand_bias, aes(x=type, y=log2(ratio), fill=type)) +
    scale_fill_manual(values=COLORS6) +
    geom_bar(colour="black", stat="identity",position = "identity") +
    scale_y_continuous(limits=c(-0.7, 0.7), breaks=seq(-1, 1, 0.2)) +
    geom_text(aes(x = type, y = log2(ratio), ymax = log2(ratio), 
                  label=significant, vjust=ifelse(sign(log2(ratio)) > 0, 0.5, 1)), 
              size = 8, position = position_dodge(width=1)) +
    facet_grid(. ~ group) +
    theme_bw()  +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank(), legend.title=element_blank()) +
    xlab("") + 
    ylab("log2(transcribed/untranscribed)") +
    scale_x_discrete(breaks=NULL)
  return(plot)
}
