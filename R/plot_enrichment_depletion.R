#' Plot enrichment depletion barplot
#' 
#' @param df Dataframe result from enrichment_depletion_test()
#' @return plot
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 element_blank
#' @export

plot_enrichment_depletion = function(df)
{
  plot = ggplot(data=df, aes(x=by, y=log2(observed/expected), fill=by)) + 
    geom_bar(colour="black" , stat="identity", position=position_dodge()) +
    scale_y_continuous(limits=c(-0.75, 0.75), breaks=seq(-1, 1, 0.25)) +
    geom_text(aes(x = by, y = log2(observed/expected), ymax = log2(observed/expected), 
                  label=significant, vjust=ifelse(sign(log2(observed/expected)) > 0, 0.5, 1)), 
              size = 8, position = position_dodge(width=1)) +
    facet_grid(~ region) +
    theme_bw()  +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank(), legend.title=element_blank()) +
    xlab("") + 
    scale_x_discrete(breaks=NULL)
  return(plot)
}
