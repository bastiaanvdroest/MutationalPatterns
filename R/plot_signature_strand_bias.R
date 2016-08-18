#' Plot signature strand bias
#' 
#' @description Plot strand bias per mutation type for each signature
#' @param signatures_strand_bias Signature matrix with 192 features
#' @return Barplot
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 position_dodge
#' @importFrom plyr adply
#' @export



plot_signature_strand_bias = function(signatures_strand_bias)
{
  # check if there are 192 features in the signatures
  if(dim(signatures_strand_bias)[1] != 192){stop("Input signature matrix does not have 192 features (96 trinucleotide * 2 strands).")}

  # aggregate by strand and type
  sum_per_type = aggregate(signatures_strand_bias, by=list(STRAND, SUBSTITUTIONS_192), FUN=sum)
  sum_per_strand = aggregate(signatures_strand_bias, by=list(STRAND), FUN=sum)
  # melt data frames
  sum_per_strand =  melt(sum_per_strand)
  colnames(sum_per_strand) = c("strand", "Signature", "value")
  sum_per_type =  melt(sum_per_type)
  colnames(sum_per_type) = c("strand", "type", "Signature", "value")
  # ratio per signature per type
  ratio = as.matrix(subset(sum_per_type,strand == "T")$value / subset(sum_per_type,strand == "U")$value)
  ratio_per_type_per_signature = cbind(subset(sum_per_type,strand == "T")[,2:3], ratio)
  
  
  # binomial test per type per signature
  size = c()
  observed = c()
  transcribed = c()
  untranscribed = c()
  for(s in unique(sum_per_type$Signature))
  {
    for(t in unique(sum_per_type$type))
    {
      sub = subset(sum_per_type, Signature==s & type==t)
      size = c(size,sum(sub$value))
      observed = c(observed,subset(sub, strand == "T")$value)
      transcribed = c(transcribed, subset(sub, strand == "T")$value)
      untranscribed = c(untranscribed, subset(sub, strand == "U")$value)
    }
  }
  
  stats_per_type = data.frame(Signature = c(rep("Signature A",6), rep("Signature B",6)), type = rep(SUBSTITUTIONS,2), size = as.integer(size), transcribed = transcribed, untranscribed = untranscribed, observed = as.integer(observed))
  stats_per_type = adply(stats_per_type, 1, function(x) binomial_test(0.5, x$size, x$observed))
  
  ratio_per_type_per_signature = cbind(ratio_per_type_per_signature, stats_per_type)
  strand_bias_per_type_df = melt(ratio_per_type_per_signature[,c(1,2,3,12)])
  
  max = round(max(abs(log2(strand_bias_per_type_df$value))))
  
  plot = ggplot(strand_bias_per_type_df, aes(x=type, y=log2(value), fill=type)) +
    geom_bar(stat="identity", position="dodge", color="black") +
    scale_y_continuous(limits=c(-max,max)) +
    scale_fill_manual(values=COLORS6) +
    facet_grid(Signature ~ .) +
    ylab("log2(transcribed/untranscribed)") +
    theme_bw() + 
    scale_x_discrete(breaks=NULL) +
    xlab("") +
    geom_text(aes(x = type, y = log2(value), ymax = log2(value), 
                  label=significant, vjust=ifelse(sign(log2(value)) > 0, 0.5, 1)), 
              size = 8, position = position_dodge(width=1)) 
  
  return(plot)
}
