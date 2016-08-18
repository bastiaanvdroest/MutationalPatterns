#' Significance test for transcriptional strand asymmetry
#' 
#' @description Performs a Poisson test for the ratio between mutations on the transcribed and untranscribed strand
#' @param strand_occurences Dataframe with mutation count per strand, result from strand_occurences()
#' @return Dataframe with poisson test P value for the ratio between the transcribed and untrascribed strand per group per base substitution type
#' @importFrom reshape2 dcast
#' @importFrom plyr .
#' @export

strand_bias_test = function(strand_occurences)
{
  # statistical test for strand ratio
  # poisson test
  df_strand = dcast(melt(strand_occurences), type + group ~ strand, sum, subset = .(variable == "no_mutations"))
  df_strand = ddply(df_strand, c("group", "type", "T", "U"), summarise, total = T+U, ratio = T/U,  p_poisson = poisson.test(c(U,T), r=1)$p.value)
  df_strand$significant[df_strand$p_poisson < 0.05] = "*"
  df_strand$significant[df_strand$p_poisson >= 0.05] = " "
  return(df_strand)
}
