#' Count occurences per base substitution type and transcriptional strand
#' 
#' @description For each base substitution type and transcriptional strand the total number of mutations
#' and the relative contribution within a group is returned
#' @param mut_mat_s 192 feature mutation count matrix, result from mut_matrix_stranded()
#' @param by Character vector with grouping info, optional
#' @return Dataframe with the total number of mutations and relative contribution within group per base substitution 
#' type and transcriptional strand (T = transcribed strand, U = untranscribed strand)
#' @importFrom reshape2 melt
#' @export

strand_occurences = function(mut_mat_s, by)
{
  df = t(mut_mat_s)
  # check if grouping parameter by was provided, if not group by all
  if(missing(by)){by = rep("all", nrow(df))}
  # sum by group
  x = aggregate(df, by=list(by), FUN=sum) 
  # add group as rownames
  rownames(x) = x[,1]
  x = x[,-1]
  # calculate relative contribution within group
  x_r = x/ rowSums(x)
  # sum per substition per strand
  substitutions = rep(SUBSTITUTIONS, each=32)
  x2 = melt(aggregate(t(x), by = list(substitutions, STRAND), FUN=sum))
  x2_r = melt(aggregate(t(x_r), by = list(substitutions, STRAND), FUN=sum))
  colnames(x2) = c("type", "strand", "group", "no_mutations")
  colnames(x2_r) = c("type", "strand", "group", "relative_contribution")
  # combine relative and absolute
  y = merge(x2, x2_r)
  # reorder group, type, strand
  y = y[,c(3,1,2,4,5)]
  y = y[order(y$group, y$type),]
  return(y)
}
