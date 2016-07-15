#' Calculate strand bias
#'
#' A function to calculate strand bias.
#' @param mut_mat_s The mutation matrix
#' @param by The column to group by
#' @return A data frame with strand bias data
#' @importFrom reshape2 melt
#' @export

strand_bias = function(mut_mat_s, by)
{
  df = t(mut_mat_s)
  # check if grouping paramter by was provided
  if(missing(by)){by = rep("all", nrow(df))}
  x = aggregate(df, by=list(by), FUN=sum) 
  rownames(x) = x[,1]
  x = x[,-1]
  x_r = x/ rowSums(x)
  substitutions = rep(SUBSTITUTIONS,each=32)
  x2 = melt(aggregate(t(x), by = list(substitutions, STRAND), FUN=sum))
  x2_r = melt(aggregate(t(x_r), by = list(substitutions, STRAND), FUN=sum))
  colnames(x2) = c("type", "strand", "group", "no_mutations")
  colnames(x2_r) = c("type", "strand", "group", "relative_contribution")
  res = merge(x2, x2_r)
  return(res)
}
