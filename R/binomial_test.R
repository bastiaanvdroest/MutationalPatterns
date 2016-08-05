#' Binomial test for enrichment or depletion testing
#' 
#' Performs lower.tail binomial test for depletion and upper tail test for enrichment
#' 
#' @param p Probability of success
#' @param n Number of trials
#' @param x Observed number of successes
#' @return data.frame With direction of effect (enrichment/depletion), P value and significance asterisks
#' @export


binomial_test = function(p, n, x)
{
  # calculate expected number of successes
  expected = p * n
  # if observed is less than expected
  if(x < expected)
  {
    # For depletion
    # do lower tail test
    pval = pbinom(x, n, p, lower.tail=TRUE)
    effect = "depletion"
  }else{
    # For enrichment
    # do upper tail test
    pval = pbinom(x-1, n, p, lower.tail=FALSE)
    effect = "enrichment"
  }
  # add significance asteriks
  if(pval < 0.05){
    significant="*"
  }else{significant=""} 
  
  res = data.frame(effect, pval, significant)
  return(res)
}
