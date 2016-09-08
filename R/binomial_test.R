#' Binomial test for enrichment or depletion testing
#'
#' This function performs lower-tail binomial test for depletion and
#' upper-tail test for enrichment
#'
#' @param p Probability of success
#' @param n Number of trials
#' @param x Observed number of successes
#' @return A Data.frame with direction of effect (enrichment/depletion),
#' P-value and significance asterisks
#' @export

binomial_test = function(p, n, x)
{
    # Calculate expected number of successes
    expected = p * n

    # Handle depletion
    if (x < expected)
    {
        # do lower tail test
        pval = pbinom(x, n, p, lower.tail=TRUE)
        effect = "depletion"
    }

    # Handle enrichment
    else
    {
        # do upper tail test
        pval = pbinom(x-1, n, p, lower.tail=FALSE)
        effect = "enrichment"
    }

    # Add significance asteriks
    if (pval < 0.05)
        significant = "*"
    else
        significant = ""

    res = data.frame(effect, pval, significant)
    return(res)
}
