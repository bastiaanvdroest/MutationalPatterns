#' Test for enrichment or depletion of mutations in genomic regions
#'
#' This function aggregates mutations per group (optional) and performs an
#' enrichment depletion test.
#' @param x Dataframe result from genomic_distribution() 
#' @param by Optional grouping variable, e.g. tissue type
#' @return Data.frame with the observed and expected number of mutations per
#' genomic region per group (by) or sample
#' @importFrom BiocGenerics cbind
#' @importFrom BiocGenerics rbind
#' @export

enrichment_depletion_test = function(x, by = c())
{
    # Handle the 'by' parameter when necessary by aggregating x
    if (length(by) > 0){
        x$by = by
        # sum the columns while aggregating rows based on unique values in 'by'
        # and 'region'.
        res2 = stats::aggregate(cbind(n_muts,
                                        surveyed_length,
                                        surveyed_region_length,
                                        observed) ~ by + region,
                                data = x, sum)
    }
    else {
        res2 = x
        # by variable is sample variable
        res2$by = res2$sample
        # select output columns
        res2 = res2[,c(9,1,3,4,6,8)]
    }

    # Calculate probability and expected number of mutations
    res2$prob = res2$n_muts / res2$surveyed_length
    res2$expected = res2$prob * res2$surveyed_region_length

    # Perform enrichment/depletion test for each row
    res3 = data.frame()
    for(i in 1:nrow(res2))
    {
        x = res2[i,]
        res3 = rbind(res3, binomial_test(x$prob, x$surveyed_region_length,  x$observed))
    }

    # Combine results into one data frame
    df = cbind(res2, res3)
    return(df)
}
