#' Estimate optimal rank for NMF decomposition
#' 
#' Find optimal number of signatures for NMF decomposition
#' @param mut_matrix 96 mutation count matrix
#' @param rank_range Range of ranks one would like to test 
#' @param nrun Number of runs to perform, default=100
#' @return NMF rank survey plot
#' @export

estimate_rank = function(mut_matrix, rank_range, nrun=100)
{
    mut_matrix = as.matrix(mut_matrix)

    # Add small pseudocount
    mut_matrix = mut_matrix + 0.0001

    # Check if rank_range is appropriate
    if (ncol(mut_matrix) < max(rank_range))
        stop("Maximum rank should be smaller than the number of samples")

    # Estimate ranks
    print("Estimating ranks...")
    estim.r = nmf(mut_matrix, rank = rank_range, method = "brunet", nrun = nrun, seed = 123456)

    # Plot result
    plot = plot(estim.r)
    return(plot)
}
