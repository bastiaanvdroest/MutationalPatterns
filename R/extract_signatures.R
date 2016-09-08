#' Extract mutational signatures from 96 mutation matrix using NMF
#' 
#' Decomposes trinucleotide count matrix into signatures and contribution of those signatures to the spectra of the samples/vcf files
#' @param mut_matrix 96 mutation count matrix 
#' @param rank Number of signatures to extract
#' @param nrun Number of iterations, default = 200
#' @return Named list of mutation matrix, signatures and signature contribution
#' @importFrom NMF nmf
#' @importFrom NMF basis
#' @importFrom NMF coef
#' @export

extract_signatures = function(mut_matrix, rank, nrun = 200)
{
    mut_matrix = as.matrix(mut_matrix)

    # add small pseudocount to avoid features with zero counts
    mut_matrix = mut_matrix + 0.0001

    # Check if rank_range is appropriate
    if (!(rank > 0 & rank == round(rank)))
        stop("Rank should be a positive integer")

    if (ncol(mut_matrix) < max(rank))
        stop("Rank should be smaller than the number of samples in the input matrix")

    # Calculate NMF
    print("Decomposing matrix using NMF...")
    res = nmf(mut_matrix,
                rank = rank,
                method = "brunet",
                nrun=nrun,
                seed = 123456)

    print(paste("Number of iterations:", nrun))

    # Find signatures and contribution of signatures
    signatures = basis(res)
    contribution = coef(res)

    # Reconstruct mutation matrix
    reconstructed = signatures %*% contribution
    return(list(signatures = signatures,
                contribution = contribution,
                reconstructed = reconstructed))
}
