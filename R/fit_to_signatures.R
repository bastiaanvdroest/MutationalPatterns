#' Find optimal nonnegative linear combination of mutation signatures to reconstruct the mutation matrix
#' 
#' @description Find linear combination of mutation signatures that most closely reconstructs the mutation matrix by solving nonnegative least-squares constraints problem
#' @param mut_matrix 96 mutation count matrix (dimensions: 96 mutations X n samples)
#' @param signatures Signature matrix (dimensions: 96 mutations X n signatures)
#' @return Named list with signature contributions and reconstructed mutation matrix
#' @importFrom pracma lsqnonneg
#' @export

fit_to_signatures = function(mut_matrix, signatures)
{
    # make sure dimensions of input matrix are correct
    if (dim(mut_matrix)[1] != 96)
        stop("Mutation count matrix input should have dimensions 96 X n samples")

    if (dim(signatures)[1] != 96)
        stop("Signatures input should have dimensions 96 X n signatures")

    n_samples = dim(mut_matrix)[2]
    n_signatures = dim(signatures)[2]
    lsq_contribution = matrix(NA, nrow=n_signatures, ncol=n_samples)
    lsq_reconstructed = matrix(NA, nrow=96, ncol=n_samples)

    # Process each sample
    for (i in 1:ncol(mut_matrix))
    {
        y = mut_matrix[,i]
        lsq = lsqnonneg(signatures, y)
        lsq_contribution[,i] = lsq$x
        lsq_reconstructed[,i] = signatures %*% as.matrix(lsq$x) 
    }

    # Add row and col names
    sample_names = colnames(mut_matrix)
    signature_names = colnames(signatures)
    mut_type_names = rownames(signatures)

    colnames(lsq_contribution) = sample_names
    rownames(lsq_contribution) = signature_names

    colnames(lsq_reconstructed) = sample_names
    rownames(lsq_reconstructed) = mut_type_names

    res = list(lsq_contribution, lsq_reconstructed)
    names(res) = c("contribution", "reconstructed")

    return(res)
}
