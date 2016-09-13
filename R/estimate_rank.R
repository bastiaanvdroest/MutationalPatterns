#' Estimate optimal rank for NMF decomposition
#' 
#' Find optimal number of signatures for NMF decomposition
#' @param mut_matrix 96 mutation count matrix
#' @param rank_range Range of ranks one would like to test 
#' @param nrun Number of runs to perform, default=100
#' @return NMF rank survey plot
#'
#' @examples
#' ## See the 'read_vcf()' example for how we obtained the following data:
#' vcfs <- readRDS(system.file("states/read_vcf_output.R",
#'                 package="MutationalPatterns"))
#' 
#' # Rename the seqlevels to the UCSC standard.
#' vcfs <- lapply(vcfs, rename_chrom)
#'
#' # Exclude mitochondrial and allosomal chromosomes.
#' autosomal = extractSeqlevelsByGroup(species="Homo_sapiens",
#'                                     style="UCSC",
#'                                     group="auto")
#'
#' vcfs = lapply(vcfs, function(x) keepSeqlevels(x, autosomal))
#'
#' ## Define the reference genome we are going to use.
#' ref_genome = "BSgenome.Hsapiens.UCSC.hg19" 
#'
#' ## If neccessary, download it from Bioconductor.
#' # source("http://bioconductor.org/biocLite.R")
#' # biocLite(ref_genome)
#'
#' ## Load the reference genome.
#' library(ref_genome, character.only = TRUE)
#'
#' test_matrix = mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
#' estimate_rank(test_matrix, rank_range = 2:5, nrun = 50)
#'
#' @seealso
#' \code{\link{read_vcf}},
#' \code{\link{rename_chrom}}
#'
#' @export

estimate_rank = function(mut_matrix, rank_range, nrun=100)
{
    mut_matrix = as.matrix(mut_matrix)

    # Add small pseudocount
    mut_matrix = mut_matrix + 0.0001

    # Check if rank_range is appropriate
    if (ncol(mut_matrix) < max(rank_range))
        stop("The maximum rank should be smaller than the number of samples")

    # Estimate ranks
    print("Estimating ranks...")
    estim.r = nmf(mut_matrix,
                  rank = rank_range,
                  method = "brunet",
                  nrun = nrun,
                  seed = 123456)

    # Plot result
    plot = plot(estim.r)
    return(plot)
}
