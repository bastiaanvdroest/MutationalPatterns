#' Estimate optimal rank for NMF decomposition
#' 
#' Find optimal number of signatures for NMF decomposition
#' @param mut_matrix 96 mutation count matrix
#' @param rank_range Range of ranks one would like to test 
#' @param nrun Number of runs to perform, default=100
#' @return NMF rank survey plot
#'
#' @examples
#' vcf_files = list.files(system.file("extdata", package="MutationalPatterns"),
#'                                    pattern = ".vcf",
#'                                    full.names = TRUE)
#' sample_names = c("colon1", "colon2", "colon3",
#'                  "intestine1", "intestine2", "intestine3",
#'                  "liver1", "liver2", "liver3")
#' vcfs = read_vcf(vcf_files, sample_names, genome = "hg19")
#'
#' # only select autosomal chromosomes, mt dna length is different for vcf and
#' # ref genome, why??
#' auto = extractSeqlevelsByGroup(species="Homo_sapiens",
#'                                style="UCSC",
#'                                group="auto")
#' vcfs = lapply(vcfs, function(x) keepSeqlevels(x, auto))
#'
#' ref_genome = "BSgenome.Hsapiens.UCSC.hg19" 
#' library(ref_genome, character.only = TRUE)
#' tissue = c("colon", "colon", "colon",
#'            "intestine", "intestine", "intestine",
#'            "liver", "liver", "liver")
#' vcfs = lapply(vcfs, function(x) rename_chrom(x))
#'
#' test_matrix = mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
#' estimate_rank(test_matrix, rank_range = 2:5, nrun = 50)
#'
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
