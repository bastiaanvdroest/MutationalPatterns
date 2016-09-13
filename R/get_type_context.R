#' Retrieve context of base substitution types
#' 
#' A function to extract the bases 3' upstream and 5' downstream of the base
#' substitution types.
#'
#' @param vcf A CollapsedVCF object
#' @param ref_genome Reference genome
#' @return Mutation types and context character vectors in a named list
#'
#' @examples
#' ## See the 'read_vcf()' example for how we obtained the following data:
#' vcfs <- readRDS(system.file("states/read_vcf_output.R",
#'                 package="MutationalPatterns"))
#'
#' ## Rename the seqlevels to the UCSC standard.
#' vcfs <- lapply(vcfs, rename_chrom)
#'
#' ## Exclude mitochondrial and allosomal chromosomes.
#' autosomal = extractSeqlevelsByGroup(species="Homo_sapiens",
#'                                     style="UCSC",
#'                                     group="auto")
#'
#' vcfs = lapply(vcfs, function(x) keepSeqlevels(x, autosomal))
#'
#' ## Load the corresponding reference genome.
#' ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' mut_context = get_type_context(vcfs[[1]], ref_genome)
#'
#' @seealso \code{\link{read_vcf}}, \code{\link{rename_chrom}},
#'          \code{\link{get_mut_context}}
#'
#' @export

get_type_context = function(vcf, ref_genome)
{
    mut_context = get_mut_context(vcf, ref_genome)
    muts = get_muts(vcf)
    types = get_types(vcf)

    # find the mutations for which the context needs to be adjusted
    x = which(muts != types)

    # subset mut_context
    y = mut_context[x]

    # Change the context of these mutations to reverse complement
    # of the context
    y = reverse(chartr('ATGC', 'TACG', y))

    # replace subset with reverse complement
    mut_context[x] = y

    # return as named list
    res = list(types, mut_context)
    names(res) = c("types", "context")

    return(res)
}
