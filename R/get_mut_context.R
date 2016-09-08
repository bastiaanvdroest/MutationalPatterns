#' Retrieve context of base substitutions
#'
#' A function to extract the bases 3' upstream and 5' downstream of the base
#' substitutions from the reference genome
#' @param vcf A Granges object
#' @param ref_genome Reference genome
#' @return Character vector with the context of the base substitutions
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom Biostrings getSeq
#' @export

get_mut_context = function(vcf, ref_genome) 
{
    # check if chromosome names are the same
    if (!(all( seqlevels(vcf) %in% seqlevels(get(ref_genome)) )))
        stop("Chromosome names (seqlevels) of vcf and reference genome object do not match. Use rename_chrom() function to rename chromosome names.")

    ranges = resize(vcf, 3, fix = "center")
    vcf_context = getSeq(get(ref_genome), ranges)
    return(vcf_context)
}
