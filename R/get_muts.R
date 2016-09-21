#' Retrieve base substitutions from vcf
#' 
#' A function to extract base substitutions of each position in vcf
#' @param vcf A CollapsedVCF object
#' @return Character vector with base substitutions
#' @import GenomicRanges
#'
#' @examples
#' ## See the 'read_vcf()' example for how we obtained the following data:
#' vcfs <- readRDS(system.file("states/read_vcf_output.rds",
#'                 package="MutationalPatterns"))
#'
#' muts = get_muts(vcfs[[1]])
#'
#' @seealso
#' \code{\link{read_vcf}}
#'
#' @export

get_muts = function(vcf) 
{
    ref = as.character(vcf$REF)
    alt = as.character(unlist(vcf$ALT))
    muts = paste(ref, alt, sep=">")
    return(muts)
}

