#' Retrieve base substitution types from vcf
#' 
#' A function to extract the base substitutions from a vcf and translate to
#' the 6 common base substitution types.
#' 
#' @param vcf A CollapsedVCF object
#' @return Character vector with base substitution types
#'
#' @examples
#' ## See the 'vcf_to_granges()' example for how we obtained the following data:
#' vcfs <- readRDS(system.file("states/vcf_to_granges_output.R",
#'                 package="MutationalPatterns"))
#'
#' get_types(vcfs[[1]])
#'
#' @seealso
#' \code{\link{vcf_to_granges}}
#'
#' @export

get_types = function(vcf) 
{
    muts = get_muts(vcf)
    types = unlist(muts)
    types = gsub('G>T', 'C>A', types)
    types = gsub('G>C', 'C>G', types)
    types = gsub('G>A', 'C>T', types)
    types = gsub('A>T', 'T>A', types)
    types = gsub('A>G', 'T>C', types)
    types = gsub('A>C', 'T>G', types)
    return(types)
}
