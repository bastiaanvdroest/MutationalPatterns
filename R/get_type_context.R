#' Retrieve context of base substitution types
#' 
#' A function to extract the bases 3' upstream and 5' downstream of the base substitution types
#' @param vcf A CollapsedVCF object
#' @param ref_genome Reference genome
#' @return Mutation types and context character vectors in a named list
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

    # change the context of these mutations to reverse complement of the context
    y = reverse(chartr('ATGC', 'TACG', y))

    # replace subset with reverse complement
    mut_context[x] = y

    # return as named list
    res = list(types, mut_context)
    names(res) = c("types", "context")

    return(res)
}
