#' Make mutation count matrix of 96 trinucleotides 
#'  
#' @description Make 96 trinucleotide mutation count matrix
#' @param vcf_list List of collapsed vcf objects
#' @param ref_genome BSGenome reference genome object 
#' @return 96 mutation count matrix
#' @import GenomicRanges
#' @export

mut_matrix = function(vcf_list, ref_genome)
{
    df = data.frame()
    for(vcf in vcf_list)
    {
        type_context = get_type_context(vcf, ref_genome)
        row = mut_96_occurences(type_context)
        df = rbind(df, row)
    }

    names(df) = names(row)
    row.names(df) = names(vcf_list)

    # transpose
    return(t(df))
}
