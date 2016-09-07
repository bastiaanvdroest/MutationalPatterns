#' Make mutation count matrix of 96 trinucleotides with transcriptional strand information
#'  
#' @description Make mutation count matrix for 96 trinucleotides for both transcribed and untranscribed strand of gene bodies. 
#' Mutations outside gene bodies are not counted.
#' @param vcf_list List of collapsed vcf objects
#' @param ref_genome BSGenome reference genome object 
#' @param genes Granges with definition of gene bodies, should include strand information
#' @return 192 mutation count matrix (96 * 2 strands)
#' @import GenomicRanges
#' @export

mut_matrix_stranded = function(vcf_list, ref_genome, genes)
{
    df = data.frame()
    for(vcf in vcf_list)
    {
        type_context = get_type_context(vcf, ref_genome)
        strand = get_strand(vcf, genes)
        row = mut_192_occurences(type_context, strand)
        df = rbind(df, row)
    }

    names(df) = names(row)
    row.names(df) = names(vcf_list)

    # transpose
    return(t(df))
}
