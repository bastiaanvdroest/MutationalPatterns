#' Count the occurences of each base substitution type
#' 
#' @param vcf_list A list of CollapsedVCF object
#' @param ref_genome Reference genome
#' @return Dataframe with counts of each base substitution type for each sample in vcf_list
#' @export
#' @import BiocGenerics

mut_type_occurences = function(vcf_list, ref_genome)
{  
    n_samples = length(vcf_list)
    df = data.frame()

    print("Counting base substitution type occurrences for sample:")

    for(i in 1:n_samples)
    {
        print(i)
        vcf = vcf_list[[i]]
        types = get_types(vcf)
        CT_muts = which(types == "C>T")
        CT_context = get_type_context(vcf[CT_muts], ref_genome)[[2]]
        CpG = c("ACG", "CCG", "TCG", "GCG")
        CT_at_CpG = sum(!(is.na(BiocGenerics::match(CT_context,CpG))))
        CT_at_other = length(CT_muts) - CT_at_CpG
        counts = as.vector(table(types))
        counts = c(counts, CT_at_CpG, CT_at_other)
        df = rbind(df,counts)
    }

    row.names(df) = names(vcf_list)
    colnames(df) = c(   "C>A", "C>G", "C>T",
                        "T>A", "T>C", "T>G",
                        "C>T at CpG", "C>T other")
    return(df)
}
