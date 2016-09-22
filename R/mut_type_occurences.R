#' Count the occurences of each base substitution type
#' 
#' @param vcf_list A list of CollapsedVCF object
#' @param ref_genome Reference genome
#' @return data.frame with counts of each base substitution type for
#' each sample in vcf_list
#'
#' @importFrom BiocGenerics rbind
#'
#' @examples
#' ## See the 'vcf_to_granges()' example for how we obtained the following data:
#' vcfs <- readRDS(system.file("states/vcf_to_granges_output.rds",
#'                     package="MutationalPatterns"))
#' 
#' ## Exclude mitochondrial and allosomal chromosomes.
#' autosomal <- extractSeqlevelsByGroup(species="Homo_sapiens",
#'                                         style="UCSC",
#'                                         group="auto")
#'
#' vcfs <- lapply(vcfs, function(x) keepSeqlevels(x, autosomal))
#'
#' ## Load a reference genome.
#' ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' ## Get the type occurrences for all VCF objects.
#' type_occurences = mut_type_occurences(vcfs, ref_genome)
#'
#' @seealso
#' \code{\link{vcf_to_granges}},
#'
#' @export

mut_type_occurences = function(vcf_list, ref_genome)
{  
    n_samples = length(vcf_list)
    df = data.frame()

    for(i in 1:n_samples)
    {
        vcf = vcf_list[[i]]
        types = get_types(vcf)
        CT_muts = which(types == "C>T")
        CT_context = get_type_context(vcf[CT_muts], ref_genome)[[2]]
        CpG = c("ACG", "CCG", "TCG", "GCG")
        CT_at_CpG = sum(!(is.na(BiocGenerics::match(CT_context,CpG))))
        CT_at_other = length(CT_muts) - CT_at_CpG

        # Construct a table and handle missing mutation types.
        column_names = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G",
                            "C>T at CpG", "C>T other")

        types_table = table(types)
        full_table = table(column_names)
        for(i in 1:length(column_names) - 2)
        {
            value = as.vector(types_table[column_names[i]])[1]

            # Set NA to 0
            if (is.na(value))
                value = 0

            full_table[column_names[i]] = value
        }

        full_table["C>T at CpG"] = CT_at_CpG
        full_table["C>T other"] = CT_at_other

        df = rbind(df,full_table)
    }

    row.names(df) = names(vcf_list)
    colnames(df) = names(full_table)
    return(df)
}
