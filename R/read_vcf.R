#' Read vcf files into list of CollapsedVCF objects
#' 
#' Function reads Variant Call Format VCF files into a GRanges object and combines them in a list object
#' @param vcf_files Character vector of vcf file names
#' @param sample_names Character vector of sample names
#' @param genome A character or Seqinfo object
#' @return List of GRanges objects
#' @importFrom BiocGenerics match
#' @importFrom VariantAnnotation readVcf
#' @importFrom SummarizedExperiment rowRanges
#'
#' @examples
#' # The example data set consists of three colon samples, three intestine
#' # samples and three liver samples.  So, to map each file to its appropriate
#' # sample name, we create a vector containing the sample names:
#' sample_names <- c("colon1", "colon2", "colon3",
#'                     "intestine1", "intestine2", "intestine3",
#'                     "liver1", "liver2", "liver3")
#'
#' # We assemble a list of files we want to load.  These files match the sample
#' # names defined above.
#' vcf_files <- list.files(system.file("extdata", package="MutationalPatterns"),
#'                                     pattern = ".vcf", full.names = TRUE)
#'
#' # This function loads the files as CollapsedVCF objects
#' vcfs <- read_vcf(vcf_files, sample_names, genome = "hg19")
#'
#' @export

read_vcf = function(vcf_files, sample_names, genome = "-")
{
    # check sample names
    if (length(vcf_files) != length(sample_names))
        stop("Provide the same number of sample names as vcf files")

    # make list with vcf objects
    vcf_list = list()
    for(i in 1:length(vcf_files))
    {
        # parse vcf file with variantannotation
        vcf = readVcf(vcf_files[i], genome)

        # only store genomic ranges info (to reduce memory usage)
        vcf = rowRanges(vcf)

        # find and exclude positions with indels or multiple alternative alleles
        rem = which(all(!( !is.na(match(vcf$ALT, DNA_BASES)) &
                           !is.na(match(vcf$REF, DNA_BASES)) &
                           (lengths(vcf$ALT) == 1) )))
        if (length(rem) > 0)
        {
            vcf = vcf[-rem]
            warning(length(rem), " position(s) with indels and multiple alternative alleles are removed.")
        }

        vcf = list(vcf)
        names(vcf) = sample_names[i]
        vcf_list = c(vcf_list, vcf)
    }

    return(vcf_list)
}
