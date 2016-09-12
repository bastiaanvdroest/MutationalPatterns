#' Test for enrichment or depletion of mutations in genomic regions
#'
#' This function aggregates mutations per group (optional) and performs an
#' enrichment depletion test.
#' @param x Data.frame result from genomic_distribution() 
#' @param by Optional grouping variable, e.g. tissue type
#' @return Data.frame with the observed and expected number of mutations per
#' genomic region per group (by) or sample
#' @importFrom BiocGenerics cbind
#' @importFrom BiocGenerics rbind
#'
#' @examples
#' # See the 'read_vcf()' example for how we obtained the following data:
#' vcfs <- readRDS(system.file("states/read_vcf_output.R",
#'                 package="MutationalPatterns"))
#' 
#' # Rename the seqlevels to the UCSC standard.
#' vcfs <- lapply(vcfs, rename_chrom)
#'
#' # Exclude mitochondrial and allosomal chromosomes.
#' autosomal = extractSeqlevelsByGroup(species="Homo_sapiens",
#'                                     style="UCSC",
#'                                     group="auto")
#'
#' vcfs = lapply(vcfs, function(x) keepSeqlevels(x, autosomal))
#'
#' # We need to retrieve data using biomaRt to do a sensible genomic
#' # distribution analysis.
#' library(biomaRt)
#' #segmentation = useEnsembl(biomart="regulation",
#' #                          dataset="hsapiens_segmentation_feature",
#' #                          GRCh = 37)
#' regulatory = useEnsembl(biomart="regulation",
#'                         dataset="hsapiens_regulatory_feature",
#'                         GRCh = 37)
#' #annotated = useEnsembl(biomart="regulation",
#' #                       dataset="hsapiens_annotated_feature",
#' #                       GRCh = 37)
#' CTCF = getBM(attributes = c('chromosome_name',
#'                             'chromosome_start',
#'                             'chromosome_end',
#'                             'feature_type_name',
#'                             'cell_type_name'),
#'              filters = "regulatory_feature_type_name", 
#'              values = "CTCF Binding Site", 
#'              mart = regulatory)
#'
#' # Make a GRanges object.
#' CTCF_g = reduce(GRanges(CTCF$chromosome_name,
#'                 IRanges(CTCF$chromosome_start,
#'                 CTCF$chromosome_end)))
#' regions = list(CTCF_g)
#' names(regions) = c("CTCF")
#' regions = lapply(regions, function(x) rename_chrom(x))
#'
#' # Get the filename with surveyed/callable regions
#' surveyed_file = list.files(system.file("extdata",
#'                                        package="MutationalPatterns"),
#'                            pattern = ".bed",
#'                            full.names = TRUE)
#' # Read the file as a GRanges object.
#' surveyed_list = bed_to_granges(surveyed_file, "surveyed_all")
#'
#' # For this example we use the same surveyed file for each sample.
#' surveyed_list= rep(surveyed_list, 9)
#' 
#' # Calculate the number of observed and expected number of mutations in
#' # each genomic regions for each sample.
#' distr = genomic_distribution(vcfs, surveyed_list, regions)
#'
#' tissue = c("colon", "colon", "colon",
#'            "intestine", "intestine", "intestine",
#'            "liver", "liver", "liver")
#'
#' distr_test = enrichment_depletion_test(distr, by = tissue)
#' distr_test2 = enrichment_depletion_test(distr)
#'
#' @seealso \code{\link{genomic_distribution}}
#'
#' @export

enrichment_depletion_test = function(x, by = c())
{
    # Handle the 'by' parameter when necessary by aggregating x
    if (length(by) > 0){
        x$by = by
        # sum the columns while aggregating rows based on unique values in 'by'
        # and 'region'.
        res2 = stats::aggregate(cbind(n_muts,
                                        surveyed_length,
                                        surveyed_region_length,
                                        observed) ~ by + region,
                                data = x, sum)
    }
    else {
        res2 = x
        # by variable is sample variable
        res2$by = res2$sample
        # select output columns
        res2 = res2[,c(9,1,3,4,6,8)]
    }

    # Calculate probability and expected number of mutations
    res2$prob = res2$n_muts / res2$surveyed_length
    res2$expected = res2$prob * res2$surveyed_region_length

    # Perform enrichment/depletion test for each row
    res3 = data.frame()
    for(i in 1:nrow(res2))
    {
        x = res2[i,]
        res3 = rbind(res3, binomial_test(x$prob, x$surveyed_region_length,  x$observed))
    }

    # Combine results into one data frame
    df = cbind(res2, res3)
    return(df)
}
