#' Find overlaps between mutations and a genomic region
#' 
#' Function finds the number of mutations that reside in genomic region and takes surveyed area of genome into account
#' @param vcf_list A list with vcf Granges objects
#' @param surveyed_list A list with Granges of regions of the genome that have been surveyed (e.g. determined using GATK CallableLoci)
#' @param region_list List with GRanges objects containing locations of genomic regions
#' @return A Data.frame containing the number observed and number of expected mutations in each genomic region.
#'
#' @examples
#' library(biomaRt)
#' mart="ensemble"
#' regulation_segmentation = useEnsembl(biomart="regulation",
#'                                      dataset="hsapiens_segmentation_feature",
#'                                      GRCh = 37)
#' regulation_regulatory = useEnsembl(biomart="regulation",
#'                                    dataset="hsapiens_regulatory_feature",
#'                                    GRCh = 37)
#' regulation_annotated = useEnsembl(biomart="regulation",
#'                                   dataset="hsapiens_annotated_feature",
#'                                   GRCh = 37)
#' CTCF = getBM(attributes = c('chromosome_name',
#'                             'chromosome_start',
#'                             'chromosome_end',
#'                             'feature_type_name',
#'                             'cell_type_name'),
#'              filters = "regulatory_feature_type_name", 
#'              values = "CTCF Binding Site", 
#'              mart = regulation_regulatory)
#' # Make a GRanges object.
#' CTCF_g = reduce(GRanges(CTCF$chromosome_name,
#'                 IRanges(CTCF$chromosome_start,
#'                 CTCF$chromosome_end)))
#'
#' # Get the filename with surveyed/callable regions
#' surveyed_file = list.files(system.file("extdata",
#'                                        package="MutationalPatterns"),
#'                            pattern = ".bed",
#'                            full.names = T)
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
#' @seealso \code{\link{bed_to_granges}}
#'
#' @export

genomic_distribution = function(vcf_list, surveyed_list, region_list)
{
    if (length(vcf_list) != length(surveyed_list))
        stop("vcf_list and surveyed_list must have the same length")

    df = data.frame()
    for(j in 1:length(region_list) )
    {
        print(paste("Region:", names(region_list)[j]))

        for(i in 1:length(vcf_list) )
        {
            print(paste("Sample:", names(vcf_list)[i]))
            res = intersect_with_region(vcf_list[[i]], surveyed_list[[i]], region_list[[j]])
            res$region = names(region_list)[j]
            res$sample = names(vcf_list)[i]
            res = res[,c(7,8,1:6)]
            df = rbind(df, res)
        }
    }

    # region as factor
    # make sure level order is the same as in region_list input (important for
    # plotting later)
    df$region = factor(df$region, levels = names(region_list))
    return(df)
}
