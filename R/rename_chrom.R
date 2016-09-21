#' Rename chromosome names to standard
#' 
#' @param granges Granges object with chromosome names to rename
#' @param style "UCSC", NCBI" or "Ensembl", default = "UCSC"
#' @return A GRanges object with renamed chromosomes according to style.
#' @importFrom GenomeInfoDb seqlevels<-
#' @importFrom GenomeInfoDb mapSeqlevels
#' @importFrom GenomeInfoDb renameSeqlevels
#'
#' @examples
#' ## After loading VCF files using 'vcf_to_granges()', the chromosome names
#' ## may differ between the reference genome and your samples.  See the
#' ## 'vcf_to_granges()' example for how we obtain the following data.
#' vcfs <- readRDS(system.file("states/vcf_to_granges_output.R",
#'                 package="MutationalPatterns"))
#'
#' ## You can standardize the naming of your VCF samples loaded into R.
#' ## Notice how the seqnames change from:
#' vcfs$liver1
#'
#' ## To the UCSC (default) standardized format:
#' rename_chrom(vcfs$liver1)
#'
#' ## You can also rename all objects returned by 'vcf_to_granges()' using:
#' vcfs <- lapply(vcfs, rename_chrom)
#'
#' @seealso
#' \code{\link{vcf_to_granges}}
#'
#' @export

rename_chrom = function(granges, style = "UCSC")
{
    # rename mitochondrial DNA manually
    seqlevels(granges)[seqlevels(granges)=="chrMT"] = "chrM"

    # get chromosome style
    chrom_style = mapSeqlevels(seqlevels(granges), style)

    # removing NA cases
    chrom_style = chrom_style[complete.cases(chrom_style)] 

    # rename chromosome names (seqlevels)
    res = renameSeqlevels(granges, chrom_style)

    return(res)
}
