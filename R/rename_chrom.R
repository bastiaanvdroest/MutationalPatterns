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
#' ## After loading VCF files using 'read_vcf()', the chromosome names
#' ## may differ between the reference genome and your samples.  See the
#' ## 'read_vcf()' example for how we obtain the following data.
#' vcfs <- readRDS(system.file("states/read_vcf_output.R",
#'                 package="MutationalPatterns"))
#'
#' ## You can standardize the naming of your VCF samples loaded into R.
#' ## Notice how the seqnames change from:
#' vcfs$liver1
#'
#' ## To the UCSC (default) standardized format:
#' rename_chrom(vcfs$liver1)
#'
#' ## You can also rename all objects returned by 'read_vcf()' using:
#' vcfs <- lapply(vcfs, rename_chrom)
#'
#' @seealso
#' \code{\link{read_vcf}}
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
