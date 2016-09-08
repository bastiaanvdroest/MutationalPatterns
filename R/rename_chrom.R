#' Rename chromosome names to standard
#' 
#' @param granges Granges object with chromosome names to rename
#' @param style "UCSC", NCBI" or "Ensembl", default = "UCSC"
#' @return A GRanges object with renamed chromosomes according to style.
#' @importFrom GenomeInfoDb seqlevels<-
#' @importFrom GenomeInfoDb mapSeqlevels
#' @importFrom GenomeInfoDb renameSeqlevels
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
