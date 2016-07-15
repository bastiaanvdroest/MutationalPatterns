#' Find overlap between mutations and a genomic region
#' 
#' Find the number of mutations that reside in genomic region and take surveyed area of genome into account
#' @param vcf CollapsedVCF object with mutations
#' @param surveyed GRanges object with regions of the genome that were surveyed
#' @param region GRanges object with genomic region(s)
#' @importFrom GenomeInfoDb seqlevelsStyle
#' 
intersect_with_region = function(vcf, surveyed, region)
{
  # number of mutations in vcf file
  n_muts = length(vcf)
  # number of bp that was surveyed
  surveyed_length = sum(as.numeric(width(surveyed)))
  
  # check if chromosome names are the same in the objects
  if (seqlevelsStyle(vcf) != seqlevelsStyle(surveyed))
    stop("Chromosome names (seqlevels) of vcf and surveyed granges object do not match.")

  if (seqlevelsStyle(region) != seqlevelsStyle(surveyed))
    stop("Chromosome names (seqlevels) of surveyed and region granges object do not match.")

  # Intersect genomic region and surveyed region
  surveyed_region = intersect(surveyed, region, ignore.strand = T)
  surveyed_region_length = sum(width(surveyed_region))
  

  # Find which mutations lie in surveyed genomic region
  overlap = findOverlaps(vcf, surveyed_region)
  muts_in_region = as.data.frame(as.matrix(overlap))$queryHits
  
  observed = length(muts_in_region)
  prob = n_muts / surveyed_length
  expected = prob * surveyed_region_length

  # output
  res = data.frame(n_muts, surveyed_length, prob, surveyed_region_length, expected, observed)
  return(res)
}
