#' Find transcriptional strand of base substitutions in vcf
#'  
#' @description For the positions that are within gene bodies it is determined whether the "C" or "T" base 
#' (since by convention we regard base substitutions as C>X or T>X) are on the same strand as the gene definition.
#' Base substitions on the same strand as the gene definitions are considered "untranscribed", and on the opposite strand of gene bodies as transcribed, 
#' since the gene definitions report the coding or sense strand, which is untranscribed.  
#' No strand information "-" is reported for base substitution that overlap with more than one gene body.
#' @param vcf Granges vcf object
#' @param genes Granges with definition of gene bodies, should include strand information
#' @return Character vector with transcriptional strand information with length of vcf: "-" for positions outside gene bodies, "U" for untranscribed/sense/coding strand, "T" for transcribed/anti-sense/non-coding strand
#' @export

get_strand = function(vcf, genes)
{
  # check chromosome names
  if (!(all( seqlevels(vcf) %in% seqlevels(genes)) ))
    stop("Chromosome names (seqlevels) of vcf and genes Granges object do not match. Use rename_chrom() function to rename chromosome names.")

  # determine overlap between vcf positions and genes
  overlap = findOverlaps(vcf, genes)
  overlap = as.data.frame(as.matrix(overlap))
  colnames(overlap) = c('vcf_id', 'gene_body_id')
  # remove mutations that overlap with multiple genes and therefore cannot be 
  # determined whether they are on transcribed or untranscribed strand
  # duplicated mutations
  dup_pos = overlap$vcf_id[duplicated(overlap$vcf_id)]
  # index of duplicated mutations
  dup_idx = which(overlap$vcf_id %in% dup_pos)
  # remove all duplicated (non-unique mapping) mutations
  if(length(dup_idx) > 0){overlap = overlap[-dup_idx,]}
  # subset of mutations in genes
  vcf_overlap = vcf[overlap$vcf_id]
  # find reference allele of mutations (+ strand of reference genome is reported in vcf file)
  ref = vcf_overlap$REF
  # Find the strand of C or T (since we regard base substitutions as C>X or T>X)
  # which mutations have ref allele C or T
  i = which(ref == "C" | ref == "T")
  # store mutation strand info in vector
  strand_muts = rep(0, nrow(overlap))
  strand_muts[i] = "+"
  strand_muts[-i] = "-"
  # find strand of gene bodies of overlaps
  strand_genebodies = as.character(strand(genes)[overlap$gene_body_id])
  # find if mut and gene_bodies are on the same strand
  same_strand = (strand_muts  == strand_genebodies)
  # subset vcf object for both untranscribed and transcribed
  # gene definition represents the untranscribed/sense/coding strand
  # if mutation is on same strand as gene, than its untranscribed
  U_index = which(same_strand == TRUE)
  # if mutation is on different strand than gene, than its transcribed
  T_index = which(same_strand == FALSE)
  strand = rep(0, nrow(overlap))
  strand[U_index] = "U"
  strand[T_index] = "T"
  # make vector with all positions in input vcf
  # for positions that do not overlap with gene bodies, report "-"
  strand2 = rep("-", length(vcf))
  strand2[overlap$vcf_id] = strand
  return(strand2)
}
