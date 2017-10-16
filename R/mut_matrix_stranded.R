#' Make mutation count matrix of 96 trinucleotides with 
#' strand information
#'
#' Make a mutation count matrix with 192 features: 96 trinucleotides and 2 strands,
#' these can be transcriptional strand
#'
#' @param vcf_list List of collapsed vcf objects
#' @param ref_genome BSGenome reference genome object 
#' @param ranges GRanges object with the genomic ranges of the two strand
#' @param mode "transcription" or "replication", default = "transcription"
#'
#' @return 192 mutation count matrix (96 * 2 strands)
#'
#' @import GenomicRanges
#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#'
#' ## Load the corresponding reference genome.
#' ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' ## You can obtain the known genes from the UCSC hg19 dataset using
#' ## Bioconductor:
#' # source("https://bioconductor.org/biocLite.R")
#' # biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' # library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#'
#' ## For this example, we preloaded the data for you:
#' genes_hg19 <- readRDS(system.file("states/genes_hg19.rds",
#'                         package="MutationalPatterns"))
#'
#' mut_mat_s = mut_matrix_stranded(vcfs, ref_genome, genes_hg19)
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#' \code{link{mut_matrix}}
#'
#' @export

mut_matrix_stranded = function(vcf_list, ref_genome, ranges, mode = "transcription")
{
  df = data.frame()
  
  # Detect number of cores available for parallelization
  num_cores = detectCores()
  if (!(.Platform$OS.type == "windows" || is.na(num_cores)))
    num_cores <- detectCores()
  else
    num_cores = 1
  
  # Transcription mode
  if(mode == "transcription")
    {
      # For each vcf in vcf_list count the 192 features
      rows <- mclapply (as.list(vcf_list), function (vcf)
      {
        type_context = type_context(vcf, ref_genome)
        strand = strand_from_vcf(vcf, ranges, mode = "transcription")
        row = mut_192_occurrences(type_context, strand)
        return(row)
      }, mc.cores = num_cores, mc.silent = FALSE)
      
      # Combine the rows in one dataframe
      for (row in rows)
        df = rbind (df, row)
    }
  
  # Replication mode
  if(mode == "replication")
  {
    # For each vcf in vcf_list count the 192 features
    rows <- mclapply (as.list(vcf_list), function (vcf)
    {
      type_context = type_context(vcf, ref_genome)
      strand = strand_from_vcf(vcf, ranges, mode = "replication")
      row = mut_192_occurrences(type_context, strand)
      return(row)
    }, mc.cores = num_cores, mc.silent = FALSE)
    
    # Combine the rows in one dataframe
    for (row in rows)
      df = rbind (df, row)
  }
  
  # Add row names to data.frame
  names(df) = names(row)
  row.names(df) = names(vcf_list)
  
  # Transpose and return
  return(t(df))
}
