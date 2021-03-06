% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mut_matrix.R
\name{mut_matrix}
\alias{mut_matrix}
\title{Make mutation count matrices of single and double substitutions and indels}
\usage{
mut_matrix(vcf_list, ref_genome, type, method = "split", num_cores)
}
\arguments{
\item{vcf_list}{List of collapsed vcf objects}

\item{ref_genome}{BSGenome reference genome object}

\item{type}{(Optional) A character vector stating which type of mutation is to be extracted: 
'snv', 'dbs' and/or 'indel'. All mutation types can also be chosen by 'type = all'.\cr
Default is 'snv'}

\item{method}{(Optional) Character stating how to use the data. 
\itemize{
  \item{"split":} { Each mutation type has seperate count matrix}
  \item{"combine":} { Combined count matrix of all mutation types}
}   
Default is "split"}

\item{num_cores}{(Optional) Number of cores used for parallel processing. If no value
is given, then the number of available cores is autodetected.}

\item{indel}{(Optional) List of mutation matrix and vectors for context, class and color of indels.
Color vector must have the same length as the class vector, 
because in plotting profiles, each class is represented by one color.\cr\cr
It is also possible to give a character for predefined variables
and counting mutations with predefined context:
\itemize{
  \item{"predefined"} {Represents a indel context of 3 classes per 
  deletion and indel: "repetitive region", "microhomology" and 
  "none" of these two. Indels have lengths 1 to 5+}
  \item{"cosmic"} {Represents the indel context according to the
  COSMIC database}
}
Default is "cosmic"}

\item{...}{Arguments passed to type_context()}
}
\value{
List with mutation count matrices of snv (96 mutations), dbs (78 
mutations) and indels (number of mutations depends on chosen context for indels)
}
\description{
Make mutation count matrices for 96 trinucleotide single base 
substitutions, 78 double base substitutions and indels. Number of indels 
depends on the indel context given by the user
}
\examples{
## See the 'read_vcfs_as_granges()' example for how we obtained the
## following data:
vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
                package="MutationalPatterns"))

## Load the corresponding reference genome.
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

## Construct a mutation matrix from the loaded VCFs in comparison to the
## ref_genome.
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome, type)

}
\seealso{
\code{\link{read_vcfs_as_granges}},
\code{\link{type_context}}
}
