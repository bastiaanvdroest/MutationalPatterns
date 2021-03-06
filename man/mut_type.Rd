% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mut_type.R
\name{mut_type}
\alias{mut_type}
\title{Retrieve mutation types from a VCF object}
\usage{
mut_type(vcf, type)
}
\arguments{
\item{vcf}{A CollapsedVCF object}

\item{type}{(Optional) A character vector stating which type of mutation is to be extracted: 
'snv', 'dbs' and/or 'indel'. All mutation types can also be chosen by 'type = all'.\cr
Default is 'snv'}
}
\value{
List with character vector for each mutation type
}
\description{
A function to extract the mutations from a vcf and translate to
the 6 common base substitution types for SNV and to the 78 strand-agnostic
mutation types of the DBS from COSMIC.
}
\examples{
## See the 'read_vcfs_as_granges()' example for how we obtained the
## following data:
vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
                package="MutationalPatterns"))

mut_type(vcfs[[1]], type)

}
\seealso{
\code{\link{read_vcfs_as_granges}}
}
