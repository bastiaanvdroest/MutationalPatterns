% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mutations_from_vcf.R
\name{mutations_from_vcf}
\alias{mutations_from_vcf}
\title{Retrieve mutations from vcf}
\usage{
mutations_from_vcf(vcf, type)
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
A function to extract mutations of each position in vcf
}
\examples{
## See the 'read_vcfs_as_granges()' example for how we obtained the
## following data:
vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
                package="MutationalPatterns"))

muts = mutations_from_vcf(vcfs[[1]])

}
\seealso{
\code{\link{read_vcfs_as_granges}}
}
