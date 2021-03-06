% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/default_colors_ggplot.R
\name{default_colors_ggplot}
\alias{default_colors_ggplot}
\title{Default colors from ggplot2 package}
\usage{
default_colors_ggplot(n)
}
\arguments{
\item{n}{Numeric stating the number of default colors wanted}
}
\value{
Character vector with colornames in HEX format
}
\description{
Get the default colors in HEX format which are used by the ggplot2 
package, such that the colors can be used for the same signatures, 
even if the order or number of signatures differs. \cr\cr
Colors are chosen by equal distribution of the range [0,360]
over the number of colors wanted. 0 yields red, 120 yields green,
240 yields blue etc.
}
\examples{
## Get the default colors of ggplot2 for 10 classes
n = 10

default_colors = default_colors_ggplot(n)

}
