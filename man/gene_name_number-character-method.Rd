% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pathway.r
\name{gene_name_number,character-method}
\alias{gene_name_number,character-method}
\title{Function to get genes names and numbers from kegg (for internal use)}
\usage{
\S4method{gene_name_number}{character}(x)
}
\arguments{
\item{x}{A \code{character} hsa identifier of the pathway for which gene 
infomation should be extracted as used in KEGG database ('hsa00100').}
}
\value{
A \code{data.frame} listing the genes included in the pathway with 
their names as well as numbers used in KEGG database.
}
\description{
This function extracts for a particular \code{\link{pathway}} all genes and 
the numbers they are represented with in the KEGG network from the 
corresponding KGML pathway file.
}
\author{
Stefanie Friedrichs, Patricia Burger
}
\keyword{internal}
