\name{SNPupdate}
\alias{SNPupdate}
\title{
SNP update
}
\description{
Gives the current position of SNP as listed in the Ensemble database. SNPs are identified via rs-numbers. Ensemble database is accessed using the package BiomaRt.
}
\usage{
SNPupdate(map, column , ...) 
}
\arguments{
  \item{map}{\code{data frame} with including SNPs (map file).}
  \item{column}{numeric, column in data frame in which rsnumbers are specified.}
}
\value{
Returns a data frame with updated SNP positions, containing the columns "chr","position" and "rsnumber". 
}

\author{
Stefanie Friedrichs
}

\keyword{SNPupdate}

