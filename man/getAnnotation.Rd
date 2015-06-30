\name{getAnnotation}
\alias{getAnnotation}
\title{
getAnnotation
}
\description{
Annotates a snp file specifieng marker positions to a list of genes contained in one or more pathways via the stated positions. Uses the package sqldf. 

}
\usage{
getAnnotation(snps, pathwayanno, ...) {
}

\arguments{
  \item{snps}{
  data frame with SNPs (map file), columns called "chr", "position" and "rsnumber" are assumed (format as returned by "updateSNP").
}
 \item{pathwayanno}{
  dataframe, pathway-gene information file as returned from getPathway.
}
}

\value{
Returns a data frame mapping each SNP to the genes in which it is contained. SNPs not falling into genes are not listed. 
}

\author{
Stefanie Friedrichs
}

\keyword{getAnnotation}

