\name{getPathway}
\alias{getPathway}
\title{
getPathway
}
\description{
Extracts the genes belonging to the pathway specified with "hsa" from the KEGG database and lits them with their start and end positions. Gene positions are extracted from the Ensemble database, which is accessed using the package BiomaRt.
}
\usage{
getPathway(hsa, ...)
}

\arguments{
  \item{hsa}{
  pathway identifier as in KEGG database, character.
}
}

\value{
Returns a data frame listing the genes included in the pathway specified with the "hsa" argument. The dataframe contains the columns "pathway","GeneLowPoint","GeneHighPoint","Chr" and "Gene". 
}

\author{
Stefanie Friedrichs
}

\keyword{getPathway}

