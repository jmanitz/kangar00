\name{getNetworkmatrix}
\alias{getNetworkmatrix}
\title{
getNetworkmatrix
}
\description{
Calculates the networkmatrix for a pathway specified with the "hsa" argument, using information from the KEGG database. 

}
\usage{
getNetworkmatrix(hsa, directed,  ...)
}

\arguments{
  \item{hsa}{
  pathway identifier as in KEGG database, character.
}
 \item{directed}{
  has to be TRUE or FALSES, stating if the networkmatrix should be returned directed (TRUE) or undirected (FALSE).
}
}

\value{
Returns a matrix representing the interactions between genes in the pathway specified with "hsa".  
}

\author{
Stefanie Friedrichs
}

\keyword{getNetworkmatrix}

