\name{pathway}
\alias{pathway}
\title{
New pathway object
}
\description{
Creates a new object of the class pathway.
}
\usage{
pathway(id, adj, sign)
}
\arguments{
  \item{id}{
  Id number defining the pathway, e.g. hsa00100 as used in the KEGG database.   
  }
  \item{adj}{
  quadratic adjacency matrix of dimension equaling the number of genes in the 
  pathway namesd in 'id'. Represents the interaction network within the pathway.   
  }
  \item{sign}{
  interaction type for each link. States, whether an interaction is of 
  activation or inhibiting kind.  
  }
}
\author{
Juliane Manitz, Stefanie Friedrichs, Patricia Burger, Benjamin Hofner
}
\seealso{
\code{\link{class}}
}

\keyword{pathway}

