\name{pathway}
\alias{pathway}
\title{
New pathway object
}
\description{
Creates a new object of the class pathway.
}
\usage{
pathway(...)
}
\arguments{
  \item{\code{id}:}{Character whicht specifies ID of pathway e.g. hsa00100.}
  \item{\code{adj}:}{Adjacency matrix}
  \item{\code{sign}:}{Vector containing interaction type for each link.}
}
\seealso{
\code{\link{class}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
structure(function (...) 
new("pathway", ...), className = structure("pathway", package = "kangar00"), package = "kangar00", class = structure("classGeneratorFunction", package = "methods"))
}
\keyword{pathway}

