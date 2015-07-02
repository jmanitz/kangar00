\name{pathway-class}
\Rdversion{1.1}
\docType{class}
\alias{pathway-class}
\alias{pathway}
\alias{plot,pathway-method}
\alias{show,pathway-method}
\alias{summary,pathway-method}

\title{Class \code{pathway}}
\description{
This class describes a pathway using an id number, a gene interaction matrix 
a vector of signs differentiation between activating and inhibiting interactions. 
}
\section{Objects from the Class}{
New objects can be created by calls of the form \code{pathway(id, adj, sign)}.
}
\section{Slots}{
  \describe{
    \item{\code{id}:}{Object of class \code{"character"} identifier for a 
    pathway, e.g. hsa00100 as used in the KEGG database.}
    \item{\code{adj}:}{Adjacency matrix of gene interactions}
    \item{\code{sign}:}{Object of class \code{"vector"} containing interaction 
    type for each link.}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(object = "pathway")}: ... }
    \item{show}{\code{signature(object = "pathway")}: Summarize pathway information and prints  pathway adjacency matrix.}
    \item{summary}{\code{signature(object = "pathway")}: Information about nodes, links, activations and inhibitions.}
	 }
}
\author{
Juliane Manitz, Stefanie Friedrichs, Patricia Burger, Benjamin Hofner
}
\seealso{
\code{\link{pathway}}
}
\examples{
showClass("pathway")
}
\keyword{classes}
\keyword{pathway class}
\keyword{pathway methods}
\keyword{pathway slots}

